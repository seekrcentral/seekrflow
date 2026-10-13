"""
modules/transfer/globus.py

Handle file transfers via Globus, which can be useful for file transfers
where password entry and two-factor authentication would otherwise be
burdensome.
"""

# TODO: see if this module can be simplified

from __future__ import annotations

import os
import time

GLOBUS_SEEKRFLOW_APP_CLIENT_ID = "683ab038-1578-4520-bfb4-57de7411102f"
GLOBUS_TRANSFER_RESOURCE_SERVER = "transfer.api.globus.org"

# How long to block per status poll, and the total budget for one attempt.
GLOBUS_TRANSFER_POLL_INTERVAL_S = 30.0
GLOBUS_TRANSFER_TIMEOUT_S = 3600.0

# Whole-transfer retries for transient conditions (busy/unreachable endpoint).
# Also the cap on queued relaunches when the copy cannot be confirmed.
GLOBUS_TRANSFER_MAX_ATTEMPTS = 3
GLOBUS_TRANSFER_RETRY_SLEEP_S = 30.0

# How many consecutive polls a task may report a transient error before we
# give up on this attempt instead of waiting out the full timeout.
GLOBUS_TRANSFER_MAX_ERROR_POLLS = 5

# ``nice_status`` values that will not resolve on their own.
# FILE_NOT_FOUND is omitted on purpose: a recursive copy of a live model tree
# routinely lists a status tempfile that is gone by the time Globus reads it.
# A missing model.json is caught by verify_destination_contents instead.
FATAL_NICE_STATUSES = frozenset({
    "PERMISSION_DENIED",
    "PATH_NOT_FOUND",
    "NO_CREDENTIALS",
    "EXPIRED_CREDENTIALS",
    "AUTHENTICATION_FAILED",
    "CONSENT_REQUIRED",
    "QUOTA_EXCEEDED",
    "ENDPOINT_TOO_BUSY_NO_RETRY",
})

# HTTP statuses that mean "try again", vs. a human has to fix credentials.
RETRYABLE_HTTP_STATUSES = frozenset({0, 408, 425, 429})
FATAL_HTTP_STATUSES = frozenset({401, 403})

# Live files the monitor writes while Globus is hashing the tree. Checksum
# mismatch and vanished temps come from these, not from model.json.
# TODO: are these even still used?
STATUS_FILE_NAME = ".seekrflow_job_status.json"
STATUS_TMP_GLOB = ".seekrflow_job_status.*"

# Files that must exist on the destination after a forward transfer, when they
# exist locally. A silently failed transfer used to look like success here,
# which then showed up much later as a remote "Model file not found" error.
# TODO: is this insufficient? What if the transfer is partially successful, 
# and the model.json file is transferred, but some other key files are not?
FORWARD_TRANSFER_SENTINELS = ("model.json",)

_ACTIVE_TASKS: dict[str, object] = {}
_CANCELLED_TASKS: set[str] = set()

class GlobusTransferError(RuntimeError):
    """A Globus transfer task did not deliver the files."""


class GlobusTransferRetryableError(GlobusTransferError):
    """Transient transfer failure; retrying may succeed."""


def volatile_exclude_globs() -> list[str]:
    """
    Names Globus must not copy. The monitor rewrites these while a transfer
    runs, which is what produced checksum failures and FILE_NOT_FOUND on
    otherwise-complete remote trees.
    """
    from seekrflow.modules.transfer.aws_s3 import LOCAL_RUNNER_DIRNAMES

    return [
        STATUS_FILE_NAME,
        STATUS_TMP_GLOB,
        *LOCAL_RUNNER_DIRNAMES,
    ]


def apply_volatile_excludes(tdata) -> None:
    """Attach exclude filter rules to a TransferData object."""
    for name in volatile_exclude_globs():
        tdata.add_filter_rule(method="exclude", name=name)


def build_transfer_data(
        source_collection_id: str,
        destination_collection_id: str,
        label: str,
        TransferData=None,
        ):
    """
    TransferData for a live model tree: checksum the payload, skip files that
    vanish mid-listing, and exclude monitor/runner state.
    """
    if TransferData is None:
        from globus_sdk import TransferData as TransferData
    tdata = TransferData(
        source_collection_id,
        destination_collection_id,
        label=label,
        sync_level="mtime",
        verify_checksum=True,
        encrypt_data=True,
        skip_source_errors=True,
    )
    apply_volatile_excludes(tdata)
    return tdata


def http_status_from_error(error: BaseException) -> int | None:
    """
    Infer an HTTP status from a globus_sdk / urllib3 / requests error.

    Real connection failures often have ``args==(0,)`` or ``str(error)=="0"``
    and no ``http_status`` attribute. Formatting used to report HTTP 0 while
    classification treated the same object as fatal.
    """
    status = getattr(error, "http_status", None)
    try:
        if status is not None:
            return int(status)
    except (TypeError, ValueError):
        pass
    args = getattr(error, "args", ()) or ()
    if len(args) == 1:
        if isinstance(args[0], int):
            return args[0]
        if isinstance(args[0], str) and args[0].strip() == "0":
            return 0
    text = str(error).strip()
    if text == "0":
        return 0
    upper = text.upper()
    marker = "HTTP "
    idx = upper.find(marker)
    if idx >= 0:
        digits: list[str] = []
        for ch in text[idx + len(marker):]:
            if ch.isdigit():
                digits.append(ch)
            else:
                break
        if digits:
            try:
                return int("".join(digits))
            except ValueError:
                return None
    return None


def format_transfer_exception(error: BaseException) -> str:
    """
    Readable message for globus_sdk errors, which often stringify as ``0`` or
    a truncated ``('POST', 'https://transfer.api.globus.org/...')`` tuple.
    """
    if isinstance(error, GlobusTransferError):
        return str(error)

    method = getattr(error, "method", None) or getattr(error, "http_method", None)
    url = getattr(error, "url", None)
    http_status = http_status_from_error(error)
    code = getattr(error, "code", None)
    message = getattr(error, "message", None)
    args = getattr(error, "args", ()) or ()

    if method is None and args:
        first = args[0]
        if isinstance(first, str) and first.upper() in {
                "GET", "POST", "PUT", "PATCH", "DELETE"}:
            method = first
            if len(args) > 1 and url is None:
                url = args[1]
        elif isinstance(first, tuple) and first:
            if isinstance(first[0], str):
                method = first[0]
            if len(first) > 1 and url is None:
                url = first[1]

    bits: list[str] = []
    if method or url:
        bits.append(f"{method or '?'} {url or ''}".strip())
    if http_status is not None:
        bits.append(f"HTTP {http_status}")
    if code:
        bits.append(str(code))
    if message:
        bits.append(str(message))
    if bits:
        return "globus-transfer: " + " -> ".join(bits)

    raw = str(error).strip()
    typename = type(error).__name__
    if raw in {"", "0"}:
        status_bit = (
            f"HTTP {http_status}" if http_status is not None else "HTTP 0")
        return f"globus-transfer: {typename} ({status_bit})"
    return f"globus-transfer: {raw}"


def sdk_error_is_retryable(error: BaseException) -> bool:
    """
    True for network blips, HTTP 0/429/5xx, and checksum mismatches.

    False for 401/403 and consent/credential failures that need a human.
    """
    if isinstance(error, GlobusTransferRetryableError):
        return True
    if isinstance(error, (TimeoutError, ConnectionError, BrokenPipeError)):
        return True
    if isinstance(error, OSError) and not isinstance(error, FileNotFoundError):
        return True

    http_status = http_status_from_error(error)
    if http_status in FATAL_HTTP_STATUSES:
        return False
    if http_status is not None and (
            http_status in RETRYABLE_HTTP_STATUSES or http_status >= 500):
        return True
    if isinstance(error, GlobusTransferError):
        return False

    code = str(getattr(error, "code", "") or "")
    text = f"{code} {error}".upper()
    if any(token in text for token in (
            "CONSENT_REQUIRED", "AUTHENTICATION_FAILED", "PERMISSION_DENIED",
            "EXPIRED_CREDENTIALS", "NO_CREDENTIALS")):
        return False
    if "CHECKSUM" in text:
        return True
    module = getattr(type(error), "__module__", "") or ""
    if "globus" in module.lower():
        return True
    if "transfer.api.globus" in str(error).lower():
        return True
    if http_status is not None:
        return http_status not in FATAL_HTTP_STATUSES | {400, 404}
    return False


def wrap_transfer_exception(error: BaseException) -> GlobusTransferError:
    """Classify a raw SDK/network error as fatal or retryable."""
    if isinstance(error, GlobusTransferRetryableError):
        return error
    if isinstance(error, GlobusTransferError):
        if sdk_error_is_retryable(error):
            return GlobusTransferRetryableError(str(error))
        return error
    message = format_transfer_exception(error)
    if sdk_error_is_retryable(error):
        return GlobusTransferRetryableError(message)
    return GlobusTransferError(message)


def get_transfer_client():
    """
    Build an authorized TransferClient using the seekrflow Globus app.
    """
    from globus_sdk import TransferClient, UserApp, GlobusAppConfig

    # The GlobusApp framework handles token storage (in ~/.globus/app/),
    # the login flow, and authorization automatically. When a required
    # token is missing, it prints the login URL and prompts for the code.
    #
    # Force the command-line login flow so it uses the redirect URI that is
    # registered on the native ("Thick Client") app
    # (https://auth.globus.org/v2/web/auth-code). The default would
    # auto-select a local-server flow that redirects to http://localhost:<port>,
    # which is not a registered redirect URI and causes a
    # "Mismatching redirect URI" error.
    config = GlobusAppConfig(
        login_flow_manager="command-line",
        request_refresh_tokens=True,
    )
    app = UserApp(
        "seekrflow",
        client_id=GLOBUS_SEEKRFLOW_APP_CLIENT_ID,
        config=config,
    )
    return TransferClient(app=app)


def classify_task_document(task: dict) -> tuple[str, str]:
    """
    Classify an in-flight or finished transfer task.

    Returns ``(verdict, detail)`` where verdict is one of ``ok`` (progressing
    or succeeded), ``retry`` (transient trouble), or ``fatal``.
    """
    status = str(task.get("status") or "").upper()
    nice_status = str(task.get("nice_status") or "").upper()
    detail = (
        task.get("nice_status_short_description")
        or task.get("nice_status_details")
        or nice_status
        or status
    )
    if status == "SUCCEEDED":
        return "ok", str(detail)
    if status == "FAILED":
        if nice_status and nice_status in FATAL_NICE_STATUSES:
            return "fatal", str(detail)
        return "retry", str(detail)
    if status == "INACTIVE":
        # Credentials or consent expired; a human has to intervene.
        return "fatal", str(detail or "task is INACTIVE")
    if task.get("is_paused"):
        return "fatal", f"task is paused ({detail})"
    if nice_status and nice_status not in {"OK", "QUEUED"}:
        if nice_status in FATAL_NICE_STATUSES:
            return "fatal", str(detail)
        return "retry", str(detail)
    return "ok", str(detail)


def await_transfer_task(
        transfer_client,
        task_id: str,
        label: str,
        timeout: float = GLOBUS_TRANSFER_TIMEOUT_S,
        poll_interval: float = GLOBUS_TRANSFER_POLL_INTERVAL_S,
        ) -> dict:
    """
    Wait for one transfer task and confirm it actually succeeded.

    Unlike a bare ``task_wait`` loop, this bounds the wait, inspects the task
    document while it runs, and raises instead of reporting success for a
    task that failed or is stuck erroring.

    Network errors (HTTP 0) do not cancel the task: the copy may already have
    finished, and cancelling would throw away a successful transfer.
    """
    deadline = time.monotonic() + max(poll_interval, timeout)
    error_polls = 0
    network_polls = 0
    last_task: dict | None = None
    while True:
        try:
            if task_id in _CANCELLED_TASKS:
                raise GlobusTransferError(
                    f"globus-transfer: task {task_id} ({label}) cancelled")
            finished = transfer_client.task_wait(
                task_id, timeout=poll_interval, polling_interval=poll_interval)
            task = dict(transfer_client.get_task(task_id))
        except Exception as e:
            classified = wrap_transfer_exception(e)
            if not isinstance(classified, GlobusTransferRetryableError):
                raise classified from e
            network_polls += 1
            print(
                f"[transfer] globus task {task_id} ({label}) API error "
                f"{classified} ({network_polls}/"
                f"{GLOBUS_TRANSFER_MAX_ERROR_POLLS}); not cancelling")
            if (network_polls >= GLOBUS_TRANSFER_MAX_ERROR_POLLS
                    or time.monotonic() >= deadline):
                raise classified from e
            continue
        last_task = task
        network_polls = 0
        verdict, detail = classify_task_document(task)
        status = str(task.get("status") or "").upper()
        if verdict == "fatal":
            raise GlobusTransferError(
                f"globus-transfer: task {task_id} ({label}) failed: {detail}")
        if finished or status in {"SUCCEEDED", "FAILED"}:
            if status == "SUCCEEDED":
                return task
            raise GlobusTransferRetryableError(
                f"globus-transfer: task {task_id} ({label}) ended in "
                f"{status or 'UNKNOWN'}: {detail}")
        if verdict == "retry":
            error_polls += 1
            print(
                f"[transfer] globus task {task_id} ({label}) reporting "
                f"{detail} ({error_polls}/{GLOBUS_TRANSFER_MAX_ERROR_POLLS})")
            if error_polls >= GLOBUS_TRANSFER_MAX_ERROR_POLLS:
                cancel_transfer_task(transfer_client, task_id)
                raise GlobusTransferRetryableError(
                    f"globus-transfer: task {task_id} ({label}) kept "
                    f"reporting {detail}")
        else:
            error_polls = 0
        if time.monotonic() >= deadline:
            if last_task is not None:
                cancel_transfer_task(transfer_client, task_id)
            raise GlobusTransferRetryableError(
                f"globus-transfer: task {task_id} ({label}) did not finish "
                f"within {timeout:.0f}s (last status {status or 'UNKNOWN'}: "
                f"{detail}); cancelled it")


def cancel_transfer_task(transfer_client, task_id: str) -> None:
    """
    Best-effort cancel so an abandoned task cannot race a retry.
    """
    try:
        transfer_client.cancel_task(task_id)
        print(f"[transfer] cancelled globus task {task_id}")
    except Exception as e:
        print(f"[transfer] could not cancel globus task {task_id}: {e}")

def cancel_running_transfers() -> None:
    """
    Cancel all running transfer tasks.
    """
    for task_id, transfer_client in list(_ACTIVE_TASKS.items()):
        _CANCELLED_TASKS.add(task_id)
        cancel_transfer_task(transfer_client, task_id)

def list_destination_names(
        transfer_client,
        collection_id: str,
        remote_path: str,
        ) -> set:
    """
    Names in a destination directory. Listing failures are retryable.
    """
    try:
        listing = transfer_client.operation_ls(collection_id, path=remote_path)
        return {entry.get("name") for entry in listing}
    except GlobusTransferError:
        raise
    except Exception as e:
        raise wrap_transfer_exception(e) from e


def destination_sentinels_present(
        transfer_client,
        collection_id: str,
        remote_path: str,
        required_names: list[str],
        ) -> bool:
    """True when every required name is on the destination."""
    if not required_names:
        return False
    present = list_destination_names(
        transfer_client, collection_id, remote_path)
    return all(name in present for name in required_names)


def verify_destination_contents(
        transfer_client,
        collection_id: str,
        remote_path: str,
        required_names: list[str],
        ) -> None:
    """
    Confirm required files landed on the destination collection.
    """
    if not required_names:
        return
    present = list_destination_names(
        transfer_client, collection_id, remote_path)
    missing = [name for name in required_names if name not in present]
    if missing:
        raise GlobusTransferRetryableError(
            f"globus-transfer: {', '.join(missing)} missing from "
            f"{remote_path} after transfer reported success")


def transfer_files_with_globus(
        name: str,
        local_path: str,
        remote_path: str,
        local_collection_id: str,
        remote_collection_id: str,
        backwards: bool = False,
        max_attempts: int = GLOBUS_TRANSFER_MAX_ATTEMPTS,
        retry_sleep: float = GLOBUS_TRANSFER_RETRY_SLEEP_S,
        ) -> None:
    """
    Transfer files to or from a remote system using Globus.

    Raises ``GlobusTransferError`` if the files did not arrive; callers rely on
    that to avoid submitting jobs against an unpopulated remote directory.

    HTTP 0 / network errors after a task is submitted retry status and listing
    only — they never start a second copy. If the destination already has the
    sentinel files, this returns success without submitting.
    """
    # Ensure trailing slashes for directories
    if not local_path.endswith("/"):
        local_path += "/"
    if not remote_path.endswith("/"):
        remote_path += "/"

    transfer_client = get_transfer_client()
    transfer_client.add_app_data_access_scope(remote_collection_id)

    if backwards:
        # Transferring from remote to local
        source_path = remote_path
        destination_path = local_path
        source_collection_id = remote_collection_id
        destination_collection_id = local_collection_id
        required_names: list[str] = []
    else:
        # Transferring from local to remote
        source_path = local_path
        destination_path = remote_path
        source_collection_id = local_collection_id
        destination_collection_id = remote_collection_id
        required_names = [
            sentinel for sentinel in FORWARD_TRANSFER_SENTINELS
            if os.path.exists(os.path.join(local_path, sentinel))
        ]

    label = f"{name} files transfer"
    last_error: Exception | None = None
    task_id: str | None = None
    max_attempts = max(1, max_attempts)

    def _destination_ready() -> bool:
        if not required_names:
            return False
        return destination_sentinels_present(
            transfer_client, destination_collection_id, destination_path,
            required_names)

    for attempt in range(1, max_attempts + 1):
        try:
            if _destination_ready():
                print(
                    f"[transfer] {label}: destination already has "
                    f"{', '.join(required_names)}; skipping copy")
                return
        except GlobusTransferRetryableError as e:
            last_error = e
            print(
                f"[transfer] could not list destination "
                f"({e}); {'polling existing task' if task_id else 'will submit'}")

        if task_id is None:
            tdata = build_transfer_data(
                source_collection_id, destination_collection_id, label)
            tdata.add_item(source_path, destination_path, recursive=True)
            try:
                submit_result = transfer_client.submit_transfer(tdata)
                task_id = submit_result["task_id"]
                print(
                    f"[transfer] globus task {task_id} submitted "
                    f"({source_path} -> {destination_path})")
            except GlobusTransferError as e:
                last_error = e
                if isinstance(e, GlobusTransferRetryableError) \
                        and attempt < max_attempts:
                    print(
                        f"[transfer] submit attempt {attempt}/{max_attempts} "
                        f"failed ({e}); retrying in {retry_sleep:.0f}s")
                    time.sleep(retry_sleep)
                    continue
                raise
            except Exception as e:
                classified = wrap_transfer_exception(e)
                last_error = classified
                if isinstance(classified, GlobusTransferRetryableError) \
                        and attempt < max_attempts:
                    print(
                        f"[transfer] submit attempt {attempt}/{max_attempts} "
                        f"raised ({classified}); retrying in {retry_sleep:.0f}s")
                    time.sleep(retry_sleep)
                    continue
                raise classified from e

        _ACTIVE_TASKS[task_id] = transfer_client
        try:
            await_transfer_task(transfer_client, task_id, label)
            verify_destination_contents(
                transfer_client, destination_collection_id, destination_path,
                required_names)
            print(f"[transfer] globus task {task_id} complete")
            return
        except GlobusTransferRetryableError as e:
            last_error = e
            try:
                if _destination_ready():
                    print(
                        f"[transfer] globus task {task_id} poll failed "
                        f"({e}) but destination has required files; "
                        "treating as complete")
                    return
            except GlobusTransferRetryableError:
                pass
            if attempt < max_attempts:
                print(
                    f"[transfer] attempt {attempt}/{max_attempts} failed "
                    f"({e}); retrying status/ls in {retry_sleep:.0f}s "
                    "(not submitting again)")
                time.sleep(retry_sleep)
                continue
            raise
        except GlobusTransferError:
            raise
        except Exception as e:
            classified = wrap_transfer_exception(e)
            last_error = classified
            try:
                if _destination_ready():
                    print(
                        f"[transfer] globus task {task_id} raised "
                        f"({classified}) but destination has required files; "
                        "treating as complete")
                    return
            except GlobusTransferRetryableError:
                pass
            if isinstance(classified, GlobusTransferRetryableError) \
                    and attempt < max_attempts:
                print(
                    f"[transfer] attempt {attempt}/{max_attempts} raised "
                    f"({classified}); retrying status/ls in {retry_sleep:.0f}s "
                    "(not submitting again)")
                time.sleep(retry_sleep)
                continue
            raise classified from e
        finally:
            _ACTIVE_TASKS.pop(task_id, None)

    if last_error is not None:
        raise last_error
