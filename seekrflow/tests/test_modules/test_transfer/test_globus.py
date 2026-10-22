"""
Tests for seekrflow/modules/transfer/globus.py
"""

from __future__ import annotations

import pytest

import seekrflow.modules.transfer.base as transfer_base
import seekrflow.modules.transfer.globus as transfer_globus
from seekrflow.modules.transfer.aws_s3 import LOCAL_RUNNER_DIRNAMES


class FakeTransferData:
    def __init__(self, source, destination, **kwargs):
        self.source = source
        self.destination = destination
        self.kwargs = kwargs
        self.filters = []

    def add_filter_rule(self, **kwargs):
        self.filters.append(kwargs)

    def add_item(self, *args, **kwargs):
        return None


class FakeSdkError(Exception):
    """Duck-typed globus_sdk-like error."""

    def __init__(
            self,
            *args,
            http_status=None,
            code=None,
            message=None,
            method=None,
            url=None,
            ):
        super().__init__(*args)
        self.http_status = http_status
        self.code = code
        self.message = message
        self.method = method
        self.url = url


FakeSdkError.__module__ = "globus_sdk.exc"


def test_build_transfer_data_excludes_volatile_files():
    tdata = transfer_globus.build_transfer_data(
        "src-id", "dst-id", "label", TransferData=FakeTransferData)
    assert tdata.kwargs["skip_source_errors"] is True
    assert tdata.kwargs["verify_checksum"] is True
    assert tdata.kwargs["encrypt_data"] is True
    assert all(rule["method"] == "exclude" for rule in tdata.filters)
    assert [rule["name"] for rule in tdata.filters] == [
        transfer_globus.STATUS_FILE_NAME,
        transfer_globus.STATUS_TMP_GLOB,
        *LOCAL_RUNNER_DIRNAMES,
    ]


def test_classify_task_document_retries_missing_files_and_fails_closed():
    assert transfer_globus.classify_task_document(
        {"status": "FAILED", "nice_status": "FILE_NOT_FOUND"})[0] == "retry"
    assert transfer_globus.classify_task_document(
        {"status": "ACTIVE", "nice_status": "FILE_NOT_FOUND"})[0] == "retry"
    assert transfer_globus.classify_task_document(
        {"status": "FAILED", "nice_status": "PERMISSION_DENIED"})[0] == "fatal"
    assert transfer_globus.classify_task_document(
        {"status": "INACTIVE"})[0] == "fatal"
    assert transfer_globus.classify_task_document(
        {"status": "ACTIVE", "is_paused": True})[0] == "fatal"


def test_http_zero_formats_as_retryable():
    attributed = FakeSdkError(0, http_status=0)
    bare = Exception(0)
    for error in (attributed, bare):
        formatted = transfer_globus.format_transfer_exception(error)
        assert formatted != "0"
        assert "HTTP 0" in formatted
        wrapped = transfer_globus.wrap_transfer_exception(error)
        assert isinstance(wrapped, transfer_globus.GlobusTransferRetryableError)
        assert transfer_base.is_retryable_transfer_error(error)
    assert transfer_base.is_retryable_transfer_error(
        transfer_globus.GlobusTransferError("globus-transfer: HTTP 0"))


def test_post_tuple_and_rate_limit_are_retryable_auth_failure_is_not():
    url = "https://transfer.api.globus.org/v0.10/transfer"
    posted = FakeSdkError("POST", url)
    formatted = transfer_globus.format_transfer_exception(posted)
    assert formatted.startswith("globus-transfer:")
    assert "POST" in formatted
    assert "transfer.api.globus.org" in formatted
    assert isinstance(
        transfer_globus.wrap_transfer_exception(posted),
        transfer_globus.GlobusTransferRetryableError)

    too_many = FakeSdkError(
        http_status=429, code="RateLimitExceeded", method="POST", url=url,
        message="Too Many Requests")
    wrapped_429 = transfer_globus.wrap_transfer_exception(too_many)
    assert isinstance(wrapped_429, transfer_globus.GlobusTransferRetryableError)
    assert "HTTP 429" in str(wrapped_429)
    assert "RateLimitExceeded" in str(wrapped_429)

    denied = FakeSdkError(
        http_status=401, code="AUTHENTICATION_FAILED", method="POST", url=url)
    wrapped_401 = transfer_globus.wrap_transfer_exception(denied)
    assert isinstance(wrapped_401, transfer_globus.GlobusTransferError)
    assert not isinstance(
        wrapped_401, transfer_globus.GlobusTransferRetryableError)
    assert not transfer_base.is_retryable_transfer_error(wrapped_401)

    checksum = FakeSdkError("Checksum mismatch on status file")
    assert isinstance(
        transfer_globus.wrap_transfer_exception(checksum),
        transfer_globus.GlobusTransferRetryableError)


def _patch_wait(monkeypatch):
    monkeypatch.setattr(transfer_globus, "GLOBUS_TRANSFER_MAX_ERROR_POLLS", 1)
    monkeypatch.setattr(transfer_globus, "GLOBUS_TRANSFER_POLL_INTERVAL_S", 0.001)


def _transfer(local, monkeypatch, client, **kwargs):
    monkeypatch.setattr(
        transfer_globus, "get_transfer_client", lambda: client)
    transfer_globus.transfer_files_with_globus(
        "sys", str(local), "/remote", "local-id", "remote-id",
        retry_sleep=0.0, **kwargs)


def test_skips_submit_when_destination_already_has_model(tmp_path, monkeypatch):
    local = tmp_path / "local"
    local.mkdir()
    (local / "model.json").write_text("{}")

    class Client:
        def __init__(self):
            self.submits = []

        def add_app_data_access_scope(self, collection_id):
            return None

        def operation_ls(self, collection_id, path=None):
            return [{"name": "model.json"}]

        def submit_transfer(self, tdata):
            self.submits.append(tdata)
            return {"task_id": "should-not-happen"}

    client = Client()
    _transfer(local, monkeypatch, client)
    assert client.submits == []


def test_http0_after_submit_succeeds_when_model_arrives(tmp_path, monkeypatch):
    _patch_wait(monkeypatch)
    local = tmp_path / "local"
    local.mkdir()
    (local / "model.json").write_text("{}")

    class Client:
        def __init__(self):
            self.submits = []
            self.cancelled = []
            self._listings = 0

        def add_app_data_access_scope(self, collection_id):
            return None

        def operation_ls(self, collection_id, path=None):
            self._listings += 1
            if self._listings == 1:
                return []
            return [{"name": "model.json"}]

        def submit_transfer(self, tdata):
            self.submits.append(1)
            return {"task_id": "t1"}

        def task_wait(self, task_id, timeout=0, polling_interval=0):
            return True

        def get_task(self, task_id):
            raise Exception(0)

        def cancel_task(self, task_id):
            self.cancelled.append(task_id)

    client = Client()
    _transfer(local, monkeypatch, client, max_attempts=2)
    assert client.submits == [1]
    assert client.cancelled == []


def test_http0_after_submit_does_not_submit_again(tmp_path, monkeypatch):
    _patch_wait(monkeypatch)
    local = tmp_path / "local"
    local.mkdir()
    (local / "model.json").write_text("{}")

    class Client:
        def __init__(self):
            self.submits = []

        def add_app_data_access_scope(self, collection_id):
            return None

        def operation_ls(self, collection_id, path=None):
            return []

        def submit_transfer(self, tdata):
            self.submits.append(1)
            return {"task_id": "t1"}

        def task_wait(self, task_id, timeout=0, polling_interval=0):
            return True

        def get_task(self, task_id):
            raise Exception(0)

        def cancel_task(self, task_id):
            return None

    client = Client()
    with pytest.raises(transfer_globus.GlobusTransferRetryableError):
        _transfer(local, monkeypatch, client, max_attempts=3)
    assert client.submits == [1]


def test_submit_http0_retries_then_raises_classified_error(tmp_path, monkeypatch):
    local = tmp_path / "empty"
    local.mkdir()

    class Client:
        def __init__(self):
            self.submits = []

        def add_app_data_access_scope(self, collection_id):
            return None

        def submit_transfer(self, tdata):
            self.submits.append(1)
            raise FakeSdkError(0, http_status=0)

    client = Client()
    with pytest.raises(
            transfer_globus.GlobusTransferRetryableError, match="HTTP 0"):
        _transfer(local, monkeypatch, client, max_attempts=2)
    assert client.submits == [1, 1]


def test_repeated_file_not_found_cancels_the_task(tmp_path, monkeypatch):
    _patch_wait(monkeypatch)
    local = tmp_path / "local"
    local.mkdir()
    (local / "model.json").write_text("{}")

    class Client:
        def __init__(self):
            self.cancelled = []

        def add_app_data_access_scope(self, collection_id):
            return None

        def operation_ls(self, collection_id, path=None):
            return []

        def submit_transfer(self, tdata):
            return {"task_id": "t1"}

        def task_wait(self, task_id, timeout=0, polling_interval=0):
            return False

        def get_task(self, task_id):
            return {"status": "ACTIVE", "nice_status": "FILE_NOT_FOUND"}

        def cancel_task(self, task_id):
            self.cancelled.append(task_id)

    client = Client()
    with pytest.raises(transfer_globus.GlobusTransferRetryableError):
        _transfer(local, monkeypatch, client, max_attempts=1)
    assert client.cancelled == ["t1"]
