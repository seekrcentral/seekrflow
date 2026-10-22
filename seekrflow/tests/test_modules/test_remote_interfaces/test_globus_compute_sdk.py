"""
Tests for seekrflow/modules/remote_interfaces/globus_compute_sdk.py
"""

import sys
import types

import pytest

from seekrflow.modules.remote_interfaces import globus_compute_sdk as gcs


def _workload(payload):
    return payload


class _HttpError(Exception):
    def __init__(self, status):
        super().__init__(f"http {status}")
        self.http_status = status


def _reset_globus_state():
    gcs._last_endpoint_warn.clear()
    gcs._function_ids.clear()
    gcs._client = None


def _install_sdk(monkeypatch, status, outcomes):
    _reset_globus_state()
    state = {
        "registers": [],
        "submits": [],
        "shutdowns": [],
        "timeouts": [],
        "outcomes": list(outcomes),
    }

    class FakeClient:
        def get_endpoint_status(self, endpoint):
            return status

    class FakeFuture:
        def result(self, timeout=None):
            state["timeouts"].append(timeout)
            outcome = state["outcomes"].pop(0)
            if isinstance(outcome, BaseException):
                raise outcome
            return outcome

    class FakeExecutor:
        def __init__(self, endpoint, client=None):
            self.endpoint = endpoint

        def register_function(self, workload, description="run"):
            function_id = f"fid-{self.endpoint}-{len(state['registers'])}"
            state["registers"].append(self.endpoint)
            return function_id

        def submit_to_registered_function(self, function_id=None, args=None):
            state["submits"].append((function_id, args))
            return FakeFuture()

        def shutdown(self, wait=False, cancel_futures=False):
            state["shutdowns"].append((wait, cancel_futures))

    fake_mod = types.ModuleType("globus_compute_sdk")
    fake_mod.Client = FakeClient
    fake_mod.Executor = FakeExecutor
    fake_ser = types.ModuleType("globus_compute_sdk.serialize")

    class ComputeSerializer:
        def __init__(self, strategy_code=None):
            pass

    class CombinedCode:
        pass

    fake_ser.ComputeSerializer = ComputeSerializer
    fake_ser.CombinedCode = CombinedCode
    monkeypatch.setitem(sys.modules, "globus_compute_sdk", fake_mod)
    monkeypatch.setitem(sys.modules, "globus_compute_sdk.serialize", fake_ser)
    return state


def _online(idle_workers):
    return {
        "status": "online",
        "details": {"idle_workers": idle_workers, "total_workers": 2},
    }


def test_warn_endpoint_dedupes_same_key(capsys):
    gcs._last_endpoint_warn.clear()
    name = "ep-dedupe-test"
    gcs._warn_endpoint(name, "online_idle0", "WARN A")
    gcs._warn_endpoint(name, "online_idle0", "WARN B")
    out = capsys.readouterr().out
    assert out.count("WARN A") == 1
    assert "WARN B" not in out


def test_warn_endpoint_prints_on_key_change(capsys):
    gcs._last_endpoint_warn.clear()
    name = "ep-key-change"
    gcs._warn_endpoint(name, "online_idle0", "WARN IDLE")
    gcs._warn_endpoint(name, "state_offline", "WARN OFF")
    out = capsys.readouterr().out
    assert "WARN IDLE" in out
    assert "WARN OFF" in out


def test_warn_endpoint_reprints_after_interval(capsys, monkeypatch):
    gcs._last_endpoint_warn.clear()
    name = "ep-interval"
    clock = {"t": 1000.0}
    monkeypatch.setattr(gcs.time, "monotonic", lambda: clock["t"])
    monkeypatch.setattr(gcs, "_ENDPOINT_WARN_INTERVAL_S", 10.0)

    gcs._warn_endpoint(name, "online_idle0", "WARN 1")
    clock["t"] = 1005.0
    gcs._warn_endpoint(name, "online_idle0", "WARN 2")
    clock["t"] = 1011.0
    gcs._warn_endpoint(name, "online_idle0", "WARN 3")

    out = capsys.readouterr().out
    assert "WARN 1" in out
    assert "WARN 2" not in out
    assert "WARN 3" in out


def test_submit_reuses_function_id_per_endpoint(monkeypatch):
    payload = {"stage": "mmvt"}
    state = _install_sdk(
        monkeypatch, _online(1), [{"success": True}] * 3)
    assert gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload) == {"success": True}
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload)
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-b", payload)

    assert state["registers"] == ["ep-a", "ep-b"]
    assert state["submits"][0][0] == state["submits"][1][0]
    assert state["submits"][2][0] != state["submits"][0][0]
    assert state["submits"][0][1] == ([payload],)
    assert state["timeouts"] == [gcs.GLOBUS_RESULT_TIMEOUT_S] * 3
    assert state["shutdowns"] == [(False, True)] * 3


def test_submit_wraps_result_timeout(monkeypatch):
    state = _install_sdk(
        monkeypatch, _online(1), [TimeoutError("slow")])
    with pytest.raises(TimeoutError, match="delta"):
        gcs.submit_remote_workload_with_globus_compute(
            "delta", _workload, "ep-a", {"stage": "mmvt"})
    assert state["timeouts"] == [gcs.GLOBUS_RESULT_TIMEOUT_S]
    assert state["shutdowns"] == [(False, True)]


def test_submit_forgets_function_id_only_for_missing_registration(monkeypatch):
    payload = {"stage": "mmvt"}
    state = _install_sdk(
        monkeypatch,
        _online(1),
        [{"ok": True}, _HttpError(404), {"ok": True}],
    )
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload)
    with pytest.raises(_HttpError):
        gcs.submit_remote_workload_with_globus_compute(
            "delta", _workload, "ep-a", payload)
    assert state["registers"] == ["ep-a"]
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload)
    assert state["registers"] == ["ep-a", "ep-a"]
    assert state["shutdowns"] == [(False, True)] * 3

    state = _install_sdk(
        monkeypatch,
        _online(1),
        [{"ok": True}, RuntimeError("boom"), {"ok": True}],
    )
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload)
    with pytest.raises(RuntimeError, match="boom"):
        gcs.submit_remote_workload_with_globus_compute(
            "delta", _workload, "ep-a", payload)
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", payload)
    assert state["registers"] == ["ep-a"]


def test_submit_warns_once_for_idle_endpoint_unless_silent(monkeypatch, capsys):
    state = _install_sdk(
        monkeypatch, _online(0), [{"success": True}] * 3)
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", {"stage": "mmvt"}, silent=True)
    assert capsys.readouterr().out == ""
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", {"stage": "mmvt"})
    gcs.submit_remote_workload_with_globus_compute(
        "delta", _workload, "ep-a", {"stage": "mmvt"})
    assert capsys.readouterr().out.count("0 idle workers") == 1
    assert len(state["shutdowns"]) == 3
