"""Tests for Globus Compute remote interface helpers."""

from seekrflow.modules.remote_interfaces import globus_compute_sdk as gcs


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


def test_globus_result_timeout_constant():
    assert gcs.GLOBUS_RESULT_TIMEOUT_S == 1200.0
