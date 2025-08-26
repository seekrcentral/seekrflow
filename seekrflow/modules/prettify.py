"""
modules/prettify.py

We can print prettier outputs for job information, etc.
"""

import sys
import shutil
import datetime

def C(
        s: str, 
        code: str, 
        use_color: bool
        ) -> str:
    if not use_color: return s
    codes = {
        "b": "\033[1m",   # bold
        "dim": "\033[2m",
        "red": "\033[31m",
        "yel": "\033[33m",
        "grn": "\033[32m",
        "rst": "\033[0m",
    }
    return f"{codes.get(code,'')}{s}{codes['rst'] if code in codes else ''}"

def prettify_job_info(
        result: dict,
        color: bool = True
        ) -> None:
    """
    Prettify the job information for display.
    """
    use_color = sys.stdout.isatty() if color is None else bool(color)
    rows = []
    max_rows = None
    for r in result["rows"]:
        rows.append({
            "JobID": str(r.get("JobID") or r.get("JobId") or r.get("Job") or ""),
            "Name": str(r.get("Name") or r.get("JobName") or ""),
            "User": str(r.get("User") or r.get("user") or ""),
            "State": str(r.get("State") or r.get("JobState") or ""),
            "Partition": str(r.get("Partition") or r.get("partition") or ""),
            "Elapsed": str(r.get("Elapsed") or r.get("time") or ""),
            "TimeLimit": str(r.get("TimeLimit") or r.get("Timelimit") or ""),
            "Nodes": str(r.get("Nodes") or r.get("NNodes") or ""),
            "CPUs": str(r.get("CPUs") or r.get("NCPUS") or ""),
            "ExitCode": str(r.get("ExitCode") or ""),
            "Reason": str(r.get("Reason") or r.get("reason") or ""),
            "Start": str(r.get("Start") or r.get("StartTime") or ""),
            "End": str(r.get("End") or r.get("EndTime") or ""),
        })
    if max_rows is not None and len(rows) > max_rows:
        rows = rows[:max_rows]

    # Decide which columns to show (hide empty ones, keep a sensible order)
    preferred = [
        "JobID","Name","User","State","Partition",
        "Elapsed","TimeLimit","Nodes","CPUs","ExitCode","Reason","Start","End"
    ]
    columns = [c for c in preferred if any(r.get(c) for r in rows)]
    if not columns:
        print("No data to display.")

    # Compute widths and fit to terminal
    term_w = shutil.get_terminal_size(fallback=(120, 24)).columns
    maxw = {
        "JobID": 12, "Name": 28, "User": 14, "State": 12, "Partition": 14,
        "Elapsed": 9, "TimeLimit": 9, "Nodes": 5, "CPUs": 5, "ExitCode": 9,
        "Reason": 45, "Start": 19, "End": 19
    }
    minw = {
        "JobID": 6, "Name": 10, "User": 6, "State": 5, "Partition": 7,
        "Elapsed": 7, "TimeLimit": 7, "Nodes": 3, "CPUs": 3, "ExitCode": 6,
        "Reason": 10, "Start": 16, "End": 16
    }
    widths = {}
    for c in columns:
        raw = max(len(c), *(len(str(r.get(c,""))) for r in rows))
        widths[c] = max(minw.get(c, 4), min(maxw.get(c, raw), raw))

    def total_width():
        # 3 spaces between columns
        return sum(widths[c] for c in columns) + 3*(len(columns)-1)

    def shrink_to(target):
        # Prefer shrinking wide text columns first
        order = [x for x in ["Reason","Name","Partition","User","Start","End"] \
                    if x in columns]
        while total_width() > target:
            changed = False
            for c in order + columns:
                if widths[c] > minw.get(c, 4):
                    widths[c] -= 1
                    changed = True
                    break
            if not changed:
                break

    if total_width() > term_w:
        shrink_to(term_w)

    def trunc(s, w):
        s = "" if s is None else str(s)
        if len(s) <= w:
            return s.ljust(w)
        if w <= 1:
            return "…"
        if w == 2:
            return s[0] + "…"
        return s[:w-1] + "…"

    # Header
    ts = result.get("ts")
    when = ""
    if ts is not None:
        try:
            when = datetime.datetime.fromtimestamp(float(ts)).strftime(
                "%Y-%m-%d %H:%M:%S")
        except Exception:
            pass
    title = (result.get("tool","slurm") or "slurm").upper()
    print(C(f"{title}" + (f" @ {when}" if when else ""), "b", use_color))

    # Table
    header = "   ".join(trunc(c, widths[c]) for c in columns)
    sep    = "   ".join("-"*widths[c] for c in columns)
    print(C(header, "b", use_color))
    print(sep)
    for r in rows:
        line = "   ".join(trunc(r.get(c,""), widths[c]) for c in columns)
        print(line)

    # Footer
    print(C(f"{len(rows)} rows", "dim", use_color))
    return