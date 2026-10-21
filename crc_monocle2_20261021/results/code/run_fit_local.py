"""Record and supervise one foreground fit in this checkout's isolated environment."""
import argparse
import datetime
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / "results/04_日志和检查"


def now():
    return datetime.datetime.now().astimezone().isoformat()


def save(path, record):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=["Monocle2", "Slingshot"])
    method = parser.parse_args().method
    LOG.mkdir(parents=True, exist_ok=True)
    with (LOG / "fit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if method == "Slingshot" and not (LOG / "CRC_Monocle2_completed.json").exists():
            raise RuntimeError("Save and review Monocle2 before the method comparison")
        command = [sys.executable, "-u", str(ROOT / "Monocle2_PC.py"), "fit",
                   "--rscript", str(ROOT / ".venv/bin/Rscript"), "--method", method]
        record = {"method": method, "command": command, "cwd": str(ROOT),
                  "started_at": now(), "status": "running", "returncode": None,
                  "formal_fit_started": True, "analysis_accepted": False}
        run_file = LOG / f"{method}_local_process.json"
        if run_file.exists():
            raise RuntimeError("Run record already exists; inspect it before a versioned retry")
        started = time.monotonic()
        with (LOG / f"{method}_launcher.log").open("x") as console:
            child = subprocess.Popen(command, cwd=ROOT, stdout=console,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            record["launcher_child_pid"] = child.pid
            save(run_file, record)
            save(LOG / "progress.json", record)
            with (ROOT / "results/RUN_SUMMARY.md").open("a") as summary:
                summary.write(f"\n{now()}：启动全量{method}，入口PID {child.pid}；"
                              "命令、目录及退出状态保存在日志目录。拟合保存不等于验收通过。\n")
            code = child.wait()
        record.update(returncode=code, finished_at=now(),
                      elapsed_seconds=time.monotonic() - started,
                      status="fit_saved_pending_acceptance" if code == 0 else "fit_failed")
        save(run_file, record)
        save(LOG / "progress.json", record)
        with (ROOT / "results/RUN_SUMMARY.md").open("a") as summary:
            summary.write(f"\n{now()}：{method}入口退出码{code}，"
                          f"用时{record['elapsed_seconds']:.1f}秒；"
                          "具体进度和警告见日志。\n")
        return code


if __name__ == "__main__":
    raise SystemExit(main())
