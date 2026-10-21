"""Observe saved process/stage records and host memory; never controls a fit."""
import datetime
import json
from pathlib import Path
import sys
import time

root = Path(__file__).resolve().parents[2]
logs = root / "results/04_日志和检查"
method = sys.argv[1]
previous = None

with (logs / f"{method}_memory_history.jsonl").open("a", buffering=1) as history:
    while True:
        process = json.loads((logs / f"{method}_local_process.json").read_text())
        memory = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            name = line.split(":", 1)[0]
            if name in {"MemTotal", "MemAvailable", "SwapTotal", "SwapFree"}:
                memory[name + "_kB"] = int(line.split()[1])
        try:
            stage = json.loads((logs / f"CRC_{method}_progress.json").read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            stage = {"stage": "preflight"}
        stamp = datetime.datetime.now().astimezone().isoformat()
        history.write(json.dumps(dict(at=stamp, memory=memory, process_status=process["status"],
                                      fit_stage=stage), ensure_ascii=False) + "\n")
        if stage["stage"] != previous:
            previous = stage["stage"]
            with (root / "results/RUN_SUMMARY.md").open("a") as summary:
                summary.write(f"\n{stamp}：实测{method}阶段为{previous}；"
                              f"可用内存{memory['MemAvailable_kB']/1024**2:.2f} GiB。\n")
            print(stamp, method, previous, flush=True)
        if process["status"] != "running":
            break
        time.sleep(15)
