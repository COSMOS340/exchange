"""Read only the current task's bound child process and resource counters."""
import argparse
import datetime
import json
import os
from pathlib import Path
import time

from pipeline_common import ROOT, LOG, BLOCK, save_json

parser = argparse.ArgumentParser()
parser.add_argument('method', choices=['Monocle2', 'Slingshot'])
method = parser.parse_args().method
record_file = LOG / (method + '_local_process.json')
samples = LOG / (method + '_resource_samples.jsonl')
while True:
    if not record_file.exists():
        time.sleep(2)
        continue
    record = json.loads(record_file.read_text())
    row = dict(time=datetime.datetime.now().astimezone().isoformat(), method=method,
               status=record['status'], stage='starting')
    progress = LOG / (BLOCK + '_' + method + '_progress.json')
    if progress.exists():
        row['stage'] = json.loads(progress.read_text()).get('stage')
    # Tool sessions use separate PID namespaces. On the host, locate only the
    # bound R script in this exact task directory, never another project.
    rpid = None
    expected = str(ROOT / ('work/fit_generated.R' if method == 'Monocle2'
                          else 'results/code/fit_slingshot_exact_entry.R'))
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit():
            continue
        try:
            args = (directory / 'cmdline').read_bytes().decode(errors='replace').split('\0')
            if args and Path(args[0]).name in ['R', 'Rscript'] and any(
                    a == expected or a == '--file=' + expected for a in args):
                if (directory / 'cwd').resolve() == ROOT and method in args:
                    if rpid is not None:
                        raise RuntimeError('More than one matching formal R child; inspect before continuing')
                    rpid = directory.name
        except OSError:
            pass
    try:
        if rpid is not None:
            info = Path('/proc/' + rpid + '/status').read_text()
            fields = {a: b.strip() for a, b in (s.split(':', 1) for s in info.splitlines())}
            row.update(R_pid=int(rpid), R_RSS_GiB=int(fields['VmRSS'].split()[0])/1024**2,
                       R_peak_RSS_GiB=int(fields['VmHWM'].split()[0])/1024**2)
            stat = Path('/proc/' + rpid + '/stat').read_text().rsplit(')', 1)[1].split()
            row['CPU_seconds'] = (int(stat[11]) + int(stat[12])) / os.sysconf('SC_CLK_TCK')
    except (OSError, KeyError, ValueError):
        pass
    memory = {a: int(b.strip().split()[0]) for a, b in (s.split(':', 1) for s in Path('/proc/meminfo').read_text().splitlines())}
    row.update(MemAvailable_GiB=memory['MemAvailable']/1024**2,
               swap_used_GiB=(memory['SwapTotal']-memory['SwapFree'])/1024**2)
    with samples.open('a') as f:
        f.write(json.dumps(row) + '\n')
    save_json(LOG / (method + '_monitor_latest.json'), row)
    print(json.dumps(row), flush=True)
    if record['status'] != 'running':
        break
    time.sleep(20)
