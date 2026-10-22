"""Record actual document/export/provenance commands after scientific acceptance."""
import datetime
import json
import subprocess
import sys
import time
from pipeline_common import ROOT, LOG, digest, save_json


def main():
    assert json.loads((LOG / 'numeric_acceptance.json').read_text())['passed']
    assert json.loads((LOG / 'visual_audit.json').read_text())['passed']
    receipt = {'scientific_acceptance_already_passed': True, 'stages': []}
    for name in ['build_deliverable_docs.py', 'export_small_review_tables.py',
                 'record_delivery_provenance.py']:
        source = ROOT / 'results/code' / name
        console = LOG / (name.removesuffix('.py') + '_delivery_console.log')
        resource = LOG / (name.removesuffix('.py') + '_delivery_GNU_time.txt')
        command = ['/usr/bin/time', '-v', '-o', str(resource), sys.executable, str(source)]
        stage = dict(command=command, cwd=str(ROOT), source_sha256=digest(source),
                     started_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
        started = time.monotonic()
        with console.open('x') as output:
            completed = subprocess.run(command, cwd=ROOT, stdout=output, stderr=subprocess.STDOUT)
        stage.update(returncode=completed.returncode, elapsed_seconds=time.monotonic()-started,
                     finished_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     console_sha256=digest(console), GNU_time_sha256=digest(resource))
        receipt['stages'].append(stage)
        save_json(LOG / 'delivery_document_commands.json', receipt)
        if completed.returncode:
            raise RuntimeError('Document/export stage failed; preserve receipt: ' + name)
        print('Completed', name, flush=True)


if __name__ == '__main__':
    main()
