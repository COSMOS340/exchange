"""Retry only transient network failures of the unchanged, resumable publisher."""
import sys,time,json,subprocess,datetime,hashlib
from pathlib import Path
root=Path.cwd();source=root/'results/code/publish_new_release.py'
command=[sys.executable,str(source),'return/pc-all-epithelial-results-v1-20261022','--tag','crc-all-epithelial-trajectory-results-v1-20261022','--notes-file','return/RELEASE_NOTES_v1.md']
receipt={'publisher_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'attempts':[],'published':False}
for attempt in range(3,9):
 console=root/f'return/publish_v1_retry{attempt:02d}_console.log'
 resource=root/f'return/publish_v1_retry{attempt:02d}_GNU_time.txt'
 actual=['/usr/bin/time','-v','-o',str(resource)]+command
 item=dict(command=actual,started_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),attempt=attempt)
 print('Starting publication attempt',attempt,flush=True);start=time.monotonic()
 with console.open('x') as output:r=subprocess.run(actual,stdout=output,stderr=subprocess.STDOUT,cwd=root)
 item.update(returncode=r.returncode,elapsed_seconds=time.monotonic()-start,finished_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),console_sha256=hashlib.sha256(console.read_bytes()).hexdigest(),GNU_time_sha256=hashlib.sha256(resource.read_bytes()).hexdigest())
 receipt['attempts'].append(item);receipt['published']=r.returncode==0
 (root/'return/PUBLICATION_RETRY_COMMANDS.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if r.returncode==0:
  print('Publication completed with all remote SHA checks.',flush=True);sys.exit(0)
 text=console.read_text().lower()
 transient=any(x in text for x in ['connection reset','i/o timeout','tls handshake timeout','unexpected eof','http 502','http 503','http 504','context deadline exceeded'])
 if not transient:raise RuntimeError('Nontransient failure; preserved logs, do not blindly retry: '+str(console))
 print('Transient network failure, retrying the saved upload state.',flush=True)
 time.sleep(min(2**(attempt-2),10))
raise RuntimeError('Network retry budget reached; preserve state for diagnosis')
