"""Recover only the empty draft demonstrably created by our interrupted upload."""
import sys,json,subprocess,datetime
from pathlib import Path
sys.path.insert(0,str(Path('results/code').resolve()))
from pipeline_common import ROOT,digest,save_json
folder=ROOT/'return/pc-all-epithelial-results-v1-20261022'
state_file=folder/'NEW_RELEASE_UPLOAD_STATE.json'
assert not state_file.exists()
returned=json.loads((folder/'RETURN_MANIFEST.json').read_text())
commit=returned['code_commit'];tag='crc-all-epithelial-trajectory-results-v1-20261022'
def run(a):return subprocess.check_output(a,text=True).strip()
assert run(['git','rev-parse','HEAD'])==commit
assert not run(['git','status','--porcelain','--untracked-files=all','--',str(ROOT)])
release=json.loads(run(['gh','api','repos/COSMOS340/exchange/releases/420893762']))
user=json.loads(run(['gh','api','user']))['login']
assert release['draft'] and release['tag_name']==tag and release['target_commitish']==commit
assert release['author']['login']==user and not release['assets']
assert release['html_url'] in (ROOT/'return/publish_v1_console.log').read_text()
notes=(ROOT/'return/RELEASE_NOTES_v1.md').read_text()
assert release['body'].replace('\r\n','\n').strip()==notes.strip()
paths=[folder/x['name'] for x in returned['assets']]
paths += [folder/'RESULT_MANIFEST.json',folder/'SIDECAR_MANIFEST.json',folder/'RETURN_MANIFEST.json',folder/'LOCAL_TRANSPORT_VERIFICATION.json',ROOT/'results/RUN_SUMMARY.md']
expected=[dict(name=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in paths]
assert all(x['bytes']<2*1024**3 for x in expected)
proof=dict(reason='New draft was created successfully, then REST releases/tags lookup returned 404; recover by verified release ID',recovered_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),release_id=release['id'],draft_author=user,draft_created_at=release['created_at'],original_console_sha256=digest(ROOT/'return/publish_v1_console.log'),original_source_sha256=digest(ROOT/'results/code/publish_new_release.py'),own_draft_url_matches_original_console=True,exact_body_target_and_empty_assets_verified=True)
save_json(folder/'DRAFT_RECOVERY_EVIDENCE.json',dict(proof=proof,release=release))
save_json(state_file,dict(tag=tag,release_id=release['id'],code_commit=commit,branch='pc-all-epithelial-results-20261022',expected_assets=expected,verified_assets=[],published=False,draft_recovery=proof))
print('Recovered our exact empty draft',release['id'],'without mutating GitHub or existing assets.')
