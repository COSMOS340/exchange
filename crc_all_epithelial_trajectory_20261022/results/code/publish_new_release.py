"""Publish a new, exact-commit result Release with per-asset remote SHA checks.

Run only after local acceptance/transport verification. A saved upload state
allows safe resumption of this newly created draft; unrelated releases cannot
be overwritten and existing asset names are never clobbered.
"""
import argparse,json,subprocess
from pathlib import Path
from pipeline_common import ROOT,digest,save_json

REPO='COSMOS340/exchange'


def run(command):
    return subprocess.check_output(command,text=True).strip()


def api(endpoint):return json.loads(run(['gh','api',endpoint]))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path)
    parser.add_argument('--tag',required=True);parser.add_argument('--notes-file',type=Path,required=True)
    args=parser.parse_args();folder=args.folder.resolve()
    assert folder.is_relative_to(ROOT/'return')
    returned=json.loads((folder/'RETURN_MANIFEST.json').read_text())
    verification=json.loads((folder/'LOCAL_TRANSPORT_VERIFICATION.json').read_text());assert verification['passed']
    commit=returned['code_commit']
    assert run(['git','rev-parse','HEAD'])==commit
    assert not run(['git','status','--porcelain','--untracked-files=all','--',str(ROOT)])
    branch=run(['git','branch','--show-current']);assert branch=='pc-all-epithelial-results-20261022'
    assert args.notes_file.is_file()
    paths=[folder/x['name'] for x in returned['assets']]
    paths += [folder/'RESULT_MANIFEST.json',folder/'SIDECAR_MANIFEST.json',folder/'RETURN_MANIFEST.json',folder/'LOCAL_TRANSPORT_VERIFICATION.json',ROOT/'results/RUN_SUMMARY.md']
    expected=[dict(name=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in paths]
    assert len({r['name'] for r in expected})==len(expected)
    assert all(r['bytes']<2*1024**3 for r in expected)
    state_file=folder/'NEW_RELEASE_UPLOAD_STATE.json'
    if state_file.exists():
        state=json.loads(state_file.read_text())
        assert state['tag']==args.tag and state['code_commit']==commit and state['expected_assets']==expected
        release=api(f'repos/{REPO}/releases/{state["release_id"]}')
        assert release['tag_name']==args.tag and release['target_commitish']==commit
    else:
        # The user has authorized the independent result branch and new Release.
        subprocess.run(['git','push','origin',branch],check=True)
        assert api(f'repos/{REPO}/commits/{branch}')['sha']==commit
        releases=api(f'repos/{REPO}/releases?per_page=100')
        if args.tag in {x['tag_name'] for x in releases}:
            raise RuntimeError('Tag already exists: choose the next unused version; never overwrite it')
        subprocess.run(['gh','release','create',args.tag,'--repo',REPO,'--target',commit,'--draft',
            '--title','All epithelial trajectory | 223,736 cells | PC results',
            '--notes-file',str(args.notes_file)],check=True)
        release=api(f'repos/{REPO}/releases/tags/{args.tag}')
        assert release['draft'] and release['target_commitish']==commit
        state=dict(tag=args.tag,release_id=release['id'],code_commit=commit,branch=branch,
            expected_assets=expected,verified_assets=[],published=False)
        save_json(state_file,state)
    # Verify or upload each complete asset serially; no --clobber.
    for path,record in zip(paths,expected):
        assert path.stat().st_size==record['bytes'] and digest(path)==record['sha256']
        release=api(f'repos/{REPO}/releases/{state["release_id"]}')
        found=[a for a in release['assets'] if a['name']==record['name']]
        if not found:
            subprocess.run(['gh','release','upload',args.tag,str(path),'--repo',REPO],check=True)
            release=api(f'repos/{REPO}/releases/{state["release_id"]}')
            found=[a for a in release['assets'] if a['name']==record['name']]
        assert len(found)==1
        asset=found[0]
        assert asset['state']=='uploaded' and asset['size']==record['bytes']
        assert asset.get('digest')=='sha256:'+record['sha256'],record['name']+' remote digest mismatch'
        state['verified_assets']=[a for a in state['verified_assets'] if a['name']!=record['name']]+[
            dict(**record,asset_id=asset['id'],remote_digest=asset['digest'],download_url=asset['browser_download_url'])]
        save_json(state_file,state);print('Remote SHA verified',record['name'],flush=True)
    release=api(f'repos/{REPO}/releases/{state["release_id"]}')
    assert {a['name'] for a in release['assets']}=={r['name'] for r in expected}
    if release['draft']:
        subprocess.run(['gh','release','edit',args.tag,'--repo',REPO,'--draft=false'],check=True)
    release=api(f'repos/{REPO}/releases/tags/{args.tag}')
    assert not release['draft'] and release['target_commitish']==commit
    assert api(f'repos/{REPO}/commits/{args.tag}')['sha']==commit
    for record in expected:
        asset=next(a for a in release['assets'] if a['name']==record['name'])
        assert asset['size']==record['bytes'] and asset.get('digest')=='sha256:'+record['sha256']
    state.update(published=True,url=release['html_url'],all_remote_names_bytes_SHA256_verified=True,
        exact_tag_code_commit_verified=True,Mac_receiving_acceptance='pending')
    save_json(state_file,state)
    save_json(folder/'GITHUB_PUBLICATION_RECEIPT.json',state)
    print('Published and fully verified:',release['html_url'],flush=True)


if __name__=='__main__':main()
