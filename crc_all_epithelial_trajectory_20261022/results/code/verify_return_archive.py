"""Stream-check transport parts, ZIP CRC and every manifest entry, without extraction."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile


def file_hash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('folder',type=Path)
    folder=parser.parse_args().folder.resolve()
    returned=json.loads((folder/'RETURN_MANIFEST.json').read_text())
    archive=folder/'results.zip'
    assert archive.stat().st_size==returned['archive_bytes']
    assert file_hash(archive)==returned['archive_sha256']
    joined=hashlib.sha256();total=0
    for asset in returned['assets']:
        part=folder/asset['name']
        assert part.parent==folder and part.stat().st_size==asset['bytes']
        assert part.stat().st_size<2*1024**3
        assert file_hash(part)==asset['sha256']
        with part.open('rb') as stream:
            for block in iter(lambda:stream.read(8*1024**2),b''):
                joined.update(block);total+=len(block)
    assert total==returned['archive_bytes'] and joined.hexdigest()==returned['archive_sha256']
    expected=json.loads((folder/'RESULT_MANIFEST.json').read_text())
    mapped={'results/'+item['path']:item for item in expected}
    count=0;bytes_checked=0
    with zipfile.ZipFile(archive) as zipped:
        names=zipped.namelist();assert len(names)==len(set(names))
        assert json.loads(zipped.read('RESULT_MANIFEST.json'))==expected
        assert {name for name in names if name.startswith('results/')}==set(mapped)
        for info in zipped.infolist():
            path=PurePosixPath(info.filename)
            assert not path.is_absolute() and '..' not in path.parts
            checksum=hashlib.sha256();size=0
            # Reading to EOF also verifies the ZIP member CRC in zipfile.
            with zipped.open(info) as stream:
                for block in iter(lambda:stream.read(8*1024**2),b''):
                    checksum.update(block);size+=len(block)
            assert size==info.file_size
            if info.filename in mapped:
                item=mapped[info.filename]
                assert size==item['bytes'] and checksum.hexdigest()==item['sha256']
                count+=1;bytes_checked+=size
    assert count==len(expected)
    report=dict(passed=True,archive_sha256=returned['archive_sha256'],
        archive_bytes=returned['archive_bytes'],transport_assets=len(returned['assets']),
        concatenated_transport_hash_matches=True,all_zip_member_CRCs_verified=True,
        result_files_verified=count,result_bytes_verified=bytes_checked,
        scientific_acceptance='Separate: inspect numeric_acceptance.json and visual_audit.json')
    (folder/'LOCAL_TRANSPORT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
