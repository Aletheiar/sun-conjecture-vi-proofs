"""Package the explicit proof whitelist and verify it in a clean directory."""
from pathlib import Path
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
import sympy

root=Path(__file__).resolve().parent
manifest_path=root/"proof_manifest_2026-10-25.json"
manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
for entry in manifest["files"]:
    p=(root/entry["path"]).resolve()
    assert p.is_relative_to(root),p
    assert digest(p)==entry["sha256"],entry["path"]
temp_parent=root/"tmp"
temp_parent.mkdir(exist_ok=True)
stage=Path(tempfile.mkdtemp(prefix="proof-review-",dir=temp_parent)).resolve()
assert stage.is_relative_to(temp_parent.resolve())
for entry in manifest["files"]:
    target=stage/entry["path"]
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root/entry["path"],target)
shutil.copy2(manifest_path,stage/manifest_path.name)
validation=stage/"ISOLATED_VALIDATION"
validation.mkdir()
records=[]
scripts=["verify_vi12_elliptic_reduction_210.py","verify_vi12_hypergeometric_211.py",
         "verify_vi3_proof_220.py","verify_vi3_literal_certificates_223.py"]
for script in scripts:
    result=subprocess.run([sys.executable,"-I",str(stage/script)],
                          cwd=stage,capture_output=True,text=True)
    (validation/(Path(script).stem+".log")).write_text(result.stdout+result.stderr,encoding="utf-8")
    records.append({"script":script,"exit_code":result.returncode})
    print(script,"exit",result.returncode,flush=True)
    if result.returncode:
        print(result.stdout+result.stderr,flush=True)
        raise SystemExit(result.returncode)
for entry in manifest["files"]:
    assert digest(stage/entry["path"])==entry["sha256"],"Verification modified "+entry["path"]
report={"generated_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "python":sys.version.split()[0],"sympy":sympy.__version__,
        "isolated_mode":True,"whitelist_files":len(manifest["files"]),
        "verification":records,
        "scope":"Reproducible exact checks only; independent AI review and response completed; professional mathematical review pending"}
(validation/"REPORT.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
archive=root/"PROOF_REVIEW_BUNDLE_2026-10-25.zip"
with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED) as zf:
    for p in sorted(stage.rglob("*")):
        if p.is_file():
            zf.write(p,p.relative_to(stage).as_posix())
report["archive"]=archive.name
report["archive_sha256"]=digest(archive)
(root/"proof_bundle_validation_224.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print("Isolated proof bundle ready:",archive.name,flush=True)
print("Files on whitelist:",len(manifest["files"]),"archive bytes:",archive.stat().st_size,flush=True)
