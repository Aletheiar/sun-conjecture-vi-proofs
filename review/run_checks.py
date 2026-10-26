"""Reproduce the nine independent AI-review scripts, retaining fresh logs."""
from pathlib import Path
import json
import subprocess
import sys

here = Path(__file__).resolve().parent
out = here/'reproduced'
out.mkdir(exist_ok=True)
records = []
for script in sorted((here/'deepseek').glob('0[1-9]_*.py')):
    p = subprocess.run([sys.executable, '-I', str(script)], cwd=here.parent,
                       capture_output=True, text=True, encoding='utf-8', errors='replace')
    (out/(script.stem+'.log')).write_text(p.stdout+p.stderr, encoding='utf-8')
    records.append({'script': script.name, 'exit_code': p.returncode})
    print(script.name, 'exit', p.returncode, flush=True)
    if p.returncode:
        break
(out/'summary.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
raise SystemExit(0 if len(records)==9 and all(r['exit_code']==0 for r in records) else 1)
