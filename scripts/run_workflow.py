from pathlib import Path
import subprocess
import sys

project = Path(__file__).resolve().parent.parent
logs = project / "logs"
logs.mkdir(exist_ok=True)

scripts = sorted((project / "scripts").glob("[0-9][0-9]_*.py"))
if len(scripts) != 19:
    raise SystemExit("Expected numbered scripts 00-18.")

for script in scripts:
    print(f"START: {script.name}", flush=True)
    log = logs / f"{script.stem}.log"

    with log.open("w") as handle:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=project,
            stdout=handle,
            stderr=subprocess.STDOUT
        )

    if result.returncode:
        print(f"FAILED: {script.name}\nLog: {log}")
        print("\n".join(
            log.read_text(errors="replace").splitlines()[-30:]
        ))
        raise SystemExit(result.returncode)

    print(f"COMPLETED: {script.name}", flush=True)

print("ALL NUMBERED SCRIPTS COMPLETED.", flush=True)
