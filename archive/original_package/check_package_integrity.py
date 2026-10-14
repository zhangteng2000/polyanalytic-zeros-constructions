"""Check the SHA-256 manifest before rerunning scripts that write reports."""
from pathlib import Path
import hashlib
import json

base = Path(__file__).resolve().parent
manifest = json.loads((base / "SHA256.json").read_text(encoding="utf-8"))
failures = []
for name, expected in manifest["files"].items():
    path = base / name
    if not path.is_file():
        failures.append(name + ": missing")
    elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        failures.append(name + ": SHA-256 mismatch")
if failures:
    raise SystemExit("\n".join(failures))
print("PASS:", len(manifest["files"]), "files match the delivery manifest.")
