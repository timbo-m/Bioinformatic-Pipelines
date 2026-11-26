"""Small local settings file: ~/.famscan/config.json (or $FAMSCAN_CONFIG).

Holds only what a run needs and should not be typed every time - currently the
contact email that the EBI's BLAST service requires. It is written on your own
machine with owner-only permissions where the OS supports it, and is never part
of the repository: keep your address out of committed files and out of shared
command lines.
"""
from __future__ import annotations

import json
import os
import re
import stat
from pathlib import Path

EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


def config_path() -> Path:
    p = Path(os.environ.get("FAMSCAN_CONFIG") or Path.home() / ".famscan" / "config.json")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def load() -> dict:
    p = config_path()
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def save(**values) -> Path:
    p = config_path()
    data = load()
    data.update({k: v for k, v in values.items() if v is not None})
    p.write_text(json.dumps(data, indent=2))
    try:
        p.chmod(stat.S_IRUSR | stat.S_IWUSR)
    except OSError:
        pass
    return p


def valid_email(value: str | None) -> bool:
    return bool(value and EMAIL_RE.match(value.strip()))


def email(explicit: str | None = None) -> str | None:
    """Contact email: the one passed in, else $FAMSCAN_EMAIL, else the saved one."""
    for candidate in (explicit, os.environ.get("FAMSCAN_EMAIL"), load().get("email")):
        if valid_email(candidate):
            return candidate.strip()
    return None
