#!/usr/bin/env python3
"""Validate and mirror one installed Skill; Git commit/push remain explicit steps."""

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

NAME = "genome-synteny-duplication"
ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "skills" / NAME
MANIFEST = ROOT / "skill-manifest.json"
DIRECTORIES = {"agents", "assets", "references", "scripts", "tests"}
EXTENSIONS = {".md", ".py", ".json", ".yaml", ".yml"}
PRIVATE = re.compile(r"/Users/[^/\s]+/|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----|github_pat_[A-Za-z0-9_]+|gh[pousr]_[A-Za-z0-9]{20,}")


def inventory(root):
    files = {}
    if not root.exists():
        return files
    if root.is_symlink():
        raise ValueError(f"Refusing symlink directory: {root}")
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root)
        if "__pycache__" in rel.parts or path.name == ".DS_Store":
            continue
        if path.is_symlink():
            raise ValueError(f"Refusing symlink: {rel}")
        if not path.is_file():
            continue
        if rel.as_posix() != "SKILL.md" and (len(rel.parts) < 2 or rel.parts[0] not in DIRECTORIES or path.suffix not in EXTENSIONS):
            raise ValueError(f"File needs a publication-scope review: {rel}")
        data = path.read_bytes()
        content = data.decode("utf-8")
        if PRIVATE.search(content):
            raise ValueError(f"Private path or credential pattern needs review: {rel}")
        files[rel.as_posix()] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    return files


def version(root):
    content = (root / "SKILL.md").read_text(encoding="utf-8")
    if not content.startswith("---\n"):
        raise ValueError("SKILL.md has no YAML frontmatter")
    front = content.split("---", 2)[1]
    if not re.search(rf"^name:\s*{NAME}\s*$", front, re.M):
        raise ValueError("Unexpected Skill name")
    match = re.search(r'^\s+version:\s*[\"\x27]?([0-9]+\.[0-9]+\.[0-9]+)[\"\x27]?\s*$', front, re.M)
    if not match:
        raise ValueError("Missing semantic metadata.version")
    return match.group(1)


def validate(root, files):
    for rel in files:
        path = root / rel
        text = path.read_text(encoding="utf-8")
        if path.suffix == ".py":
            ast.parse(text, filename=rel)
        elif path.suffix == ".json":
            json.loads(text)
        elif path.suffix == ".md":
            for link in re.findall(r"\[[^\]]*\]\(([^)]+)\)", text):
                target = link.split("#", 1)[0]
                if not target or "://" in target or target.startswith("mailto:"):
                    continue
                resolved = (path.parent / target).resolve()
                if not resolved.is_relative_to(root.resolve()) or not resolved.exists():
                    raise ValueError(f"Broken or external local Markdown link in {rel}: {target}")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg")
    result = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=root, env=env, text=True, capture_output=True)
    sys.stderr.write(result.stdout + result.stderr)
    if result.returncode:
        raise ValueError("Snapshot tests failed; repository copy was not changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path.home() / ".codex" / "skills" / NAME)
    parser.add_argument("--check", action="store_true", help="Read-only diff; no tests, writes, commits or network")
    args = parser.parse_args()
    source = args.source.expanduser().absolute()
    if source == DEST or source.is_symlink() or not source.is_dir():
        raise ValueError("Source must be an existing, distinct, nonsymlink Skill directory")
    before = inventory(source)
    skill_version = version(source)
    current = inventory(DEST)
    payload = {"skill": NAME, "version": skill_version, "files": before}
    changed = sorted(k for k, v in before.items() if current.get(k) != v)
    removed = sorted(set(current) - set(before))
    manifest_changed = not MANIFEST.exists() or json.loads(MANIFEST.read_text()) != payload
    status = {"status": "changed" if changed or removed or manifest_changed else "unchanged", "version": skill_version, "changed_files": changed, "removed_files": removed, "manifest_changed": manifest_changed}
    if args.check or status["status"] == "unchanged":
        print(json.dumps(status, indent=2))
        return
    with tempfile.TemporaryDirectory(prefix="synteny-publish-") as tmp:
        snapshot = Path(tmp) / NAME
        for rel in before:
            dest = snapshot / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source / rel, dest)
        if inventory(snapshot) != before:
            raise ValueError("Source changed while copying; retry after editing finishes")
        validate(snapshot, before)
        if inventory(source) != before or inventory(snapshot) != before:
            raise ValueError("Source or snapshot changed during validation; no publication prepared")
        for rel in changed:
            target = DEST / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(snapshot / rel, target)
        for rel in removed:
            (DEST / rel).unlink()
        MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if inventory(DEST) != before:
        raise ValueError("Repository copy does not match validated source")
    status["status"] = "prepared_and_validated"
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        sys.exit(f"ERROR: {error}")
