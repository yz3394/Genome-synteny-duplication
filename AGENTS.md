# Skill publication and maintenance

This repository publishes `skills/genome-synteny-duplication`. Read its entrypoint and the affected references before changing scientific behavior. The user's 2026-09-28 request authorizes publishing this Skill and synchronizing future local Skill updates to `https://github.com/yz3394/Genome-synteny-duplication.git`.

- The owner's installed `~/.codex/skills/genome-synteny-duplication` is the source for the configured local-to-GitHub sync. Historical copies inside study projects remain frozen.
- Follow `MAINTENANCE.md`. `tools/sync_installed_skill.py --check` reports content differences; without `--check` it validates a snapshot and prepares the repository copy. It does not commit or push.
- Synchronize only the Skill and its necessary repository documentation/validation files. Keep raw genomes, project outputs, local logs, personal paths and credentials outside this repository.
- Run relevant tests before publishing; report skips and external-tool limitations. For scientific or code changes, update the Skill version and describe the actual change. Packaging changes alone do not change the Skill version.
- Fetch and inspect the branch before a push. Preserve remote work, uncommitted edits and existing release tags. If the remote or local changes conflict, stop the sync and report the specific issue; do not force-push or reset.
- Confirm the remote commit and the published Skill hashes. Local files, validation, commit and successful remote publication are separate states.

These instructions do not opt the research project into project-state maintenance.
