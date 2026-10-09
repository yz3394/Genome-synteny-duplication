# Local skill and GitHub synchronization

## Authorization and scope

The current version has been published. The user subsequently cancelled hourly update checks and requires an explicit request before future GitHub updates. Each subsequent synchronization requires a new explicit request; local skill edits alone do not trigger checks or uploads.

- **Local source:** `~/.codex/skills/genome-synteny-duplication/`.
- **Repository destination:** `skills/genome-synteny-duplication/`.
- **Only remote target:** `https://github.com/yz3394/Genome-synteny-duplication.git`.
- **Direction:** installed local skill → this repository; do not automatically overwrite the local installation with remote content.

Synchronization covers only the skill directory above and its package manifest. The README, maintenance documentation, repository configuration and synchronization tools are not copied from the local skill directory; review any changes to them separately. Do not package or upload raw data, analysis outputs, run logs, temporary files, credentials or private links from neighboring research projects. Existing validation summaries within the skill are package documentation.

## Manual initiation

The former hourly heartbeat `skill-github` has been deleted. No periodic update checks are scheduled.

Validate and synchronize using the steps below only when the user explicitly requests a GitHub update or an equivalent action. Routine analysis or local skill editing does not trigger publication. On completion, report the actual commit and remote verification results. On failure, explain the specific cause and retain local changes.

## Synchronization steps

1. Verify the local source, working copy and remote URL, and inspect the current Git status. Do not include changes of unknown origin or unfinished work from other tasks in the synchronization commit.
2. Fetch the remote state and inspect history for changes. If there are remote-only edits, diverged histories or changes that conflict with the local package, stop automatic overwriting and report the issue for review and resolution.
3. Compare the source with the repository package. Check that added, modified and removed files are within the skill's scope, and review them for sensitive content or research data.
4. Run the synchronization tool to prepare the repository copy and SHA256 package manifest and execute its validation. Then review the actual Git diff and confirm that it matches the local updates.
5. Run tests relevant to the changes and record both completed checks and explicit skips. Stop before committing or pushing if script structure, references, content scope or tests fail validation.
6. Commit only the reviewed package and manifest changes, then push normally to the designated repository. Do not force-push or automatically create, move or delete tags.
7. Read the remote state. Report successful GitHub synchronization only after confirming that the remote branch points to the intended local commit. If pushing fails, retain local changes and commits and state clearly that synchronization is incomplete.

## Synchronization tool

Run from the repository root:

```bash
python3 tools/sync_installed_skill.py --source "$HOME/.codex/skills/genome-synteny-duplication" --check
python3 tools/sync_installed_skill.py --source "$HOME/.codex/skills/genome-synteny-duplication"
```

`--check` compares the source and repository package; without `--check`, the tool prepares and validates the package copy in the repository. **The script does not push to GitHub.** The agent performing the task executes the Git commit and push separately after validation passes, the diff has been reviewed and the remote has been checked for conflicts.

The SHA256 manifest records package file contents to verify copying and detect changes. Matching hashes establish content identity only; they do not replace sensitive-content checks, scientific-method review or execution-based validation.

## Preservation and conflict handling

- The installed local skill is the synchronization source. Preserve its updates; do not automatically revert it because the remote version differs.
- Existing versions remain in Git history. Do not rewrite history or resolve conflicts by deleting the working copy.
- Inspect unexplained file removals or extensive changes first. Package removals must match the actual local updates and reviewed changes.
- Do not expose authentication tokens or write credentials to configuration, logs or commits. Use the existing GitHub authentication mechanism.
- Obtain new explicit authorization before changing the synchronization direction or target repository, or expanding the scope to other data directories.
