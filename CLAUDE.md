@PROJECT.md

Read `handoffs/CURRENT.md` before starting.

# Claude Code instructions

- Operate from this repository root.
- Access Dropbox only through:
  - `workspace/research-input/`
  - `workspace/research-output/`
  - `workspace/large-input/`
  - `workspace/large-output/`
- Use `workspace/scratch/` only for local, reproducible temporary files.
- Do not inspect Dropbox roots, parent directories, sibling repositories, or other projects.
- Do not recursively enumerate all workspace links unless the task explicitly requires a complete inventory.
- Treat `workspace/large-input/` as read-only unless explicitly instructed otherwise.
- Put persistent outputs in `workspace/research-output/` or `workspace/large-output/`, preferably below a task-specific subdirectory.
- Do not leave the only copy of a required artifact in local scratch.

## Agent switching

This worktree may be continued later by Codex. Before making changes, inspect `git status` and `git diff`. Before handing off, update `handoffs/CURRENT.md` with completed work, remaining work, validation, and shared output paths. Do not run concurrently with Codex in the same worktree.
