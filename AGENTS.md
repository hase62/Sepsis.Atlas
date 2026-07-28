# Codex instructions

Read `PROJECT.md` and `handoffs/CURRENT.md` before starting.

## Scope

- Operate from this repository root.
- Access Dropbox only through these project links:
  - `workspace/research-input/`
  - `workspace/research-output/`
  - `workspace/large-input/`
  - `workspace/large-output/`
- Use `workspace/scratch/` only for local, reproducible temporary files.
- Do not inspect Dropbox roots, parent directories, sibling repositories, or other projects.
- Do not recursively enumerate all workspace links unless the task explicitly requires a complete inventory. Start from the named subdirectory or file relevant to the task.

## Storage rules

- `workspace/research-input/`: shared lightweight inputs, documents, metadata, and manually exchanged files; read/write.
- `workspace/research-output/`: shared lightweight or normal-sized outputs; write outputs into a task-specific subdirectory when practical.
- `workspace/large-input/`: shared large inputs; treat as read-only unless explicitly instructed otherwise.
- `workspace/large-output/`: shared large outputs; write into a task-specific subdirectory when practical.
- `workspace/scratch/`: local disposable scratch; never leave the only copy of a required artifact here.
- Commit persistent text, code, instructions, manifests, and small metadata to Git.

Follow all mandatory rules in `PROJECT.md`.

## Agent switching

This worktree may be continued later by Claude Code. Before making changes, inspect `git status` and `git diff`. Before handing off, update `handoffs/CURRENT.md` with completed work, remaining work, validation, and shared output paths. Do not run concurrently with Claude Code in the same worktree.
