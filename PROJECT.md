# Sepsis.Atlas

## Goal

TODO: この研究の目的を数行で記載する。

## Current phase

TODO: 現在の作業段階を記載する。

## Repository and storage scope

- Code, tests, configuration, Agent instructions, handoff notes, manifests, and small metadata are stored in Git.
- Dropbox access is limited to four fixed project directories exposed under `workspace/`.
- `workspace/` is a local link layer; the authoritative shared files remain in Dropbox.
- Local temporary files are placed only under `workspace/scratch/`.

## Fixed workspace paths

| Workspace path | Dropbox target | Intended use |
|---|---|---|
| `workspace/research-input/` | `Research/aicode/input/Sepsis.Atlas` | Lightweight inputs, documents, metadata, and shared staging |
| `workspace/research-output/` | `Research/aicode/output/Sepsis.Atlas` | Lightweight or normal-sized shared outputs |
| `workspace/large-input/` | `ForShareLargeData/aicode/input/Sepsis.Atlas` | Large shared inputs; read-only by default |
| `workspace/large-output/` | `ForShareLargeData/aicode/output/Sepsis.Atlas` | Large shared outputs |
| `workspace/scratch/` | Local `~/scratch/...` | Disposable and reproducible temporary files |

Subdirectories may be created within these project directories. Prefer task- or run-specific output subdirectories to avoid collisions between worktrees and computers.

## Mandatory rules

- Work only inside this repository and the five paths listed above.
- Do not access or search Dropbox roots directly.
- Do not search parent directories, sibling repositories, or unrelated projects.
- Do not recursively inventory every linked directory unless explicitly requested.
- Treat `workspace/large-input/` as read-only unless explicitly instructed otherwise.
- Do not read credentials, tokens, secrets, or unrelated `.env` files.
- Do not commit, push, merge, rebase, or delete shared data unless explicitly requested.
- Do not leave the only copy of a persistent artifact in `workspace/scratch/`.
- Inspect `git status` and `git diff` before completing a task.
- Codex and Claude Code may continue the same task sequentially in one worktree, but must not edit it at the same time.
- When handing off, update `handoffs/CURRENT.md`; Git files and shared paths are the source of truth, not chat history.
