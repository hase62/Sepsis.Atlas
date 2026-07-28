# Sepsis.Atlas

Research code and documentation for **Sepsis.Atlas**.

## Local setup

```bash
setup-workspace
research-doctor Sepsis.Atlas
```

`setup-workspace` creates five local links:

```text
workspace/research-input  -> Research/aicode/input/Sepsis.Atlas
workspace/research-output -> Research/aicode/output/Sepsis.Atlas
workspace/large-input     -> ForShareLargeData/aicode/input/Sepsis.Atlas
workspace/large-output    -> ForShareLargeData/aicode/output/Sepsis.Atlas
workspace/scratch         -> local ~/scratch/... directory
```

Agents must not search outside these paths. Required files should be copied into the appropriate project directory in advance or created inside it.

## Open in Visual Studio Code

```bash
cd ~/src/Sepsis.Atlas
code .
```

Use the Codex or Claude Code VS Code extension as the primary interface. The project settings exclude linked data and output trees from automatic VS Code search and file watching; files can still be opened directly.

## Long task worktree

```bash
new-worktree Sepsis.Atlas shared metadata-audit
cd ~/worktrees/Sepsis.Atlas/shared-metadata-audit
code .
```

Use one unique task name per logical task. Continue using that name until merge, then delete the worktree and branch. `shared` means Codex and Claude Code may use the same worktree sequentially, not simultaneously.

## Continue on another computer

Before leaving the current computer:

```bash
git add <reviewed-files>
git commit -m "WIP: checkpoint current task"
git push -u origin work/metadata-audit
```

On the next computer:

```bash
cd ~/src/Sepsis.Atlas
git fetch --all --prune
new-worktree Sepsis.Atlas shared metadata-audit
```

GitHub carries code and instructions. Dropbox carries the four shared project directories. Uncommitted changes, environments, local scratch, and Agent chat sessions do not move between computers.
