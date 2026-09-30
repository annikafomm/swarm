# CLAUDE.md

Instructions for any Claude Code session working in this repository.

## NFS I/O — this workspace is NFS-mounted, and shared

`/workspaces/swarm` is an NFS mount, not local disk, and multiple Claude Code sessions run
against it concurrently. A broad filesystem scan here doesn't just cost the session that ran
it — it saturates the NFS server and stalls I/O for every other session on the box. Observed
damage: one session ran `bfs / -name '*.h5ad' -size +10M`; another session's reads of a 128 MB
h5ad that normally take seconds did not complete in ten minutes, and its processes parked in
uninterruptible D-state (unkillable until the I/O returns) with load average hitting 232.
Treat every wide scan as expensive and shared, not free and local.

**Never:**
- `find` or `bfs` rooted at `/`, `$HOME`, `/workspaces`, or the repo root as a whole.
- `grep -r` / `grep -R` / `ugrep -r` / a shelled-out ripgrep rooted at `.` from a large
  directory (e.g. `frontend/src/app`) or the repo root — including reflexively adding `.` as
  "just to be safe" when a specific file or a couple of files would answer the question.
- A recursive shell glob over a large tree (`**/*.ts`, `data/**/*.h5ad`).
- `git log --all` / `git log --diff-filter=A --name-only` over the whole history when looking
  for one file or fact — scope with a path or `-- <path>` first.
- Reading an entire AnnData (`sc.read_h5ad()`) over NFS to use a fraction of it — slice the
  needed columns with `h5py`, or stage the file to `/tmp` (overlay, not NFS) once and re-read
  it from there.
- Repeating the same expensive scan more than once in a session on the theory that the first
  result "might be stale" — cache what you already found.
- Handing a subagent (Explore or otherwise) an open-ended "search the whole repo for X" task —
  it inherits the same problem. Give it the narrowest plausible directory to start from.

**Instead:**
- Scope every search to the narrowest directory that can hold the answer, and prefer
  `-maxdepth`.
- Prefer the built-in Grep/Glob tools (already scoped, and typically faster than a shelled
  walk) over `find`/`grep -r`. When a genuinely broad, open-ended search is needed, delegate to
  an agent with an explicit narrow starting path rather than running a raw recursive shell
  command yourself.
- For a repo-wide code search, prefer `git grep <pattern> -- <path-prefix>` over a working-tree
  walk — it searches the (much cheaper) git index rather than statting every file on disk.
- To locate a known data file, read its path out of `backend/data/*/config.json` instead of
  searching for it.
- Load average lies on this box: D-state (uninterruptible I/O wait) processes count toward it,
  so 200+ usually means a saturated NFS mount, not CPU contention — don't read it as compute
  pressure and don't launch a duplicate job assuming the first one died.
