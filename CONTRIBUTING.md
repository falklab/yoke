# Contributing to yoke

## Development setup

### Option 1: VS Code devcontainer (recommended)
Open the repository in VS Code and choose **Dev Containers: Reopen in Container**.
The container has Python 3.12, [uv](https://docs.astral.sh/uv/) and `rtk`. On first start it creates the
environment (`.venv`) and installs the pre-commit hooks. Files are formatted on save.

### Option 2: Local
Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then from the repository root
```
uv sync                          # create .venv with yoke and the dev tools
uv run pre-commit install        # run checks on each commit
git config blame.ignoreRevsFile .git-blame-ignore-revs
```
`rtk` is only needed for `yoke module_score` and its tests; see the README for installing it.

### Everyday commands
| Task | Command |
|---|---|
| Run the CLI | `uv run yoke ...` |
| Run tests | `uv run pytest` |
| Tests with coverage | `uv run pytest --cov` (add `--cov-report=html` for a browsable report in `htmlcov/`) |
| Lint (and fix what can be fixed) | `uv run ruff check --fix .` |
| Format | `uv run ruff format .` |
| Run all pre-commit checks | `uv run pre-commit run --all-files` |
| Add a dependency | `uv add <package>` (dev tools: `uv add --group dev <package>`) |

Always change dependencies with `uv add` / `uv remove` so that `pyproject.toml` and `uv.lock` stay in sync, and
commit both.

## Code style
Style is enforced by [ruff](https://docs.astral.sh/ruff/), configured in `pyproject.toml`:
PEP 8 with a 79 character line length, double quotes, sorted imports, modern type hints and Google-style
docstrings. The formatter is the final word on layout, so there is no need to format by hand.
More detailed conventions are written down in [CLAUDE.md](CLAUDE.md), which applies to people as well as
AI agents.

## Pre-commit hooks
On each `git commit` the hooks:
- fix whitespace and file endings, and check YAML/TOML files,
- run `ruff check --fix` and `ruff format`,
- run the fast tests for the Python files being committed.

If a hook modifies a file, the commit stops: review the changes, `git add` them and commit again.

The test hook maps changed files to tests by name: `src/yoke/<module>.py` runs `tests/test_<module>.py`,
a changed test file runs itself, and `tests/conftest.py` runs everything. Tests marked `slow`, `local_data`
or `requires_rtk` are not run by the hook. CI always runs the full suite.

In an emergency `git commit --no-verify` skips the hooks; CI will still check the change.

## Tests
- Tests live in `tests/test_<module>.py`, mirroring `src/yoke/<module>.py`.
- Small input files go in `tests/data/`, expected outputs in `tests/expected/<test_name>/`. Prefer building
  small DataFrames in fixtures over adding files.
- Mark tests that need extra resources, so they skip cleanly elsewhere:
  - `@pytest.mark.local_data`: needs large data that is not in the repository. These run only when
    `YOKE_LOCAL_DATA` points to a directory with that data
    (`YOKE_LOCAL_DATA=/path/to/data uv run pytest`). Currently expected layout:
    `DIME/{MAGvsGC.txt.gz, DIAass_KGM.srt.gzgeneAss.gz, KEGG.descr, KGML0.txt}` and
    `output/KEGG/MGS.*.tsv.gz`.
  - `@pytest.mark.requires_rtk`: needs `rtk` on `PATH`.
  - `@pytest.mark.slow`: takes more than a couple of seconds.
- CI fails if total coverage drops below the `fail_under` value in `pyproject.toml`. Raise it as coverage
  improves.

## Branching
We use [GitHub Flow](https://docs.github.com/en/get-started/using-github/github-flow):

1. `main` is always in a working, releasable state. Don't commit to it directly.
2. Create a short-lived branch from `main` for each piece of work, named by type:
   `feat/<topic>`, `fix/<topic>`, `docs/<topic>`, `chore/<topic>`.
3. Open a pull request early (a draft PR is fine). CI runs lint and tests on every push.
4. Get one review from another contributor. A quick look together in person is fine; leave a comment or
   approval on the PR so there is a record.
5. **Squash and merge**. The PR title becomes the single commit on `main`, so it must be a Conventional
   Commit (checked by CI). Delete the branch afterwards.

There are no `develop` or long-lived release branches. Urgent fixes are just a `fix:` PR.
Keep branches short-lived and rebase on or merge `main` if they fall behind.

## Commit messages and PR titles
Follow [Conventional Commits](https://www.conventionalcommits.org/): `type(scope): summary`.

| Type | Use for | Version bump (before 1.0) |
|---|---|---|
| `feat` | New functionality, e.g. a new subcommand or option | minor (0.**x**.0) |
| `fix` | Bug fixes | patch (0.0.**x**) |
| `perf` | Performance improvements | patch |
| `docs`, `test`, `refactor`, `build`, `ci`, `chore` | Everything else | none |

The scope is usually the subcommand or module, e.g. `feat(extract_contigs): support multiple binners` or
`fix(annotations): handle empty MGS`. Mark breaking changes (e.g. renamed CLI options, changed output format)
with `!` (`feat!: ...`) or a `BREAKING CHANGE:` footer.

Commits on feature branches can be informal, because only the squashed PR title lands on `main`.

## Versioning and releases
yoke follows [Semantic Versioning](https://semver.org/). While the version is `0.x`, breaking changes bump
the minor version. Move to `1.0.0` once the CLI and its outputs are considered stable.

Releases are automated with [release-please](https://github.com/googleapis/release-please):

1. After each merge to `main`, release-please updates an open PR titled `chore(main): release X.Y.Z`.
   It contains the next version number (from the commit types since the last release) and the
   `CHANGELOG.md` entries.
2. When you want to release, merge that PR. Several features can be batched into one release.
3. Merging tags `vX.Y.Z`, creates a GitHub Release with the changelog, and attaches the built wheel.

Never edit the version in `pyproject.toml` or `CHANGELOG.md` by hand; release-please owns them.
Users install a release with `pip install git+https://github.com/falklab/yoke.git@vX.Y.Z`.

## Repository settings (admin, one-off)
- **Settings → General → Pull Requests**: allow squash merging only, with the default commit message set to
  "Pull request title"; enable "Automatically delete head branches".
- **Settings → Actions → General → Workflow permissions**: enable "Allow GitHub Actions to create and approve
  pull requests" (needed by release-please).
- **Settings → Rules → Rulesets** (requires a paid GitHub plan for private repositories): protect `main`,
  require a pull request, and require the `lint`, `test` and `conventional-title` checks. On the free plan,
  rely on the conventions above.
- Release PRs are opened by the Actions bot, so CI does not run on them. This is fine, because they only
  change the version and changelog.
