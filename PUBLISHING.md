# Publishing gapsmith

Two distribution channels, independent:

1. **GitHub Releases** — pre-built binaries, triggered by a git tag
2. **crates.io** — `cargo install gapsmith-cli`

## GitHub Release (pre-built binaries)

Triggered automatically by `.github/workflows/release.yml` when you
push a `v*` tag.

```bash
# 1. Bump versions in Cargo.toml workspace.package (0.1.0 → 0.2.0 etc.)
# 2. Commit the bump + any release-notes edits.
# 3. Tag + push.
git tag v0.1.0
git push origin v0.1.0
```

The workflow cross-builds for:
- `x86_64-unknown-linux-gnu`
- `aarch64-unknown-linux-gnu`
- `x86_64-apple-darwin`
- `aarch64-apple-darwin`

Each tarball contains `gapsmith` binary, `data/` curation tables,
README, LICENSE, COMPARISON. A SHA256 checksum sits next to each
tarball. Generated release notes are populated from commit messages
between tags.

## crates.io

**Irreversible once published** — crate names are permanently claimed,
versions can only be yanked (not deleted).

### One-time setup

1. Create a crates.io account (GitHub OAuth): https://crates.io/
2. Generate an API token: https://crates.io/settings/tokens
3. Save locally:

   ```bash
   cargo login <your-token>
   ```

### Publish order (strict — dependents can't be published before their deps)

```bash
# Tier 0: no workspace deps
cargo publish -p gapsmith-core
cargo publish -p gapsmith-io        # depends on: core
cargo publish -p gapsmith-sbml      # depends on: core
cargo publish -p gapsmith-db        # depends on: core

# Wait ~60s between each — the index needs to propagate before the
# next cargo publish can resolve the just-uploaded crate.

# Tier 1
cargo publish -p gapsmith-align     # depends on: core
cargo publish -p gapsmith-medium    # depends on: core, db

# Tier 2
cargo publish -p gapsmith-draft     # depends on: core, db, sbml, io
cargo publish -p gapsmith-transport # depends on: core, db, align
cargo publish -p gapsmith-find      # depends on: core, db, align

# Tier 3
cargo publish -p gapsmith-fill      # depends on: core, db, draft

# Tier 4 (binary)
cargo publish -p gapsmith-cli       # depends on: all of the above
```

End state: `cargo install gapsmith-cli` installs the `gapsmith` binary.

### Before the first publish

- Verify every crate dry-runs cleanly:

  ```bash
  for crate in gapsmith-core gapsmith-io gapsmith-sbml gapsmith-db \
               gapsmith-align gapsmith-medium gapsmith-draft \
               gapsmith-transport gapsmith-find gapsmith-fill \
               gapsmith-cli; do
    cargo publish --dry-run -p "$crate" --allow-dirty
  done
  ```

  Only the first tier-0 crate can be fully verified offline; later ones
  need their upstream deps already on crates.io.

- Make sure all commits are pushed to `main` — crates.io embeds the git
  commit SHA in the published artefact.

- Check that no crate name is already taken by someone else:

  ```bash
  for c in core io sbml db align medium draft transport find fill cli; do
    curl -s -o /dev/null -w "gapsmith-$c: %{http_code}\n" \
      "https://crates.io/api/v1/crates/gapsmith-$c"
  done
  ```

  `404` means available; `200` means taken.

### Yank a broken version

```bash
cargo yank --vers 0.1.0 gapsmith-cli
```

Doesn't delete; flags the version so new resolvers don't pick it.

## Release checklist

- [ ] All 160+ tests pass (`cargo test --workspace`)
- [ ] ATP-cycle regression test green
  (`GAPSEQ_ROOT=/path/to/gapseq cargo test -p gapsmith-fill --test atp_cycle`)
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo fmt --check` clean
- [ ] Workspace version bumped in `Cargo.toml`
- [ ] CHANGELOG updated (TODO: add a CHANGELOG.md for the first release)
- [ ] Tag pushed → release workflow produces tarballs
- [ ] (Optional) `cargo publish` chain completes
