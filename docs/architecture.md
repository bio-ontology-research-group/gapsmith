# Architecture

A bird's-eye view of how gapseq-rs is organised, how data flows through
it, and what each crate is responsible for.

---

## 1. Crate dependency graph

```
             gapseq-core            ← Types only, no I/O / DB / solver
              ↑  ↑  ↑
              │  │  │
    ┌─────────┘  │  └────────────────┐
    │            │                   │
  gapseq-io   gapseq-db           gapseq-sbml
    │            │                   ↑
    └─┬──────────┘                   │
      │                              │
      ▼                              │
  gapseq-align ──┐                   │
      │          │                   │
      ▼          │                   │
  gapseq-find    │                   │
      │          │                   │
      ├─▶ gapseq-transport           │
      │                              │
      ▼                              │
  gapseq-draft ──────────────────────┤
      │                              │
      ▼                              │
  gapseq-fill ◀─── gapseq-medium    │
      │                              │
      └──────────────▶ gapseq-cli ◀──┘
```

Arrows point from "depended on" toward "depends on". `gapseq-core` is
the leaf: every other crate builds on its `Model` / `Reaction` /
`Metabolite` / `StoichMatrix` types.

---

## 2. Data flow for `gapseq doall`

```text
  genome.faa.gz
        │
        ▼
  ┌──────────────┐                      ┌──────────────────────┐
  │ gapseq find  │ ──── aligner (shell) │  dat/seq/<tax>/      │
  │ -p all       │          ─ BLAST/    │  rev/ unrev/ rxn/    │
  └──────┬───────┘            diamond/  └──────────────────────┘
         │                    mmseqs2
         │   *-Reactions.tbl  + *-Pathways.tbl
         ▼
  ┌──────────────────────┐
  │ gapseq find-transport│ ─── same aligner, subex.tbl + tcdb.fasta
  └──────┬───────────────┘
         │   *-Transporter.tbl
         ▼
  ┌────────────────────┐
  │ gapseq draft       │ ─── seed_reactions_corrected.tsv, biomass JSON
  └──────┬─────────────┘
         │   *-draft.gmod.cbor  (CBOR) + *-draft.xml  (SBML)
         ▼
  ┌──────────────────┐
  │ gapseq medium    │ ─── medium_prediction_rules.tsv + *-Pathways.tbl
  │ (auto)           │
  └──────┬───────────┘
         │   *-medium.csv
         ▼
  ┌────────────────────┐
  │ gapseq fill        │ ─── HiGHS LP (in-process via good_lp)
  │ 4-phase suite      │
  └──────┬─────────────┘
         │   *-filled.gmod.cbor  +  *-filled.xml  +  *-filled-added.tsv
         ▼
  (COBRApy, COBRAToolbox, cobrar, …)
```

Every arrow between stages is plain filesystem I/O. No daemon, no IPC.
Each subcommand is individually re-runnable; if you change the medium
you only re-run `fill`.

---

## 3. Model representation

[`gapseq_core::Model`] is the single source of truth:

```rust
pub struct Model {
    pub annot: ModelAnnot,           // provenance metadata
    pub compartments: Vec<Compartment>,
    pub mets: Vec<Metabolite>,
    pub rxns: Vec<Reaction>,
    pub genes: Vec<GeneId>,
    pub s: StoichMatrix,              // sprs CsMat in CSC
}
```

Serialised via serde as:

- **CBOR** (`.gmod.cbor`) — native, compact, fast load. The gapseq-rs
  internal format; replaces upstream's `.RDS`.
- **JSON** (`.json`) — human-readable; `gapseq convert` converts both
  ways. Same semantic content as CBOR.
- **SBML** (`.xml`) — Level 3 Version 1 + FBC 2 + groups 1. Written by
  [`gapseq_sbml::write_sbml`]; loads cleanly in COBRApy / COBRAToolbox
  / cobrar.

Metabolite ids use the `cpd00001_c0` convention (compartment baked into
the id) for SBML SId compliance. Reaction ids likewise — `rxn00001_c0`,
`EX_cpd00001_e0`, `bio1`, etc.

---

## 4. LP plumbing

The gap-filler does a lot of LPs — dozens per Step 1 / 2 / 2b / 3 / 4
run, sometimes thousands in a `--full-suite`. Key design decisions:

### Split-flux encoding

Every reaction `r` becomes two variables `vp_r, vn_r ≥ 0`. The net
flux is `v_r = vp_r − vn_r`. This makes `|v_r| = vp_r + vn_r` linear —
essential for pFBA.

### Bound translation

| model `[lb, ub]`           | `vp_r` upper | `vn_r` upper |
|----------------------------|---|---|
| `lb ≥ 0, ub ≥ 0`           | `ub`           | `0`          |
| `lb ≤ 0, ub ≥ 0`           | `ub`           | `−lb`        |
| `lb ≤ 0, ub ≤ 0`           | `0`            | `−lb`        |

Strict-direction reactions additionally get `vp_r ≥ lb` (forward-forced)
or `vn_r ≥ −ub` (backward-forced).

### Mass balance

Walk the CSC matrix once column-by-column (`O(nnz)` total). For each
non-zero `(row, col, coef)`, accumulate `coef · (vp[col] − vn[col])`
onto row `row`'s expression. The naive per-row scan would be
`O(m · n)` — ~45 B iterations on a 3k×5k model. The `O(nnz)` walk
solves it in under 100 ms on ecoli.

### Solver

`good_lp` 1.15 with HiGHS backend. HiGHS is statically linked via
CMake at build time — no runtime HiGHS dependency. Optional CBC
fallback behind `--features cbc`: `pfba_heuristic` retries once with
CBC after HiGHS exhausts the tolerance / coefficient ladder.

---

## 5. Alignment layer

`gapseq-align` exposes a single [`Aligner`] trait:

```rust
pub trait Aligner {
    fn run(&self, query: &Path, targets: &[&Path], opts: &AlignOpts)
        -> Result<Vec<Hit>, AlignError>;
}
```

Implementors — one struct each — shell out to the matching binary and
normalise its output to a [`Hit`] with the gapseq-canonical column
set (qseqid / pident / evalue / bitscore / qcovs / stitle / sstart /
send).

- [`BlastpAligner`] builds an indexed BLAST DB in a temp dir per call.
- [`DiamondAligner`] mirrors gapseq's `--more-sensitive` default.
- [`Mmseqs2Aligner`] uses the explicit `createdb → search →
  convertalis` 4-command pipeline (not `easy-search` — see
  [porting-notes](porting-notes.md) for why).
- [`PrecomputedTsvAligner`] skips all shelling out and reads a
  user-supplied TSV.
- [`BatchClusterAligner`] is the Rust-only feature: mmseqs2-clusters
  N input genomes, runs one alignment, expands cluster membership
  back to per-genome hits.

---

## 6. Candidate-pool / gap-fill

The 4-phase gap-fill suite (`gapseq-fill::run_suite`) is the heaviest
computation in the pipeline. Sketch of one fill call:

```
  draft model  +  medium CSV  +  Reactions.tbl (bitscore weights)
                         │
                         ▼
          ┌──────────────────────────────┐
          │ apply_medium(draft)          │
          │ attach EX_<target>_c0 obj    │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ build_full_model             │  ← every approved SEED rxn
          │   draft + candidates         │    not already present,
          │   deduped by stoich hash     │    sorted by core status
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ (optional) detect_futile_    │  ← parallel LP probe,
          │   cycles  +  drop            │    `--prune-futile`
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ pfba_heuristic               │  ← tolerance ladder,
          │   min_growth ≥ k             │    `1e-6 → 1e-9`
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ extract utilized candidates  │
          │ add to draft                 │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ KO essentiality loop         │  ← zero bounds, recheck FBA,
          │   serial, core-first         │    drop non-essential
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ Step 2: biomass components   │  ← per substrate, MM_glu
          ├──────────────────────────────┤
          │ Step 2b: aerobic/anaerobic   │
          ├──────────────────────────────┤
          │ Step 3: ESP1-5 energy screen │  (`--full-suite`)
          ├──────────────────────────────┤
          │ Step 4: fermentation screen  │  (`--full-suite`)
          └──────────────────────────────┘
                         │
                         ▼
                    filled model
```

---

## 7. Testing surface

Every algorithmic component has unit tests (160 total, `cargo test
--workspace`). On top, three integration layers:

1. **Shell-parity tests** (`crates/gapseq-align/tests/*_parity.rs`)
   — run the real BLAST / diamond / mmseqs2 binaries and diff against
   our wrappers' output.
2. **R-parity tests** (`crates/gapseq-find/tests/complex_parity.rs`,
   `pipeline_parity.rs`, `crates/gapseq-transport/tests/parity.rs`)
   — run the actual R gapseq via `Rscript` and diff column-by-column.
3. **SBML validator** (`tools/validate_sbml.py`) — uv-managed venv
   with `python-libsbml` + `cobra`; runs libSBML consistency checks +
   COBRApy round-trip on emitted SBML. 0 errors on every model.

---

## 8. File-system layout

```
gapseq-rs/
  crates/
    gapseq-core/               # types
    gapseq-io/                 # CBOR/JSON + path resolver
    gapseq-db/                 # dat/*.tsv loaders
    gapseq-sbml/               # SBML writer
    gapseq-align/              # aligner trait + 5 backends
    gapseq-find/               # pathway + reaction detection
    gapseq-transport/          # transporter detection
    gapseq-draft/              # draft model builder
    gapseq-medium/             # rule-based medium inference
    gapseq-fill/               # FBA / pFBA / gap-filler / suite
    gapseq-cli/                # clap dispatch + every subcommand
  tools/
    bench/                     # R-vs-rs benchmark runner
    validate_sbml.py           # libSBML + COBRApy validator
    r_complex_detection.R      # R parity harness
    .sbml-validate/            # uv venv (python-libsbml + cobra)
  docs/
    user-guide.md              # end-user walk-through
    cli-reference.md           # exhaustive flag list
    feature-matrix.md          # R source → Rust module pointers
    porting-notes.md           # intentional deviations from upstream
    architecture.md            # this file
  Cargo.toml                   # workspace manifest
  README.md                    # status + benchmarks
```
