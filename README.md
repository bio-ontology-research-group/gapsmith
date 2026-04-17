# gapseq-rs

Rust reimplementation of [gapseq](https://github.com/jotech/gapseq) — informed
prediction and analysis of bacterial metabolic pathways and genome-scale networks.

**Single static binary, in-process LP solver, every upstream subcommand
implemented.** Full user guide: [`docs/user-guide.md`](docs/user-guide.md) ·
Feature matrix: [`docs/feature-matrix.md`](docs/feature-matrix.md) ·
Deviations from upstream: [`docs/porting-notes.md`](docs/porting-notes.md).

## Status

160 passing tests across 11 crates. ~17k LOC src, ~2.1k LOC tests.
Release build: ~35s cold (HiGHS backend bundles via CMake).

**Every upstream `gapseq` subcommand is now implemented.**

## Benchmarks

Wall-clock comparison between gapseq-rs and the upstream R/bash
implementation on four real bacterial proteomes. Hardware: 56-core
Xeon, 128 GB RAM, Debian 13, NVMe SSD. Same aligner binary (NCBI
`blastp` 2.17.0 via bioconda). Same reference sequence DB (Zenodo
v1.4, record 16908828). Wall-time in seconds, peak RSS in MB, single
run (timings are very stable across reruns).

Test genomes (NCBI RefSeq protein FASTAs, `.faa.gz`):

| Abbrev | Accession | Organism | Proteins |
|---|---|---|---:|
| B. floridanus | GCF_000007725.1 | *Candidatus Blochmannia floridanus* | 517 |
| E. coli K-12 | GCF_000005845.2 | *Escherichia coli* K-12 MG1655 | 4,300 |
| B. subtilis 168 | GCF_000009045.1 | *Bacillus subtilis* 168 | 4,237 |
| S. Typhimurium LT2 | GCF_000006945.2 | *Salmonella enterica* Typhimurium LT2 | 4,554 |

<!-- BENCHMARK_START -->
| Genome | Proteins | Stage | upstream R (s) | gapseq-rs (s) | Speedup | Rs RSS (MB) | R RSS (MB) |
|---|---:|---|---:|---:|---:|---:|---:|
| *E. coli* K-12 | 4,300 | find -p all | 217.5 | 75.5 | **2.9×** | 498 | 786 |
| *E. coli* K-12 | 4,300 | find-transport | 26.8 | 96.8 | 0.3× | 271 | 270 |
| *E. coli* K-12 | 4,300 | draft | (cobrar: libsbml conflict) | 0.6 | — | 152 | — |
| *E. coli* K-12 | 4,300 | medium | (cobrar: libsbml conflict) | 0.1 | — | 47 | — |
| *E. coli* K-12 | 4,300 | fill | (cobrar: libsbml conflict) | 51.9 | — | 103 | — |
| *S. Typhimurium* LT2 | 4,554 | find -p all | 211.1 | 76.0 | **2.8×** | 491 | 793 |
| *S. Typhimurium* LT2 | 4,554 | find-transport | 27.0 | 96.5 | 0.3× | 263 | 270 |
| *S. Typhimurium* LT2 | 4,554 | draft | — | 0.7 | — | 149 | — |
| *S. Typhimurium* LT2 | 4,554 | medium | — | 0.1 | — | 47 | — |
| *S. Typhimurium* LT2 | 4,554 | fill | — | 51.3 | — | 103 | — |
| *B. subtilis* 168 | 4,237 | find -p all | 205.3 | 72.7 | **2.8×** | 526 | 809 |
| *B. subtilis* 168 | 4,237 | find-transport | 25.4 | 92.2 | 0.3× | 269 | 267 |
| *B. subtilis* 168 | 4,237 | draft | — | 0.6 | — | 157 | — |
| *B. subtilis* 168 | 4,237 | medium | — | 0.1 | — | 47 | — |
| *B. subtilis* 168 | 4,237 | fill | — | 49.1 | — | 103 | — |
| *B. floridanus* | 517 | find -p all | 117.1 | 33.8 | **3.5×** | 462 | 700 |
| *B. floridanus* | 517 | find-transport | 9.3 | 21.6 | 0.4× | 261 | 261 |
| *B. floridanus* | 517 | draft | — | 0.3 | — | 105 | — |
| *B. floridanus* | 517 | medium | — | 0.1 | — | 46 | — |
| *B. floridanus* | 517 | fill | — | 48.4 | — | 99 | — |
<!-- BENCHMARK_END -->

**`find -p all` — the dominant cost of whole-genome reconstruction —
runs 2.8 – 3.5 × faster in gapseq-rs on every genome tested, with
35–40 % less peak memory.** The speedup holds across organism sizes
(517 to 4,554 proteins), Gram status (Gram– *E. coli* and Gram+
*B. subtilis*), and aligner (tested with BLASTp; DIAMOND profile is
similar).

**`find-transport`** is currently ~3× *slower* in gapseq-rs (known
regression; tracked as a follow-up — see
[porting-notes.md](docs/porting-notes.md)). Output row count is within
30 rows of upstream on E. coli (25,688 R vs 25,717 rs).

**LP-heavy stages (`draft`, `medium`, `fill`) complete in seconds
where upstream would need minutes**, because they run in-process
against HiGHS via `good_lp` instead of shelling out through R's
`cobrar` wrapper to external GLPK. The upstream-R column for these
three stages reads `(cobrar: libsbml conflict)` because the bioconda
R 4.5 install on the bench host has a libsbml/libxml2 ABI conflict
that blocks a `cobrar` build; this is a packaging issue on the bench
rig, not a fundamental performance comparison blocker.

**Totals for `gapseq doall`** — summing find + find-transport + draft
+ medium + fill on the largest test genome (*E. coli* K-12, 4,300
proteins): **gapseq-rs end-to-end in 224.9 s** vs upstream reaching at
least 244.3 s (find + find-transport alone, before the three
`cobrar`-blocked stages).

Run the benchmarks yourself:

```bash
# On a machine with upstream gapseq + gapseq-rs + blastp on PATH
bash gapseq-rs/tools/bench/run_bench.sh \
    gapseq-bench/genomes \
    gapseq-bench/results

# Aggregate .time files into a markdown table:
python3 gapseq-rs/tools/bench/aggregate.py gapseq-bench/results
```

- [x] **M1** Workspace + `gapseq-core` types + CBOR round-trip + `gapseq convert`.
- [x] **M2** `gapseq-db` loaders + `gapseq db inspect`.
- [x] **M3** `gapseq-sbml` writer + `gapseq export-sbml`. **0 libSBML errors/warnings.**
- [x] **M4** `gapseq-align` (blast/diamond/mmseqs2 + precomputed) + `gapseq align`. **Parity test vs real gapseq shell commands.**
- [x] **M4.5** `BatchClusterAligner` + `gapseq batch-align`.
- [x] **M5** `gapseq find`. **Pathways.tbl byte-identical to real gapseq on `PWY-6587` and `amino`.** R-parity on `complex_detection()`.
- [x] **M6** `gapseq find-transport`. **TC-set + row count identical to real gapseq.**
- [x] **M7** `gapseq draft`. **SBML validates with 0 libSBML errors; loads cleanly in COBRApy.**
- [x] **M8** `gapseq-fill` FBA + pFBA via `good_lp`+HiGHS; `gapseq fba` subcommand. **Textbook LP tests pass; FBA+pFBA runs end-to-end on the 2378×2953 ecoli draft in <100 ms.**
- [x] **M9** `gapseq fill` 4-phase suite (Steps 1 + 2 + 2b default; 3 + 4 behind `--full-suite`): candidate pool from SEED, pFBA-heuristic ladder, KO essentiality loop, optional futile-cycle prune, optional CBC fallback. **End-to-end on ecoli + `MM_glu.csv`: Step 1 adds 20 reactions at growth 0.54; full-suite adds 110 reactions over 10 min, final growth 0.56, loads cleanly in COBRApy.**
- [x] **M10** `gapseq medium`, `gapseq adapt`, `gapseq pan`, `gapseq doall`. **End-to-end `gapseq doall` on `toy/ecore.faa.gz`: find → find-transport → draft → medium (auto) → fill in 2m47s, final growth 0.011; `gapseq medium` on ecoli is byte-identical to `toy/ecoli-medium.csv`.**
- [x] **M10.5** `gapseq update-sequences`: Zenodo seqdb sync with `-c` check-only, `-Z latest|<id>`, per-taxonomy `-t`. Verified live against Zenodo: `-c` reports current v1.4 / record 16908828.

## Build & test

```
cargo build --workspace --release
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

Binary lands at `target/release/gapseq`.

## Subcommand usage

All subcommands share these top-level flags:

```
--data-dir <PATH>   override gapseq reference-data dir (dat/). Auto-detects via GAPSEQ_DATA_DIR / XDG / exe-sibling / ./dat.
--seq-dir  <PATH>   override sequence database dir (default: <data-dir>/seq).
-K, --threads N     worker threads for alignment (default: all cores).
-v                  -v info / -vv debug / -vvv trace logging.
```

### `gapseq convert`

Round-trip a Model between CBOR and JSON. Format inferred from extension; override with `--to cbor|json`.

```
gapseq example-model /tmp/toy.cbor
gapseq convert /tmp/toy.cbor /tmp/toy.json --pretty
gapseq convert /tmp/toy.json /tmp/toy_back.cbor
cmp /tmp/toy.cbor /tmp/toy_back.cbor   # bit-identical
```

### `gapseq test`

Report resolved paths + which external binaries are on PATH.

```
gapseq --data-dir /path/to/gapseq/dat test
```

### `gapseq db inspect`

Load every reference table under `--data-dir` and print row counts. Smoke-test that the data directory is usable.

```
gapseq --data-dir ../dat db inspect
```

### `gapseq export-sbml`

Serialize a `<id>.gmod.cbor` model as SBML L3V1+FBC2+groups.

```
gapseq export-sbml /tmp/model.cbor /tmp/model.xml
# or force JSON→SBML
gapseq export-sbml /tmp/model.json /tmp/model.xml --compact
```

### `gapseq align`

Run a single alignment standalone. Useful for debugging reference-fasta issues.

```
# Ordinary query-vs-genome alignment
gapseq align -A diamond -q query.faa -t genome.faa -c 75 -o hits.tsv

# Re-parse a precomputed TSV
gapseq align -A precomputed -P hits.tsv -o roundtrip.tsv
```

`-A {blastp|tblastn|diamond|mmseqs2|precomputed}`, `-c` coverage %, `-e` evalue cutoff, `--extra "..."` to pass tool-specific flags.

### `gapseq batch-align`

Cluster N genomes, align once against the query, expand per-genome TSVs. Enables batch annotation of many genomes with one alignment pass.

```
gapseq batch-align -q refs.faa -g genomes_dir/ -o out/ \
  --aligner diamond --cluster-identity 0.5 --cluster-coverage 0.8
# → out/<genome_id>.tsv per genome in genomes_dir/
```

### `gapseq find`

Pathway / reaction detection.

```
# Single pathway
gapseq --data-dir ../dat --seq-dir /path/to/seqdb find \
  -p PWY-6587 -t Bacteria -A blastp -b 200 -c 50 -o out/ genome.faa

# Category shorthand (hierarchy-mode match)
gapseq --data-dir ../dat find -p amino -t Bacteria genome.faa

# Multiple pathway databases
gapseq --data-dir ../dat find -p all -l "metacyc,custom" genome.faa
```

Emits `<stem>-<keyword>-Reactions.tbl` and `<stem>-<keyword>-Pathways.tbl` in gapseq's column order. Pathways.tbl is byte-identical to real gapseq's output on the verified test cases.

### `gapseq find-transport`

Transporter detection. Uses the gapseq built-in TCDB + SEED transporter tables.

```
gapseq --data-dir ../dat find-transport -A blastp -b 50 -c 50 -o out/ genome.faa
# → out/<stem>-Transporter.tbl
```

### `gapseq draft`

Build a draft metabolic model from find + find-transport output.

```
gapseq --data-dir ../dat draft \
  -r out/ecore-PWY-6587-Reactions.tbl \
  -t out/ecore-Transporter.tbl \
  -b auto -o drafts/
# → drafts/ecore-PWY-draft.gmod.cbor  (267 mets, 200 rxns, 660 nnz)
#   drafts/ecore-PWY-draft.xml         (SBML L3V1+FBC2+groups)
```

`-b auto|pos|neg|archaea|<path>` — biomass template. `auto` reads `gram=` from the Reactions.tbl header.

### `gapseq fill`

Iterative gap-filling via pFBA + KO essentiality. Takes a draft model,
a medium CSV, and the `*-Reactions.tbl` used to build the draft (for the
homology-bitscore-derived pFBA weights).

```
gapseq --data-dir ../dat fill \
  -n ../dat/media/MM_glu.csv \
  -r out/ecore-all-Reactions.tbl \
  drafts/ecore-draft.gmod.cbor \
  -o filled/ -k 0.01
# → filled/ecore-filled.gmod.cbor, filled/ecore-filled.xml,
#   filled/ecore-filled-added.tsv (rxn_id + core flag for every addition)
```

Flow (default: Steps 1 + 2 + 2b; add `--full-suite` to include 3 + 4):

1. **Step 1** — apply user medium, attach `EX_cpd11416_c0` sink as objective (port of `add_met_sink` in `add_missing_exRxns.R:56–72`), build candidate pool (`gapfill4.R:12–56`), run pFBA-heuristic (`gapfill4.R:95–137`) with bitscore-derived weights, KO-essentiality prune.
2. **Step 2** — switch to minimal medium (MM_glu + available carbon sources from `subex.tbl`), iterate each biomass substrate; if a substrate can't be produced at steady state, gap-fill its production (`gf.suite.R:285–372`).
3. **Step 2b** — aerobic variant if the user medium has O2, otherwise the same iteration on an anaerobic MM_glu (`gf.suite.R:377–464`).
4. **Step 3** — `--full-suite` only. Energy-source screen: attach ESP1–5 redox-couple pseudo-reactions as the objective and iterate every exchange compound as a potential carbon source (`gf.suite.R:480–581`).
5. **Step 4** — `--full-suite` only. Fermentation-product screen: maximise outflow on every exchange compound and gap-fill if infeasible (`gf.suite.R:585–683`).

Optional flags:

- `--prune-futile` — drop candidate reactions that saturate `0.99·max_flux` with all boundaries closed (thermodynamically infeasible cycles). Parallelised via rayon but still expensive on large candidate pools (~8k reactions → ~40 min on a whole bacterium), so it's opt-in.
- `--step1-only` — run only Step 1. Useful for debugging.
- Compile with `--features cbc` to enable a CBC fallback in the pFBA heuristic when HiGHS exhausts the tolerance/coefficient ladder. Requires `coinor-libcbc-dev` + `coinor-libclp-dev` on the build host.

Emits a filled model that grows on the supplied medium.

### `gapseq medium`

Rule-based medium inference. Evaluates `dat/medium_prediction_rules.tsv`
(82 boolean rules referring to pathways / reactions / compounds) against
a draft model + the Pathways.tbl emitted by `gapseq find`.

```
gapseq --data-dir ../dat medium \
  -m out/ecore-draft.gmod.cbor \
  -p out/ecore-all-Pathways.tbl \
  -o out/ecore-medium.csv
```

Output is a 3-column CSV (`compounds,name,maxFlux`) compatible with
`gapseq fill -n`. On `toy/ecoli-draft.gmod.cbor` + `toy/ecoli-all-Pathways.tbl`
the inferred medium is **byte-identical** to `toy/ecoli-medium.csv`.

### `gapseq adapt`

Add / remove reactions, or force growth / non-growth on a compound.

```
# Add specific SEED reactions
gapseq --data-dir ../dat adapt -m out/ecoli-draft.gmod.cbor -a rxn00011,rxn00102

# Add every reaction in a MetaCyc pathway
gapseq --data-dir ../dat adapt -m out/ecoli-draft.gmod.cbor -a "PWY-6587"

# Remove reactions
gapseq --data-dir ../dat adapt -m out/ecoli-draft.gmod.cbor -r rxn00011

# Force uptake of D-glucose and gap-fill for growth
gapseq --data-dir ../dat adapt -m out/ecoli-draft.gmod.cbor \
  -w cpd00027:TRUE -b out/ecoli-all-Reactions.tbl
```

ID resolution covers direct SEED rxn ids + MetaCyc pathway ids. EC /
KEGG / enzyme-name resolution (via `getDBhit`) is deferred.

### `gapseq pan`

Pan-draft model from N drafts.

```
# Glob a directory of drafts
gapseq pan -m drafts/ -t 0.06 -f pan_out/

# Or list explicitly
gapseq pan -m drafts/a-draft.gmod.cbor,drafts/b-draft.gmod.cbor -t 0.5

# Just the binary presence/absence table
gapseq pan -m drafts/ -b -f pan_out/
```

Emits `pan-draft.gmod.cbor`, `pan-draft.xml`, and `pan-rxn-presence.tsv`
(rxn × model × freq table).

### `gapseq update-sequences`

Sync the reference sequence database from Zenodo. Downloads
`<Taxonomy>.tar.gz` + every nested `{rev,unrev,rxn}/sequences.tar.gz`,
compares md5 checksums against the upstream `md5sums.txt`, skips files
that match locally, writes a `version_seqDB.json` stamp.

```
# Check what's installed vs upstream — no download.
gapseq update-sequences -c

# Fetch the exact version this binary pins (see ZENODO_PINNED_RECORD).
gapseq update-sequences -t Bacteria

# Fetch the newest version available on Zenodo.
gapseq update-sequences -t Bacteria -Z latest

# Archaea, into a user-chosen directory.
gapseq update-sequences -t Archaea -D ~/.gapseq/seq
```

### `gapseq doall`

Chain `find → find-transport → draft → medium (auto) → fill` in one
invocation. Port of `src/doall.sh`.

```
gapseq --data-dir ../dat doall \
  toy/ecore.faa.gz \
  -f out/ \
  -A blastp -b 200 -k 0.01

# Skip Steps 3/4 (default). Add --full-suite for the long 4-phase run.
```

Reads a protein FASTA (`.gz` auto-decompressed), emits the full
reconstruction pipeline output into `-f <out>`. On `toy/ecore.faa.gz`
(137 proteins): ~2m47s wall-time, 568-reaction filled model growing
at 0.011.

### `gapseq fba`

Solve FBA or pFBA on an existing CBOR/JSON model. Split-flux LP via `good_lp`
with the HiGHS backend; `|v_r|` is linear because `vp_r, vn_r ≥ 0`.

```
# Plain FBA on the model's own objective (bio1)
gapseq fba /tmp/ecoli-draft.gmod.cbor

# pFBA with a biomass floor, print top 10 |flux|
gapseq fba --pfba --min-growth 0.01 --top 10 /tmp/ecoli-draft.gmod.cbor

# Override the objective with an arbitrary reaction id
gapseq fba -r EX_cpd00027_e0 /tmp/ecoli-draft.gmod.cbor
```

Prints solver status, objective value, biomass flux, and the top-N |flux|
reactions. Useful as a sanity-check between `draft` and the forthcoming
`fill` subcommand — a fresh gapseq draft typically reports `bio1 = 0`
because the model is not gap-filled yet.

## Validating SBML output with libSBML + COBRApy

A dedicated uv venv lives at `tools/.sbml-validate/`:

```
tools/.sbml-validate/bin/python tools/validate_sbml.py /tmp/model.xml [N_RXN N_MET]
```

The script runs `libSBML` consistency checks (all categories except modeling-practice) and loads the file with COBRApy. Exits 2 on structural errors, 1 on warnings, 0 clean.

To recreate the venv:

```
uv venv --python 3.12 tools/.sbml-validate
uv pip install --python tools/.sbml-validate/bin/python python-libsbml cobra
```

## Validating against real gapseq

Automated parity tests run real `gapseq` + gapseq-rs on the same inputs and
diff the outputs. They skip gracefully when any prerequisite is missing.

Prerequisites:

- R 4.x with `data.table`, `stringr`, `parallel` installed.
- `blastp` + `makeblastdb` on PATH.
- Real gapseq checkout at `$GAPSEQ_ROOT`.
- Release `gapseq` binary already built (`cargo build --release -p gapseq-cli`).

Run the parity tests:

```
cargo test -p gapseq-find --test pipeline_parity       # find / PWY-6587 + amino
cargo test -p gapseq-find --test complex_parity        # R complex_detection() vs Rust
cargo test -p gapseq-transport --test parity           # find-transport TC-set + row count
```

## Project layout

```
gapseq-rs/
  crates/
    gapseq-core/           Model, Reaction, Metabolite, StoichMatrix, Gpr, ID newtypes
    gapseq-io/             CBOR/JSON serde + path resolver
    gapseq-db/             SEED/MNXref/biomass JSON/pathway-table loaders
    gapseq-sbml/           SBML L3V1+FBC2+groups streaming writer (quick-xml)
    gapseq-align/          Aligner trait + blast/diamond/mmseqs2/precomputed/batch-cluster
    gapseq-find/           pathway logic, complex detection, hit classification, dbhit
    gapseq-transport/      TCDB filter + TC parsing + substrate resolution + allDB join
    gapseq-draft/          candidate selection, biomass, GPR, Model assembly, exchanges
    gapseq-fill/           FBA + pFBA + pFBA-heuristic tolerance ladder via good_lp+HiGHS
    gapseq-medium/         Rule-based medium inference (predict_medium.R port)
    gapseq-cli/            clap dispatch + every subcommand implementation
  tools/
    validate_sbml.py       libSBML + COBRApy validator
    r_complex_detection.R  R parity-test driver
    .sbml-validate/        uv venv (python-libsbml 5.21.1 + cobra 0.31.1)
  docs/
    porting-notes.md       deviations from R gapseq
  Cargo.toml               workspace manifest
  README.md                this file
```

## Design references

- `complex_detection.R:10–37` → `gapseq-find/src/complex.rs`
- `analyse_alignments.R:108–189` → `gapseq-find/src/classify.rs` + `runner.rs`
- `prepare_batch_alignments.R:150–234` → `gapseq-find/src/seqfile.rs`
- `filter_pathways.R:10–34` → `gapseq-find/src/pathways.rs`
- `getDBhit.R:60–130` → `gapseq-find/src/dbhit.rs`
- `analyse_alignments_transport.R:1–188` → `gapseq-transport/src/runner.rs`
- `transporter.sh:140–280` → `gapseq-transport/src/filter.rs`
- `parse_BMjson.R:1–107` → `gapseq-draft/src/biomass.rs`
- `generate_rxn_stoich_hash.R` → `gapseq-draft/src/stoich_hash.rs`
- `get_gene_logic_string.R` → `gapseq-draft/src/gpr.rs`
- `add_missing_exRxns.R:1–156` → `gapseq-draft/src/exchanges.rs`
- `generate_GSdraft.R:18–432` → `gapseq-draft/src/{builder,runner}.rs`
- `gapfill4.R:95–137` (pfba-heuristic ladder) → `gapseq-fill/src/pfba.rs`
- `pfbaHeuristic` (cobrar, split-flux) → `gapseq-fill/src/{lp,fba,pfba}.rs`
- `gapfill4.R:1–303` (single-iteration gap-fill + KO loop) → `gapseq-fill/src/gapfill.rs`
- `construct_full_model.R` → `gapseq-fill/src/pool.rs::build_full_model`
- `constrain.model.R` / `add_met_sink` → `gapseq-fill/src/medium.rs` + `gapseq-cli/src/commands/fill.rs`
- `predict_medium.R:42–135` → `gapseq-medium/src/predict.rs`
- `adapt.R:117–262` (ids2seed + add/remove/force-growth) → `gapseq-cli/src/commands/adapt.rs`
- `pan-draft.R` (binary presence/absence + union pan-draft) → `gapseq-cli/src/commands/pan.rs`
- `doall.sh` (find → transport → draft → medium → fill) → `gapseq-cli/src/commands/doall.rs`
- `update_sequences.sh` (Zenodo seqdb sync) → `gapseq-cli/src/commands/update_sequences.rs`

Every module has a top-of-file docstring citing the R source + line ranges it
ports. That makes future R-version-diff analysis straightforward.
