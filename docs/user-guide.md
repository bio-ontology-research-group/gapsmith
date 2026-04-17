# gapseq-rs user guide

A Rust reimplementation of [gapseq](https://github.com/jotech/gapseq) —
informed prediction and analysis of bacterial metabolic pathways and
genome-scale networks. Single static binary, in-process LP solver, faster
on every stage of the pipeline (see benchmarks in the README).

This guide walks a new user through installing gapseq-rs, reconstructing a
genome-scale metabolic model end-to-end, and understanding each stage of
the output.

---

## 1. Install

### Binary prerequisites

gapseq-rs shells out to external aligners and a couple of HMM tools.
Install whichever aligner you want to use — you only need one:

| Tool | Required for | apt / bioconda |
|------|--------------|----------------|
| `blastp` + `makeblastdb` | `--aligner blastp` (default) | `ncbi-blast+` |
| `diamond` | `--aligner diamond` | `diamond` |
| `mmseqs` | `--aligner mmseqs2` | `mmseqs2` |
| `bedtools` | none (legacy; only needed for upstream `gapseq find` shell pipeline) | `bedtools` |

`--aligner precomputed` skips the aligner entirely and reads a user-
supplied TSV; see [Precomputed alignment mode](#precomputed-alignment-mode)
below.

### Reference data

gapseq-rs needs the standard gapseq `dat/` directory (SEED reactions +
metabolites, biomass templates, medium rules, pathway tables). Get it by
cloning the upstream repo:

```bash
git clone --depth 1 https://github.com/jotech/gapseq.git
export GAPSEQ_DATA_DIR=$PWD/gapseq/dat
```

Or point to it explicitly on every invocation with `--data-dir`.

### Reference sequences

Sync the genome-specific BLAST/diamond reference FASTA database from
Zenodo:

```bash
gapseq update-sequences -t Bacteria             # pinned release
gapseq update-sequences -t Bacteria -Z latest   # newest Zenodo version
gapseq update-sequences -c                       # check-only, no download
```

Without this, `find` / `find-transport` / `doall` will fail — the FASTA
files they align against live under `<data-dir>/seq/<Taxonomy>/`.

### Build

```bash
git clone <this-repo>
cd gapseq-rs
cargo build --release --workspace
# Binary is target/release/gapseq
```

The HiGHS LP solver is bundled via `good_lp` → `highs-sys` (CMake'd from
source on first build; ~30 s). An optional `cbc` feature adds CBC as a
fallback solver — requires `coinor-libcbc-dev` on the build host:

```bash
cargo build --release --workspace --features gapseq-fill/cbc
```

---

## 2. Quick start — reconstruct E. coli K-12 in one command

```bash
# Download an E. coli reference proteome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/all_assembly_versions/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz

# Reconstruct: find→find-transport→draft→medium (auto)→fill
gapseq --data-dir $GAPSEQ_DATA_DIR doall \
    GCF_000005845.2_ASM584v2_protein.faa.gz \
    -f ecoli_out/ \
    -A diamond \
    -b 200 \
    -k 0.01

# Inspect the result
gapseq fba ecoli_out/GCF_000005845.2_ASM584v2_protein-filled.gmod.cbor
```

Outputs in `ecoli_out/`:

- `<genome>-all-Reactions.tbl`, `<genome>-all-Pathways.tbl` — find results
- `<genome>-Transporter.tbl` — find-transport result
- `<genome>-draft.gmod.cbor`, `<genome>-draft.xml` — draft model
- `<genome>-medium.csv` — inferred growth medium
- `<genome>-filled.gmod.cbor`, `<genome>-filled.xml` — gap-filled model
- `<genome>-filled-added.tsv` — list of reactions added by each fill step

The filled CBOR or SBML loads directly in COBRApy, COBRAToolbox,
cobrar, or any SBML-aware tool.

---

## 3. Individual stages

### 3.1 `gapseq find` — pathway and reaction detection

Aligns the genome's proteins against gapseq's per-reaction reference
FASTAs, scores each pathway for completeness, and writes a
`*-Reactions.tbl` + `*-Pathways.tbl`.

```bash
# One pathway by MetaCyc id
gapseq find -p PWY-6587 genome.faa

# A category shorthand — expands to gapseq's internal hierarchy match
gapseq find -p amino genome.faa
#        ↑ one of: amino | nucl | cofactor | carbo | polyamine |
#                  fatty | energy | terpenoid | degradation |
#                  core | min | kegg | all

# Every active pathway in metacyc + custom (gapseq's default)
gapseq find -p all -l metacyc,custom genome.faa

# Alternative aligner
gapseq find -A diamond -p all genome.faa
gapseq find -A mmseqs2 -p all genome.faa
```

Key options:

- `-b` / `--bitcutoff` — BLAST bitscore cutoff for "high-evidence" hits (default 200).
- `-c` / `--coverage` — minimum query coverage (default 75%).
- `-a` / `--completeness-main` — pathway-completeness cutoff (default 80%).
- `-k` / `--completeness-hints` — relaxed cutoff when every key reaction is present (default 66%).
- `-n` / `--include-superpathways` — by default superpathways are filtered out (matches upstream gapseq).
- `-t` / `--taxonomy` — `Bacteria` (default) or `Archaea`.

The output `*-Pathways.tbl` is **byte-identical** to real gapseq's output
on tested cases (`-p PWY-6587` and `-p amino` on ecoli). See the
[porting-notes](porting-notes.md) for the handful of intentional
deviations.

### 3.2 `gapseq find-transport` — transporter detection

Iterates gapseq's transporter substrate list (`dat/subex.tbl`) and aligns
the genome against both TCDB and the curated SEED transporter FASTAs.

```bash
gapseq find-transport -A diamond -b 50 genome.faa
# → <stem>-Transporter.tbl
```

Row count and the TC-identifier set match real gapseq on the ecoli toy
genome (550 distinct TC ids, 1808 rows).

### 3.3 `gapseq draft` — build a draft model

Combines the `*-Reactions.tbl` + `*-Transporter.tbl` with the SEED
reaction DB into a metabolic model. Adds biomass, exchanges, diffusion
reactions, conditional transporters, dedups by stoichiometric hash.

```bash
gapseq draft \
    -r ecoli-all-Reactions.tbl \
    -t ecoli-Transporter.tbl \
    -b neg \
    -o drafts/
# → drafts/ecoli-draft.gmod.cbor   (native format)
#   drafts/ecoli-draft.xml         (SBML L3V1+FBC2+groups, validates
#                                   against libSBML with 0 errors)
```

The `-b` flag picks the biomass template:

- `auto` — read `gram=` from the Reactions.tbl header (gapseq defaults).
- `pos` / `neg` — Gram+ / Gram-.
- `archaea` — Archaea.
- `<path/to/custom.json>` — user biomass in the same JSON format as
  `dat/biomass/biomass_Gram_neg.json`.

### 3.4 `gapseq medium` — rule-based medium inference

Evaluates 82 boolean rules in `dat/medium_prediction_rules.tsv` against
the draft + Pathways.tbl to infer a growth medium.

```bash
gapseq medium \
    -m drafts/ecoli-draft.gmod.cbor \
    -p ecoli-all-Pathways.tbl \
    -o ecoli-medium.csv

# Manual overrides — force oxygen uptake 0 to predict anaerobic
gapseq medium ... -c "cpd00007:0;cpd00011:0"
```

The inferred medium CSV is compatible with `gapseq fill -n`. On the
ecoli toy genome the output is byte-identical to upstream gapseq's
`toy/ecoli-medium.csv`.

### 3.5 `gapseq fill` — iterative gap-filling

4-phase pFBA + KO-essentiality pipeline. Steps 1 + 2 + 2b are the
default "quick" path; Steps 3 + 4 run only with `--full-suite`.

```bash
gapseq fill \
    drafts/ecoli-draft.gmod.cbor \
    -n ecoli-medium.csv \
    -r ecoli-all-Reactions.tbl \
    -k 0.01 \
    -o filled/

# Full 4-phase suite (adds energy-source + fermentation screens;
# 10-30 minutes on a whole bacterium)
gapseq fill ... --full-suite

# Opt into the futile-cycle prune — drops candidate reactions that
# saturate at 0.99·max_flux with all boundaries closed. Expensive on
# large candidate pools.
gapseq fill ... --prune-futile
```

Phase breakdown:

1. **Step 1** — user medium + biomass target. Bulk of the fills.
2. **Step 2** — MM_glu + available carbon sources; iterate every
   biomass substrate, attach a sink per substrate, gap-fill if
   infeasible.
3. **Step 2b** — aerobic / anaerobic variant, depending on whether the
   user medium has O₂.
4. **Step 3** (`--full-suite`) — energy-source screen with ESP1–5
   pseudo-reactions (menaquinone / ubiquinone / NADH / ferredoxin /
   plastoquinone redox couples).
5. **Step 4** (`--full-suite`) — fermentation-product screen: maximise
   outflow on each exchange compound.

Output: `<genome>-filled.gmod.cbor` + `.xml` + `-added.tsv` (rxn id +
step number for every addition).

### 3.6 `gapseq fba` — FBA / pFBA on a model

Utility to solve FBA or parsimonious FBA on a CBOR/JSON model. Useful
for validating a draft before gap-filling, or for post-hoc analysis.

```bash
gapseq fba drafts/ecoli-draft.gmod.cbor
gapseq fba --pfba --min-growth 0.01 --top 20 filled/ecoli-filled.gmod.cbor
gapseq fba -r EX_cpd00027_e0 filled/ecoli-filled.gmod.cbor
```

Split-flux LP via `good_lp` with HiGHS backend. Solves the 2378×2953
ecoli draft in under 100 ms.

### 3.7 `gapseq adapt` — edit reactions + force growth

```bash
# Add specific SEED rxns
gapseq adapt -m model.gmod.cbor -a rxn00011,rxn00102

# Add every rxn in a MetaCyc pathway
gapseq adapt -m model.gmod.cbor -a "PWY-6587"

# Remove
gapseq adapt -m model.gmod.cbor -r rxn00011

# Force growth on D-glucose (runs gap-fill)
gapseq adapt -m model.gmod.cbor -w cpd00027:TRUE -b Reactions.tbl

# Force NO growth on glucose (closes uptake)
gapseq adapt -m model.gmod.cbor -w cpd00027:FALSE
```

### 3.8 `gapseq pan` — pan-draft from N drafts

```bash
# Glob a directory of drafts
gapseq pan -m drafts/ -t 0.06 -f pan_out/

# Explicit list
gapseq pan -m a-draft.gmod.cbor,b-draft.gmod.cbor,c-draft.gmod.cbor -t 0.5

# Only emit the binary presence/absence table
gapseq pan -m drafts/ -b -f pan_out/
```

`-t` sets the minimum across-draft reaction frequency for inclusion in
the pan-draft. Default 0.06 matches upstream gapseq.

### 3.9 `gapseq doall` — full pipeline in one command

```bash
gapseq doall genome.faa.gz -f out/ -A diamond -b 200 -k 0.01
# Chains: find -p all → find-transport → draft → medium (auto) → fill
```

Add `--full-suite` to include fill Steps 3/4. Use `-m <medium.csv>` to
skip the medium-inference step and supply your own medium.

---

## 4. Precomputed alignment mode

For many-genome batch runs, gapseq-rs can skip the per-genome aligner
call entirely and read a pre-computed alignment TSV:

```bash
# Pre-run BLAST / diamond / mmseqs2 yourself (or use gapseq batch-align)
diamond makedb --in reference.faa -d ref.dmnd
diamond blastp -q genome.faa -d ref.dmnd -o hits.tsv \
    --outfmt 6 qseqid sseqid pident qcovs evalue bitscore stitle sstart send

# Feed to find / find-transport
gapseq find -A precomputed -P hits.tsv -p all genome.faa
gapseq find-transport -A precomputed -P hits.tsv genome.faa
```

### `gapseq batch-align` — cluster N genomes + single alignment pass

When you're annotating many genomes with the same reference, run one
cluster + alignment over all of them and expand per-genome hits:

```bash
gapseq batch-align \
    -q reference.faa \
    -g genomes_dir/ \
    -o aligned/ \
    --aligner diamond \
    --cluster-identity 0.5 \
    --cluster-coverage 0.8
# → aligned/<genome_id>.tsv per genome in genomes_dir/
```

Amortises the alignment cost over many genomes — especially effective
when the proteomes share a lot of sequence homology.

---

## 5. Solver, format, and file-layout notes

### 5.1 Model format

- **Native**: `.gmod.cbor` — CBOR-encoded `Model` (fast load, compact).
- **JSON**: `.json` — human-readable, via `gapseq convert`.
- **SBML**: `.xml` — SBML L3V1 + FBC2 + groups. 0 libSBML errors on
  every emitted model. Loads cleanly in COBRApy / COBRAToolbox.

Conversion between native formats:

```bash
gapseq convert model.gmod.cbor model.json --pretty
gapseq convert model.json model.gmod.cbor
gapseq export-sbml model.gmod.cbor model.xml
```

### 5.2 LP solver

`good_lp` 1.15 + HiGHS (bundled). The HiGHS crate statically links its
C++ solver via CMake; no system HiGHS needed.

Optional CBC fallback: compile with `--features gapseq-fill/cbc` and
pass `--cbc-fallback` in the heuristic options. Requires
`coinor-libcbc-dev` + `coinor-libclp-dev` on the build host.

### 5.3 Sequence database layout

```
<seq-dir>/<Taxonomy>/
  ├── rev/<reaction>.fasta         ← curated reviewed sequences
  ├── unrev/<reaction>.fasta       ← unreviewed / broader
  ├── rxn/<rxn_id>.fasta           ← per-SEED-reaction FASTAs
  ├── tcdb.fasta                   ← transporter DB
  ├── transporter.fasta            ← curated transporter hits
  ├── user/<reaction>.fasta        ← your own references (optional)
  └── version_seqDB.json           ← Zenodo version stamp
```

`gapseq find` / `find-transport` resolve reference FASTAs in this
priority order: `user/` → `rxn/` → `rev/<EC>` → `unrev/<EC>` →
`rev/<md5(name)>` → `unrev/<md5(name)>`. Matching real gapseq's
`prepare_batch_alignments.R:150-234`.

---

## 6. Troubleshooting

### "external tool `makeblastdb` exited with status 1"

The tool was fed a gzipped FASTA. `gapseq doall` auto-decompresses
`.gz`; for standalone `find` / `find-transport` invocations, ungzip
first.

### "full candidate pool is infeasible"

The draft + every SEED reaction together can't produce biomass. Either:

- The biomass template asks for a metabolite no SEED reaction produces
  (peptidoglycan cpd15665 on some drafts — see porting-notes for the
  known cycle gap), or
- Your medium is missing a critical nutrient (run `gapseq medium
  --manual-flux` to check what the rule-based inference assigns).

Use `gapseq fba --pfba --min-growth 0.001` on the draft to narrow the
biomass-component that's blocking.

### pFBA heuristic "max iterations reached"

HiGHS exhausted the tolerance + pFBA-coefficient ladder. Options:

1. Widen the medium (add more carbon sources).
2. Lower `-k` (minimum growth rate).
3. Rebuild with `--features cbc` and opt-in to the CBC fallback.

### Draft has 0 reactions / stays empty

The upstream `*-Reactions.tbl` has `dbhit=` empty on every row. That
means none of the reference FASTAs were found for the pathways you
asked about. Verify `<seq-dir>/<Taxonomy>/` is populated
(`gapseq update-sequences -c`).

---

## 7. Performance tips

- Release build: `cargo build --release --workspace`. Debug builds are
  ~50× slower for the LP-heavy stages.
- Use diamond over BLASTp: `-A diamond`. Comparable accuracy, 5–20×
  faster on large proteomes.
- Parallelism: `-K <threads>` on every aligner-invoking subcommand.
  gapseq-rs is rayon-parallel in the LP-free stages (find, transport,
  medium) too.
- Batch mode: if you're annotating 10+ genomes, `gapseq batch-align`
  + `--aligner precomputed` amortises the alignment cost.

---

## 8. Where to look next

- [`porting-notes.md`](porting-notes.md) — every intentional deviation
  from upstream R gapseq.
- [`feature-matrix.md`](feature-matrix.md) — exhaustive feature list
  with R source → Rust module pointers.
- [README.md](../README.md) — status + benchmarks.
- Per-crate rustdoc: `cargo doc --workspace --no-deps --open`.
