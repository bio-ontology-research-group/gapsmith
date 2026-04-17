# Vendored reference tables

This directory contains gapseq-authored curation tables bundled with
gapsmith. They are GPL-3 derivative works of upstream
[gapseq](https://github.com/jotech/gapseq)'s `dat/` directory.

## What's here

| File | Purpose |
|------|---------|
| `subex.tbl` | Substrate → exchange-reaction mapping (used by `find-transport`) |
| `tcdb_custom.tbl` | Custom TCDB substrate entries |
| `complex_subunit_dict.tsv` | Synonym dictionary for complex subunit detection |
| `exception.tbl` | EC/reaction exceptions that override default classification |
| `medium_prediction_rules.tsv` | Boolean-rule table for `medium` subcommand |
| `custom_pwy.tbl` | gapseq-authored custom pathway definitions |
| `seed_pwy.tbl` | gapseq's curated SEED pathway table |
| `corrections_seed_reactionDB.tsv` | gapseq's reversibility + stoichiometry corrections to SEED |
| `diffusion_mets.tsv` | Metabolites that diffuse freely |
| `ec_conflicts.tsv` | Known EC-number conflict resolution |
| `nutrients.tsv` | Nutrient classification table |
| `taxonomy.tbl` | Taxonomy string → Bacteria/Archaea mapping |
| `altec.csv` | Alternative EC number table |
| `importantCIreactions.list` | Important carbon-intake reactions |
| `reference_genomes_edited.tbl` | Curated reference genome list |
| `biomass/` | Per-taxonomy biomass template JSONs |
| `env/` | Environment constraint tables |
| `media/` | Growth-medium CSVs (MM_glu, LB, etc.) |

## What's NOT here (must come from `--data-dir`)

The large public reference tables:

- `seed_reactions_corrected.tsv`, `seed_reactions.tsv`, `seed_metabolites.tsv`,
  `seed_Enzyme_*.tsv`, `seed_transporter.tbl` — [ModelSEED](https://github.com/ModelSEED), public domain
- `mnxref_reac_xref.tsv`, `mnxref_seed.tsv`, `mnxref_*-other.tsv` —
  [MetaNetX](https://www.metanetx.org), CC0

License-restricted tables (require `--accept-license <db>` to load):

- `meta_pwy.tbl`, `meta_rea.tbl`, `meta_rea_pwy-gapseq.tbl`, `meta_genes.csv` —
  MetaCyc (requires user license acceptance)
- `kegg_pwy.tbl` — KEGG (academic / commercial license)
- `bigg_reactions.tbl`, `mnxref_bigg-other.tsv`, `SEED2VMH_translation*.csv`,
  `vmh_reactions.tsv` — BiGG Models (non-commercial academic)
- `brenda_ec.csv`, `brenda_ec_edited.csv` — BRENDA (CC-BY 4.0)
- `tcdb_substrates.tbl` — TCDB (free for academic)

## ATP-cycle integrity invariant

The vendored correction tables (`corrections_seed_reactionDB.tsv`,
`seed_pwy.tbl`) exist **solely** to keep the combined reaction database
free of ATP-producing futile cycles. The one invariant that defines a
valid gapsmith DB:

    max_v { v[atp_hydrolysis] : S·v = 0, lb ≤ v ≤ ub, all exchanges closed } ≈ 0

This is enforced by `atp_cycle::test_no_free_atp` in `gapsmith-fill`'s
test suite. **Any change to the reaction DB, corrections, or vendored
tables must keep this test green.** If it fails, trace the flux loop,
identify the reversibility or stoichiometry bug, fix it, and re-run.

## Updating

When upstream gapseq publishes new curation (via
[jotech/gapseq](https://github.com/jotech/gapseq)):

1. Diff upstream `dat/` against this directory.
2. Cherry-pick the relevant changes.
3. Run the ATP-cycle regression test (`cargo test -p gapsmith-fill atp_cycle`).
4. Commit with a reference to the upstream commit SHA.

Last sync: 2026-04-17 from gapseq master branch.
