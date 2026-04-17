//! ATP-cycle regression test — the one invariant that defines a valid
//! gapsmith reaction database.
//!
//! Build a universal model (every gapseq-approved SEED reaction, no
//! exchanges, no biomass), set the objective to ATP hydrolysis
//! (`rxn00062_c0`: `atp + h2o → adp + h + pi`), run FBA. The optimum
//! must be ≈ 0. Anything above tolerance means some combination of
//! reactions generates ATP from nothing — almost always a reversibility
//! bug, occasionally a stoichiometry typo. The entire purpose of
//! `corrections_seed_reactionDB.tsv` and the gapseq-curated `seed_pwy.tbl`
//! is to keep this invariant true.
//!
//! Skipped when the large SEED table isn't accessible. To run:
//!
//! ```ignore
//! GAPSMITH_DATA_DIR=/path/to/dat cargo test -p gapsmith-fill --test atp_cycle -- --nocapture
//! ```

use gapsmith_core::Model;
use gapsmith_db::{load_seed_reactions, SeedRxnRow};
use gapsmith_fill::{fba, FbaOptions, SolveStatus};
use std::path::PathBuf;

/// Tolerance for "ATP production ≈ 0". HiGHS default primal feasibility
/// is 1e-7; we add margin because the split-flux encoding multiplies
/// slack.
const ATP_MAX_TOL: f64 = 1e-4;

/// Reaction id for ATP hydrolysis in SEED (atp + h2o → adp + h + pi).
const ATP_HYDROLYSIS_ID: &str = "rxn00062";

fn seed_rxns_path() -> Option<PathBuf> {
    let candidates = [
        std::env::var_os("GAPSMITH_DATA_DIR").map(PathBuf::from),
        std::env::var_os("GAPSEQ_ROOT")
            .map(|r| PathBuf::from(r).join("dat")),
        Some(PathBuf::from("/opt/gapseq/dat")),
        Some(PathBuf::from("dat")),
        Some(PathBuf::from("../dat")),
    ];
    candidates.into_iter().flatten().find_map(|d| {
        let p = d.join("seed_reactions_corrected.tsv");
        if p.is_file() {
            Some(p)
        } else {
            None
        }
    })
}

/// Build a universal model: every usable SEED reaction added with its
/// own reversibility-derived bounds, no exchanges, no biomass. This is
/// the closed-system configuration — any ATP produced must come from
/// internal stoichiometry alone.
fn build_universal_model(seed_rxns: &[SeedRxnRow]) -> Model {
    let mut m = Model::default();
    for row in seed_rxns.iter().filter(|r| r.gapseq_status.is_usable()) {
        gapsmith_draft::builder::add_seed_reaction(&mut m, row, None);
    }
    gapsmith_draft::builder::rebuild_s_matrix(&mut m);
    m
}

#[test]
fn no_free_atp_from_universal_model() {
    let path = match seed_rxns_path() {
        Some(p) => p,
        None => {
            eprintln!(
                "SKIP: seed_reactions_corrected.tsv not found. Set GAPSMITH_DATA_DIR \
                 or GAPSEQ_ROOT, or run from a checkout that has dat/ alongside."
            );
            return;
        }
    };

    eprintln!("loading SEED reactions from {}", path.display());
    let seed_rxns = load_seed_reactions(&path).expect("parse SEED reactions TSV");
    eprintln!("loaded {} reactions", seed_rxns.len());

    let mut model = build_universal_model(&seed_rxns);
    eprintln!(
        "universal model: {} reactions × {} metabolites",
        model.rxn_count(),
        model.met_count()
    );

    // Find ATP hydrolysis. We look for exactly `rxn00062_c0` first; if
    // absent (shouldn't happen), search by id prefix.
    let atp_idx = model
        .rxns
        .iter()
        .position(|r| r.id.as_str() == format!("{ATP_HYDROLYSIS_ID}_c0"))
        .or_else(|| {
            model.rxns.iter().position(|r| {
                r.id.as_str().starts_with(ATP_HYDROLYSIS_ID)
                    && r.id.as_str()[ATP_HYDROLYSIS_ID.len()..].starts_with('_')
            })
        })
        .expect("rxn00062 (ATP hydrolysis) missing from SEED DB");

    // Zero every objective coefficient, set ATP hydrolysis to +1, and
    // ensure its forward direction is allowed.
    for r in &mut model.rxns {
        r.obj_coef = 0.0;
    }
    model.rxns[atp_idx].obj_coef = 1.0;
    if model.rxns[atp_idx].ub < 1000.0 {
        model.rxns[atp_idx].ub = 1000.0;
    }

    let sol = fba(&model, &FbaOptions::default()).expect("FBA ran");
    eprintln!(
        "ATP-cycle FBA: status={:?} obj={:.6e}",
        sol.status, sol.objective
    );

    assert_eq!(
        sol.status,
        SolveStatus::Optimal,
        "universal-model FBA failed"
    );
    assert!(
        sol.objective.abs() <= ATP_MAX_TOL,
        "ATP-cycle invariant violated: FBA of universal model produces ATP = {:.6e} \
         without uptake (tol={:e}). Some reaction combination forms an ATP-producing \
         futile cycle — almost always a reversibility bug in seed_reactions_corrected.tsv \
         or a missing row in corrections_seed_reactionDB.tsv. Trace `sol.fluxes` to \
         identify the loop.",
        sol.objective,
        ATP_MAX_TOL
    );
}
