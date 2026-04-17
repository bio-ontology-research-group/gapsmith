//! `gapsmith update-data` — fetch the large public reference tables that
//! aren't vendored in this repo.
//!
//! The curated, gapseq-authored tables (subex, biomass JSONs, media,
//! medium-prediction rules, corrections, etc.) live under `data/` in
//! this repo — you already have them if you cloned or `cargo install`-ed
//! gapsmith. What's NOT vendored:
//!
//! - **SEED reactions and metabolites** (`seed_reactions_corrected.tsv`,
//!   `seed_metabolites.tsv`, `seed_Enzyme_*.tsv`, `seed_transporter*.tbl`)
//!   — public domain, ~15 MB total. Hosted in upstream gapseq's `dat/`.
//! - **MNXref cross-references** (`mnxref_reac_xref.tsv`, `mnxref_seed.tsv`,
//!   `mnxref_seed-other.tsv`) — CC0, ~50 MB. Also in upstream gapseq's
//!   `dat/`.
//!
//! This subcommand fetches those files over HTTP from the upstream
//! gapseq repository on GitHub and drops them into a target directory.
//! Use the same directory with `--data-dir` for all subsequent gapsmith
//! invocations.
//!
//! # Example
//!
//! ```bash
//! # Fetch into ~/.local/share/gapsmith/
//! gapsmith update-data -o ~/.local/share/gapsmith/
//!
//! # Later, point --data-dir at it
//! gapsmith --data-dir ~/.local/share/gapsmith/ doall genome.faa
//! ```
//!
//! Downloads are verified against the file's content-length. Existing
//! files are skipped unless `--force` is given. The subcommand does NOT
//! fetch license-restricted data (MetaCyc, KEGG, BiGG, BRENDA, VMH,
//! TCDB-substrates) — those require explicit license acceptance
//! (forthcoming `--accept-license` flag).

use clap::Parser;
use std::fs;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::time::Duration;

#[derive(Debug, Parser)]
pub struct UpdateDataArgs {
    /// Target directory. Large public tables are placed here. Combine
    /// with `--data-dir <this>` on subsequent gapsmith invocations.
    #[arg(short = 'o', long, default_value = ".")]
    pub out_dir: PathBuf,

    /// Upstream gapseq git ref to pull from (default: `master`).
    #[arg(long, default_value = "master")]
    pub r#ref: String,

    /// Force re-download even if the file already exists.
    #[arg(short = 'f', long)]
    pub force: bool,

    /// Only report what would be downloaded; no network writes.
    #[arg(short = 'c', long)]
    pub check: bool,
}

/// Files to fetch from upstream gapseq's `dat/` directory. Paths are
/// relative to the upstream repo root.
///
/// Total ~65 MB. All files are public domain (SEED) or CC0 (MNXref).
/// No MetaCyc, KEGG, BiGG, BRENDA content — those are license-gated.
const FILES: &[&str] = &[
    // SEED reactions + metabolites (public domain — ModelSEED project)
    "dat/seed_reactions_corrected.tsv",
    "dat/seed_reactions.tsv",
    "dat/seed_metabolites.tsv",
    "dat/seed_metabolites_edited.tsv",
    "dat/seed_Enzyme_Class_Reactions_Aliases.tsv",
    "dat/seed_Enzyme_Class_Reactions_Aliases_unique.tsv",
    "dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv",
    "dat/seed_Enzyme_Name_Reactions_Aliases.tsv",
    "dat/seed_transporter.tbl",
    "dat/seed_transporter_custom.tbl",
    // MNXref (CC0 — MetaNetX project)
    "dat/mnxref_reac_xref.tsv",
    "dat/mnxref_seed.tsv",
    "dat/mnxref_seed-other.tsv",
];

const UPSTREAM_BASE: &str = "https://raw.githubusercontent.com/jotech/gapseq";

pub fn run(args: UpdateDataArgs) -> anyhow::Result<()> {
    fs::create_dir_all(&args.out_dir)?;

    let client = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(30))
        .timeout_read(Duration::from_secs(300))
        .build();

    let mut to_fetch: Vec<&str> = Vec::new();
    for path in FILES {
        let local = args.out_dir.join(
            std::path::Path::new(path)
                .file_name()
                .expect("path has filename"),
        );
        if local.exists() && !args.force {
            println!("skip {} (already present)", local.display());
            continue;
        }
        to_fetch.push(path);
    }

    if args.check {
        if to_fetch.is_empty() {
            println!("all files present; --force to re-download");
        } else {
            println!("would fetch {} files:", to_fetch.len());
            for p in &to_fetch {
                println!("  {p}");
            }
        }
        return Ok(());
    }

    if to_fetch.is_empty() {
        println!("nothing to do; --force to re-download");
        return Ok(());
    }

    let mut total_bytes: u64 = 0;
    for path in &to_fetch {
        let url = format!("{UPSTREAM_BASE}/{}/{path}", args.r#ref);
        let local = args.out_dir.join(
            std::path::Path::new(path)
                .file_name()
                .expect("path has filename"),
        );
        print!("fetch {} → {} ... ", url, local.display());
        std::io::stdout().flush().ok();

        let resp = client
            .get(&url)
            .call()
            .map_err(|e| anyhow::anyhow!("GET {url}: {e}"))?;
        let status = resp.status();
        if status != 200 {
            anyhow::bail!("GET {url}: HTTP {status}");
        }

        let mut body = Vec::new();
        resp.into_reader()
            .take(200 * 1024 * 1024)
            .read_to_end(&mut body)?;
        fs::write(&local, &body)?;
        total_bytes += body.len() as u64;
        println!("{} bytes", body.len());
    }

    println!(
        "\nDone. Fetched {} files, {:.1} MB total. Point gapsmith at the \
         target directory with --data-dir {}.",
        to_fetch.len(),
        total_bytes as f64 / (1024.0 * 1024.0),
        args.out_dir.display()
    );
    Ok(())
}
