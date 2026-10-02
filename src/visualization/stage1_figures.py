"""Stage 1 visualization suite — programmatic generation of all required figures.

All figures use semantic colors from configs/visualization.yaml.
All figures are exported as PNG (300 dpi) and SVG.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import yaml

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")
COLORS_CONFIG = yaml.safe_load((PROJECT_ROOT / "configs" / "visualization.yaml").read_text())

# Extract semantic colors for use
COLORS = {
    "vhh": COLORS_CONFIG["colors"]["semantic"]["vhh_sdab"]["hex"],
    "vhh_light": COLORS_CONFIG["colors"]["semantic"]["vhh_sdab_light"]["hex"],
    "conv_ab": COLORS_CONFIG["colors"]["semantic"]["conventional_antibody"]["hex"],
    "conv_ab_light": COLORS_CONFIG["colors"]["semantic"]["conventional_antibody_light"]["hex"],
    "conv_vh": COLORS_CONFIG["colors"]["semantic"]["conventional_vh"]["hex"],
    "antigen": COLORS_CONFIG["colors"]["annotation"]["protein_antigen"]["hex"],
    "grey": COLORS_CONFIG["colors"]["annotation"]["missing_excluded"]["hex"],
}

# Ensure figure directories exist
for subdir in ["qc", "main", "supplementary"]:
    (PROJECT_ROOT / "results" / "figures" / "stage1" / subdir).mkdir(parents=True, exist_ok=True)

def save_figure(fig, fig_id, title, subdir="main"):
    """Save figure as both PNG (300 dpi) and SVG."""
    base_path = PROJECT_ROOT / "results" / "figures" / "stage1" / subdir / fig_id
    png_path = base_path.with_suffix(".png")
    svg_path = base_path.with_suffix(".svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    print(f"  saved {png_path.name}")
    print(f"  saved {svg_path.name}")
    return str(png_path), str(svg_path)

def load_data():
    """Load the Stage-1 resolved metadata."""
    resolved = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet")
    return resolved

def fig_s1_01_curation_flow(resolved):
    """FIG_S1_01: Dataset and curation flow."""
    print("Generating FIG_S1_01: Dataset and curation flow...")

    # Collect key counts
    sdh = resolved[resolved["source_cohort"] == "SD-H"]
    fab = resolved[resolved["source_cohort"] == "FAB"]
    fv = resolved[resolved["source_cohort"] == "FV"]

    sdh_natural = (sdh["sdh_classification"] == "natural_camelid_candidate").sum()
    sdh_eng = (sdh["sdh_classification"] == "engineered_or_humanized_camelid_candidate").sum()
    sdh_syn = (sdh["sdh_classification"] == "synthetic_or_de_novo_candidate").sum()
    sdh_human = (sdh["sdh_classification"] == "human_sdVH_candidate").sum()
    sdh_other = (sdh["sdh_classification"] == "other_noncamelid_sdH").sum()
    sdh_unc = (sdh["sdh_classification"] == "uncertain_sdH").sum()

    fv_conv = (fv["construct_bucket"] == "CONVENTIONAL").sum()
    fv_scfv = (fv["construct_bucket"] == "SCFV").sum()
    fv_dia = (fv["construct_bucket"] == "DIABODY").sum()

    strict_nat = int(sdh["primary_vhh_candidate_strict"].sum())
    mixed_prov = int((sdh["sequence_provenance_ambiguity"] == True).sum())

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: source cohorts and SD-H classification
    categories = ["Natural\nCamelid", "Engineered/\nHumanized", "Synthetic", "Human sdVH", "Other", "Uncertain"]
    counts = [sdh_natural, sdh_eng, sdh_syn, sdh_human, sdh_other, sdh_unc]
    colors_l = [COLORS["vhh"], COLORS["conv_ab_light"], COLORS["grey"], COLORS["grey"], COLORS["grey"], COLORS["grey"]]

    ax1.barh(categories, counts, color=colors_l)
    ax1.set_xlabel("Number of SD-H rows", fontsize=11)
    ax1.set_title("SD-H Classification Breakdown\n(n=1597 total)", fontsize=12, fontweight="bold")
    ax1.grid(axis="x", alpha=0.3)
    for i, (cat, cnt) in enumerate(zip(categories, counts)):
        ax1.text(cnt + 20, i, str(cnt), va="center", fontsize=10)

    # Right: FV construct and primary cohort tiers
    tiers = ["Strict Natural\nVHH (primary)", "Mixed Provenance\n(structural QC)", "Engineered/\nHumanized", "Manual\nReview"]
    tier_counts = [strict_nat, mixed_prov, sdh_eng, sdh_unc - 2]  # -2 for the 2 uncertain in mixed-prov
    tier_colors = [COLORS["vhh"], COLORS["vhh_light"], COLORS["conv_ab_light"], COLORS["grey"]]

    ax2.barh(tiers, tier_counts, color=tier_colors)
    ax2.set_xlabel("Number of SD-H rows", fontsize=11)
    ax2.set_title("Final VHH Eligibility Tiers (Stage 1.5)", fontsize=12, fontweight="bold")
    ax2.grid(axis="x", alpha=0.3)
    for i, (tier, cnt) in enumerate(zip(tiers, tier_counts)):
        ax2.text(cnt + 20, i, str(cnt), va="center", fontsize=10)

    fig.suptitle("Stage 1: Dataset Curation and Classification", fontsize=14, fontweight="bold", y=1.00)
    plt.tight_layout()
    return save_figure(fig, "FIG_S1_01_curation_flow", "Dataset and curation flow", subdir="main")

def fig_s1_02_sdh_classification(resolved):
    """FIG_S1_02: SD-H classification distribution."""
    print("Generating FIG_S1_02: SD-H classification distribution...")
    sdh = resolved[resolved["source_cohort"] == "SD-H"]

    class_counts = sdh["sdh_classification"].value_counts()
    class_order = ["natural_camelid_candidate", "human_sdVH_candidate", "synthetic_or_de_novo_candidate",
                   "engineered_or_humanized_camelid_candidate", "uncertain_sdH", "other_noncamelid_sdH"]
    class_order = [c for c in class_order if c in class_counts.index]
    class_counts = class_counts[class_order]

    class_labels = ["Natural Camelid", "Human sdVH", "Synthetic", "Engineered/Humanized", "Uncertain", "Other"]
    class_labels = [l for i, l in enumerate(class_labels) if class_order[i] in class_counts.index]

    colors_use = [COLORS["vhh"] if c == "Natural Camelid" else COLORS["grey"] for c in class_labels]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(class_labels, class_counts.values, color=colors_use)
    ax.set_xlabel("Number of rows", fontsize=11)
    ax.set_title("Stage 1 SD-H Classification Distribution", fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    for i, (bar, cnt) in enumerate(zip(bars, class_counts.values)):
        pct = 100 * cnt / sdh.shape[0]
        ax.text(cnt + 10, i, f"{cnt}\n({pct:.1f}%)", va="center", fontsize=10)

    plt.tight_layout()
    return save_figure(fig, "FIG_S1_02_sdh_classification", "SD-H classification distribution", subdir="main")

def fig_s1_03_species_distribution(resolved):
    """FIG_S1_03: Primary VHH species distribution."""
    print("Generating FIG_S1_03: Primary VHH species distribution...")
    sdh = resolved[resolved["source_cohort"] == "SD-H"]
    primary = sdh[sdh["primary_vhh_candidate_strict"] == True]

    species_map = {
        "9838": "C. dromedarius",
        "9837": "C. bactrianus",
        "419612": "C. ferus",
        "9844": "L. glama",
        "30538": "V. pacos",
    }

    species_counts = []
    species_ids = []
    for taxid, label in sorted(species_map.items()):
        sub = primary[primary["heavy_taxid"] == taxid]
        if len(sub) > 0:
            species_counts.append(len(sub))
            species_ids.append(sub["SABDAB_ID"].nunique())
            species_counts.append(sub["SABDAB_ID"].nunique())
            species_ids.append(label)

    # More direct approach
    species_list = ["C. dromedarius", "C. bactrianus", "L. glama", "V. pacos", "Other/Unresolved"]
    taxids = ["9838", "9837", "9844", "30538", None]

    row_counts = []
    id_counts = []
    for taxid in taxids:
        if taxid:
            sub = primary[primary["heavy_taxid"] == taxid]
        else:
            sub = primary[~primary["heavy_taxid"].isin(["9838", "9837", "9844", "30538"])]
        row_counts.append(len(sub))
        id_counts.append(sub["SABDAB_ID"].nunique())

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(species_list))
    width = 0.35
    ax.bar(x - width/2, row_counts, width, label="Instance rows", color=COLORS["vhh"])
    ax.bar(x + width/2, id_counts, width, label="Unique SABDAB_IDs", color=COLORS["vhh_light"])

    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("Primary VHH Candidates: Species Distribution (n=1328 instances)", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(species_list)
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    for i, (rc, ic) in enumerate(zip(row_counts, id_counts)):
        ax.text(i - width/2, rc + 10, str(rc), ha="center", fontsize=9)
        ax.text(i + width/2, ic + 10, str(ic), ha="center", fontsize=9)

    plt.tight_layout()
    return save_figure(fig, "FIG_S1_03_species_distribution", "Primary VHH species distribution", subdir="main")

def fig_s1_04_fv_construct(resolved):
    """FIG_S1_04: FV construct architecture."""
    print("Generating FIG_S1_04: FV construct architecture...")
    fv = resolved[resolved["source_cohort"] == "FV"]

    construct_counts = fv["construct_bucket"].value_counts()
    construct_order = ["CONVENTIONAL", "SCFV", "DIABODY", "CODV", "DVD", "SCDB"]
    construct_order = [c for c in construct_order if c in construct_counts.index]
    construct_counts = construct_counts[construct_order]

    colors_use = [COLORS["vhh"] if c == "CONVENTIONAL" else COLORS["grey"] for c in construct_order]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(construct_order, construct_counts.values, color=colors_use)
    ax.set_ylabel("Number of FV rows", fontsize=11)
    ax.set_title("FV Construct Architecture Distribution (n=536 total)\nOnly CONVENTIONAL (n=294) is primary paired-antibody cohort", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    for bar, cnt in zip(bars, construct_counts.values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 5, str(int(cnt)), ha="center", fontsize=10, fontweight="bold")

    plt.tight_layout()
    return save_figure(fig, "FIG_S1_04_fv_construct", "FV construct architecture", subdir="main")

def fig_s1_05_provenance_tiers(resolved):
    """FIG_S1_05: Final VHH provenance tiers."""
    print("Generating FIG_S1_05: Final VHH provenance tiers...")
    sdh = resolved[resolved["source_cohort"] == "SD-H"]

    status_counts = sdh["provenance_resolution_status"].value_counts()
    status_order = ["RESOLVED_NATURAL_CAMELID", "MIXED_PROVENANCE_SAME_SEQUENCE",
                    "RESOLVED_ENGINEERED_OR_HUMANIZED", "MANUAL_REVIEW"]
    status_labels = ["Strict Natural\nVHH", "Mixed-Provenance\nSequence", "Engineered/\nHumanized", "Manual\nReview"]

    counts = [status_counts.get(s, 0) for s in status_order]
    colors_use = [COLORS["vhh"], COLORS["vhh_light"], COLORS["conv_ab_light"], COLORS["grey"]]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(status_labels, counts, color=colors_use)
    ax.set_ylabel("Number of SD-H rows", fontsize=11)
    ax.set_title("Final VHH Structural-QC Tiers (Stage 1.5 Derived Eligibility)", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    for bar, cnt in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 10, str(int(cnt)), ha="center", fontsize=10, fontweight="bold")

    # Add annotation about 91 rows
    fig.text(0.5, 0.02, "Note: 91 of the 108 mixed-provenance rows were previously in the natural-camelid primary cohort; they were downgraded, not excluded.",
             ha="center", fontsize=9, style="italic", wrap=True)
    plt.subplots_adjust(bottom=0.12)

    plt.tight_layout()
    return save_figure(fig, "FIG_S1_05_provenance_tiers", "Final VHH provenance tiers", subdir="main")

def fig_s1_06_redundancy(resolved):
    """FIG_S1_06: Redundancy and pseudoreplication."""
    print("Generating FIG_S1_06: Redundancy and pseudoreplication...")
    sdh = resolved[resolved["source_cohort"] == "SD-H"]
    primary_vhh = sdh[sdh["primary_vhh_candidate_strict"] == True]

    # Count instances per SABDAB_ID
    instances_per_id_vhh = primary_vhh.groupby("SABDAB_ID").size()

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.hist(instances_per_id_vhh, bins=range(1, instances_per_id_vhh.max() + 2),
            color=COLORS["vhh"], edgecolor="black", alpha=0.7)
    ax.set_xlabel("Instances per SABDAB_ID", fontsize=11)
    ax.set_ylabel("Number of SABDAB_IDs", fontsize=11)
    ax.set_title(f"Pseudoreplication: VHH Strict Primary (n={len(instances_per_id_vhh)} unique IDs, {len(primary_vhh)} total rows)", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    # Add statistics
    stats_text = f"Median: {instances_per_id_vhh.median():.0f}\nMean: {instances_per_id_vhh.mean():.2f}\nMax: {instances_per_id_vhh.max():.0f}"
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    return save_figure(fig, "FIG_S1_06_redundancy", "Redundancy and pseudoreplication", subdir="qc")

def fig_s1_07_metadata_quality(resolved):
    """FIG_S1_07: Metadata quality (resolution, R-factor)."""
    print("Generating FIG_S1_07: Metadata quality...")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # Resolution by cohort
    for cohort, color in [("SD-H", COLORS["vhh"]), ("FAB", COLORS["conv_ab"]), ("FV", COLORS["conv_vh"])]:
        sub = resolved[resolved["source_cohort"] == cohort]
        res_data = pd.to_numeric(sub["resolution"], errors="coerce").dropna()
        ax1.hist(res_data, bins=30, alpha=0.5, label=f"{cohort} (n={len(res_data)})", color=color)
    ax1.set_xlabel("Resolution (Å)", fontsize=11)
    ax1.set_ylabel("Count", fontsize=11)
    ax1.set_title("Resolution Distribution by Cohort", fontsize=12, fontweight="bold")
    ax1.legend()
    ax1.grid(axis="y", alpha=0.3)
    ax1.set_xlim(1, 3.5)

    # R-factor presence
    r_factor_data = resolved["r_factor"].apply(lambda x: pd.notna(x))
    by_cohort = resolved.groupby("source_cohort")["r_factor"].apply(lambda x: pd.notna(x).sum())
    total_by_cohort = resolved["source_cohort"].value_counts()
    completeness = {c: by_cohort.get(c, 0) / total_by_cohort[c] for c in total_by_cohort.index}

    cohorts = list(completeness.keys())
    values = [completeness[c] * 100 for c in cohorts]
    colors_use = [COLORS["vhh"], COLORS["conv_ab"], COLORS["conv_vh"]]
    ax2.bar(cohorts, values, color=colors_use)
    ax2.set_ylabel("% with r_factor value", fontsize=11)
    ax2.set_title("r_factor Data Completeness", fontsize=12, fontweight="bold")
    ax2.set_ylim(0, 105)
    ax2.grid(axis="y", alpha=0.3)
    for i, (c, v) in enumerate(zip(cohorts, values)):
        ax2.text(i, v + 2, f"{v:.0f}%", ha="center", fontsize=10)

    # R-free presence
    r_free_data = resolved["r_free"].apply(lambda x: pd.notna(x))
    by_cohort_free = resolved.groupby("source_cohort")["r_free"].apply(lambda x: pd.notna(x).sum())
    completeness_free = {c: by_cohort_free.get(c, 0) / total_by_cohort[c] for c in total_by_cohort.index}
    values_free = [completeness_free[c] * 100 for c in cohorts]
    ax3.bar(cohorts, values_free, color=colors_use)
    ax3.set_ylabel("% with r_free value", fontsize=11)
    ax3.set_title("r_free Data Completeness", fontsize=12, fontweight="bold")
    ax3.set_ylim(0, 105)
    ax3.grid(axis="y", alpha=0.3)
    for i, (c, v) in enumerate(zip(cohorts, values_free)):
        ax3.text(i, v + 2, f"{v:.0f}%", ha="center", fontsize=10)

    # Sample sizes
    sizes = [len(resolved[resolved["source_cohort"] == c]) for c in ["SD-H", "FAB", "FV"]]
    ax4.bar(["SD-H", "FAB", "FV"], sizes, color=[COLORS["vhh"], COLORS["conv_ab"], COLORS["conv_vh"]])
    ax4.set_ylabel("Number of rows", fontsize=11)
    ax4.set_title("Sample Sizes by Cohort", fontsize=12, fontweight="bold")
    ax4.grid(axis="y", alpha=0.3)
    for i, (c, s) in enumerate(zip(["SD-H", "FAB", "FV"], sizes)):
        ax4.text(i, s + 30, str(s), ha="center", fontsize=10, fontweight="bold")

    fig.suptitle("Stage 1 Metadata Quality Assessment", fontsize=13, fontweight="bold")
    plt.tight_layout()
    return save_figure(fig, "FIG_S1_07_metadata_quality", "Metadata quality", subdir="qc")

def fig_s1_08_mixed_provenance_audit(resolved):
    """FIG_S1_08: Mixed-provenance audit."""
    print("Generating FIG_S1_08: Mixed-provenance audit...")
    sdh = resolved[resolved["source_cohort"] == "SD-H"]

    mixed_ids = ["sabdab2_H00QEL0000", "sabdab2_H021EL0000", "sabdab2_H02DAL0000",
                 "sabdab2_H02Y5L0000", "sabdab2_H06BKL0000"]

    fig, axes = plt.subplots(1, 5, figsize=(16, 4))
    for idx, mixed_id in enumerate(mixed_ids):
        ax = axes[idx]
        sub = sdh[sdh["SABDAB_ID"] == mixed_id]

        class_dist = sub["sdh_classification"].value_counts()
        class_labels = class_dist.index.tolist()
        class_colors = [COLORS["vhh"] if "natural" in c else COLORS["grey"] for c in class_labels]

        ax.barh(range(len(class_labels)), class_dist.values, color=class_colors)
        ax.set_yticks(range(len(class_labels)))
        ax.set_yticklabels([l.replace("_", "\n") for l in class_labels], fontsize=8)
        ax.set_xlabel("Rows", fontsize=9)
        ax.set_title(f"{mixed_id.replace('sabdab2_', '')}\n(n={len(sub)})", fontsize=10, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)

        for i, val in enumerate(class_dist.values):
            ax.text(val + 0.2, i, str(val), va="center", fontsize=8)

    fig.suptitle("Mixed-Provenance SAbDab2 IDs: Classification Composition\n(Note: Sequences are identical; disagreement is metadata only)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    return save_figure(fig, "FIG_S1_08_mixed_provenance_audit", "Mixed-provenance audit", subdir="main")

def generate_all_stage1_figures(resolved):
    """Generate all Stage-1 figures."""
    print("\n=== GENERATING STAGE 1 FIGURES ===\n")

    fig_list = [
        fig_s1_01_curation_flow(resolved),
        fig_s1_02_sdh_classification(resolved),
        fig_s1_03_species_distribution(resolved),
        fig_s1_04_fv_construct(resolved),
        fig_s1_05_provenance_tiers(resolved),
        fig_s1_06_redundancy(resolved),
        fig_s1_07_metadata_quality(resolved),
        fig_s1_08_mixed_provenance_audit(resolved),
    ]

    return fig_list

if __name__ == "__main__":
    resolved = load_data()
    fig_list = generate_all_stage1_figures(resolved)
    print(f"\n=== COMPLETED: {len(fig_list)} Stage-1 figures generated ===")
