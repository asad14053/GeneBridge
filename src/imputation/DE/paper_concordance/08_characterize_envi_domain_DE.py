#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

GLOBAL_FILE = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "envi_full_N23_domain_DE"
    / "ENVI_imputed_domain_adjusted_DE_SCZ_vs_NTC.csv"
)

DOMAIN_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "domain_stratified_N23_DE"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance"
    / "envi_domain_characterization_N23"
)

OUT.mkdir(parents=True, exist_ok=True)


DOMAINS = {
    "spd01": "WMtz",
    "spd02": "L3/4",
    "spd03": "L6",
    "spd04": "WM",
    "spd05": "L5",
    "spd06": "L2/3",
    "spd07": "L1/M",
}

FDR10 = 0.10
FDR05 = 0.05


def section(title: str) -> None:
    print("\n" + "=" * 110)
    print(title)
    print("=" * 110)


def require_columns(df: pd.DataFrame, label: str) -> None:
    required = {
        "gene",
        "logFC_SCZ",
        "t_stat_SCZ",
        "p_value_SCZ",
        "fdr_SCZ",
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"{label}: missing columns: {sorted(missing)}"
        )


def pearson(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")

    ok = np.isfinite(x) & np.isfinite(y)

    if ok.sum() < 3:
        return np.nan

    return float(
        np.corrcoef(
            x[ok],
            y[ok],
        )[0, 1]
    )


def spearman(x: pd.Series, y: pd.Series) -> float:
    xr = pd.to_numeric(x, errors="coerce").rank()
    yr = pd.to_numeric(y, errors="coerce").rank()

    return pearson(xr, yr)


# ======================================================================================
# LOAD GLOBAL DE
# ======================================================================================

section("LOAD GLOBAL DOMAIN-ADJUSTED DE")

if not GLOBAL_FILE.exists():
    raise FileNotFoundError(GLOBAL_FILE)

global_df = pd.read_csv(GLOBAL_FILE)

require_columns(
    global_df,
    "Global DE",
)

global_df["gene"] = global_df["gene"].astype(str)

if global_df["gene"].duplicated().any():
    raise RuntimeError("Duplicate genes in global DE table.")

global_df = global_df.set_index("gene", drop=False)

print("Global genes tested:", len(global_df))
print("Global FDR < 0.10:", int((global_df["fdr_SCZ"] < FDR10).sum()))
print("Global FDR < 0.05:", int((global_df["fdr_SCZ"] < FDR05).sum()))

if len(global_df) != 31995:
    raise RuntimeError(
        f"Expected 31,995 imputed genes; found {len(global_df)}"
    )


# ======================================================================================
# LOAD DOMAIN DE
# ======================================================================================

section("LOAD DOMAIN-STRATIFIED DE")

domain_tables: dict[str, pd.DataFrame] = {}

for domain, layer in DOMAINS.items():

    path = (
        DOMAIN_ROOT
        / domain
        / "ENVI_imputed_N23_domain_stratified_DE.csv"
    )

    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    require_columns(
        df,
        f"{domain} / {layer}",
    )

    df["gene"] = df["gene"].astype(str)

    if df["gene"].duplicated().any():
        raise RuntimeError(
            f"Duplicate genes in {domain}"
        )

    df = df.set_index("gene", drop=False)

    if set(df.index) != set(global_df.index):
        raise RuntimeError(
            f"Gene universe mismatch: {domain}"
        )

    df = df.loc[global_df.index]

    domain_tables[domain] = df

    print(
        f"{domain:5s} {layer:5s} | "
        f"genes={len(df):5d} | "
        f"FDR10={(df['fdr_SCZ'] < FDR10).sum():4d} | "
        f"FDR05={(df['fdr_SCZ'] < FDR05).sum():4d}"
    )


# ======================================================================================
# BUILD MASTER GENE-BY-DOMAIN TABLE
# ======================================================================================

section("BUILD MASTER TABLE")

master = pd.DataFrame(
    index=global_df.index
)

master["gene"] = global_df.index

master["global_logFC"] = global_df["logFC_SCZ"]
master["global_t"] = global_df["t_stat_SCZ"]
master["global_p"] = global_df["p_value_SCZ"]
master["global_fdr"] = global_df["fdr_SCZ"]

master["global_FDR10"] = (
    global_df["fdr_SCZ"] < FDR10
)

master["global_FDR05"] = (
    global_df["fdr_SCZ"] < FDR05
)


for domain, layer in DOMAINS.items():

    df = domain_tables[domain]

    master[f"{domain}_logFC"] = df["logFC_SCZ"]
    master[f"{domain}_t"] = df["t_stat_SCZ"]
    master[f"{domain}_p"] = df["p_value_SCZ"]
    master[f"{domain}_fdr"] = df["fdr_SCZ"]

    master[f"{domain}_FDR10"] = (
        df["fdr_SCZ"] < FDR10
    )

    master[f"{domain}_FDR05"] = (
        df["fdr_SCZ"] < FDR05
    )


domain_fdr10_cols = [
    f"{d}_FDR10"
    for d in DOMAINS
]

domain_fdr05_cols = [
    f"{d}_FDR05"
    for d in DOMAINS
]


master["n_domains_FDR10"] = (
    master[domain_fdr10_cols]
    .sum(axis=1)
    .astype(int)
)

master["n_domains_FDR05"] = (
    master[domain_fdr05_cols]
    .sum(axis=1)
    .astype(int)
)


master.to_csv(
    OUT / "08_gene_by_domain_master.csv",
    index=False,
)

print("Master table:", master.shape)


# ======================================================================================
# DOMAIN OVERLAP SUMMARY
# ======================================================================================

section("GLOBAL VS DOMAIN OVERLAP")

global10 = set(
    global_df.index[
        global_df["fdr_SCZ"] < FDR10
    ]
)

global05 = set(
    global_df.index[
        global_df["fdr_SCZ"] < FDR05
    ]
)


summary_rows = []

for domain, layer in DOMAINS.items():

    df = domain_tables[domain]

    sig10 = set(
        df.index[
            df["fdr_SCZ"] < FDR10
        ]
    )

    sig05 = set(
        df.index[
            df["fdr_SCZ"] < FDR05
        ]
    )

    shared_global10 = sig10 & global10
    domain_not_global10 = sig10 - global10

    other_domains = set()

    for other in DOMAINS:

        if other == domain:
            continue

        odf = domain_tables[other]

        other_domains |= set(
            odf.index[
                odf["fdr_SCZ"] < FDR10
            ]
        )

    unique_domain10 = sig10 - other_domains

    if sig10:
        direction_match_global = sum(
            np.sign(df.loc[g, "logFC_SCZ"])
            ==
            np.sign(global_df.loc[g, "logFC_SCZ"])
            for g in shared_global10
        )
    else:
        direction_match_global = 0

    n_shared_global = len(shared_global10)

    direction_fraction = (
        direction_match_global / n_shared_global
        if n_shared_global > 0
        else np.nan
    )

    summary_rows.append(
        {
            "domain": domain,
            "layer": layer,
            "FDR10": len(sig10),
            "FDR05": len(sig05),
            "shared_with_global_FDR10": len(shared_global10),
            "domain_FDR10_not_global": len(domain_not_global10),
            "unique_to_this_domain_FDR10": len(unique_domain10),
            "shared_global_direction_match_n":
                direction_match_global,
            "shared_global_direction_match_fraction":
                direction_fraction,
        }
    )


summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    OUT / "08_domain_overlap_summary.csv",
    index=False,
)

print(
    summary_df.to_string(
        index=False
    )
)


# ======================================================================================
# HOW MANY DOMAINS PER GENE?
# ======================================================================================

section("NUMBER OF SIGNIFICANT DOMAINS PER GENE")

multiplicity = (
    master["n_domains_FDR10"]
    .value_counts()
    .sort_index()
    .rename_axis("n_domains_FDR10")
    .reset_index(name="n_genes")
)

multiplicity.to_csv(
    OUT / "08_number_of_significant_domains_per_gene.csv",
    index=False,
)

print(
    multiplicity.to_string(
        index=False
    )
)


union_domain10 = master[
    master["n_domains_FDR10"] > 0
].copy()

union_domain10.to_csv(
    OUT / "08_union_domain_FDR10_genes.csv",
    index=False,
)


multi_domain = master[
    master["n_domains_FDR10"] >= 2
].copy()

multi_domain = multi_domain.sort_values(
    [
        "n_domains_FDR10",
        "global_fdr",
    ],
    ascending=[
        False,
        True,
    ],
)

multi_domain.to_csv(
    OUT / "08_multi_domain_FDR10_genes.csv",
    index=False,
)

print(
    "\nGenes significant in >=1 domain:",
    len(union_domain10)
)

print(
    "Genes significant in >=2 domains:",
    len(multi_domain)
)


# ======================================================================================
# PAIRWISE DOMAIN OVERLAP / JACCARD
# ======================================================================================

section("PAIRWISE DOMAIN JACCARD")

domains = list(DOMAINS)

jaccard = pd.DataFrame(
    np.nan,
    index=domains,
    columns=domains,
)

overlap_n = pd.DataFrame(
    0,
    index=domains,
    columns=domains,
    dtype=int,
)


sig_sets = {}

for d in domains:

    df = domain_tables[d]

    sig_sets[d] = set(
        df.index[
            df["fdr_SCZ"] < FDR10
        ]
    )


for d1 in domains:
    for d2 in domains:

        a = sig_sets[d1]
        b = sig_sets[d2]

        inter = len(a & b)
        union = len(a | b)

        overlap_n.loc[d1, d2] = inter

        jaccard.loc[d1, d2] = (
            inter / union
            if union > 0
            else np.nan
        )


jaccard.to_csv(
    OUT / "08_pairwise_domain_Jaccard_FDR10.csv"
)

overlap_n.to_csv(
    OUT / "08_pairwise_domain_overlap_counts_FDR10.csv"
)

print("\nOverlap counts:")
print(overlap_n)

print("\nJaccard:")
print(
    jaccard.round(4)
)


# ======================================================================================
# GLOBAL VS DOMAIN EFFECT CONCORDANCE
# ======================================================================================

section("GLOBAL VS DOMAIN EFFECT CONCORDANCE")

corr_rows = []

for domain, layer in DOMAINS.items():

    df = domain_tables[domain]

    corr_rows.append(
        {
            "domain": domain,
            "layer": layer,
            "logFC_Pearson_all_genes":
                pearson(
                    global_df["logFC_SCZ"],
                    df["logFC_SCZ"],
                ),
            "logFC_Spearman_all_genes":
                spearman(
                    global_df["logFC_SCZ"],
                    df["logFC_SCZ"],
                ),
            "t_Pearson_all_genes":
                pearson(
                    global_df["t_stat_SCZ"],
                    df["t_stat_SCZ"],
                ),
            "t_Spearman_all_genes":
                spearman(
                    global_df["t_stat_SCZ"],
                    df["t_stat_SCZ"],
                ),
        }
    )


corr_df = pd.DataFrame(corr_rows)

corr_df.to_csv(
    OUT / "08_global_vs_domain_effect_concordance.csv",
    index=False,
)

print(
    corr_df.to_string(
        index=False
    )
)


# ======================================================================================
# TOP GENES PER DOMAIN
# ======================================================================================

section("TOP GENES PER DOMAIN")

top_all = []

for domain, layer in DOMAINS.items():

    df = domain_tables[domain].copy()

    df["abs_logFC"] = (
        df["logFC_SCZ"].abs()
    )

    df["domain"] = domain
    df["layer"] = layer

    top = (
        df.sort_values(
            [
                "fdr_SCZ",
                "p_value_SCZ",
                "abs_logFC",
            ],
            ascending=[
                True,
                True,
                False,
            ],
        )
        .head(25)
        .copy()
    )

    top["global_logFC"] = (
        global_df.loc[
            top.index,
            "logFC_SCZ"
        ].values
    )

    top["global_fdr"] = (
        global_df.loc[
            top.index,
            "fdr_SCZ"
        ].values
    )

    top["global_FDR10"] = (
        top["global_fdr"] < FDR10
    )

    top["direction_matches_global"] = (
        np.sign(top["logFC_SCZ"])
        ==
        np.sign(top["global_logFC"])
    )

    keep_cols = [
        "gene",
        "domain",
        "layer",
        "logFC_SCZ",
        "t_stat_SCZ",
        "p_value_SCZ",
        "fdr_SCZ",
        "global_logFC",
        "global_fdr",
        "global_FDR10",
        "direction_matches_global",
    ]

    top_all.append(
        top[keep_cols]
    )


top_df = pd.concat(
    top_all,
    ignore_index=True,
)

top_df.to_csv(
    OUT / "08_top25_genes_per_domain.csv",
    index=False,
)


# ======================================================================================
# SIGNIFICANT GENE LONG TABLE
# ======================================================================================

section("CREATE LONG SIGNIFICANT-GENE TABLE")

long_rows = []

for domain, layer in DOMAINS.items():

    df = domain_tables[domain]

    sig = df[
        df["fdr_SCZ"] < FDR10
    ].copy()

    for gene, row in sig.iterrows():

        long_rows.append(
            {
                "gene": gene,
                "domain": domain,
                "layer": layer,
                "logFC_SCZ": row["logFC_SCZ"],
                "t_stat_SCZ": row["t_stat_SCZ"],
                "p_value_SCZ": row["p_value_SCZ"],
                "fdr_SCZ": row["fdr_SCZ"],
                "global_logFC":
                    global_df.loc[
                        gene,
                        "logFC_SCZ",
                    ],
                "global_fdr":
                    global_df.loc[
                        gene,
                        "fdr_SCZ",
                    ],
                "global_FDR10":
                    bool(
                        global_df.loc[
                            gene,
                            "fdr_SCZ",
                        ]
                        < FDR10
                    ),
                "n_domains_FDR10":
                    int(
                        master.loc[
                            gene,
                            "n_domains_FDR10",
                        ]
                    ),
            }
        )


long_df = pd.DataFrame(long_rows)

if not long_df.empty:

    long_df = long_df.sort_values(
        [
            "domain",
            "fdr_SCZ",
        ]
    )


long_df.to_csv(
    OUT / "08_domain_FDR10_genes_long.csv",
    index=False,
)


# ======================================================================================
# FIGURE 1: DOMAIN DEG COUNTS
# ======================================================================================

section("PLOTS")

plot_df = summary_df.copy()

x = np.arange(
    len(plot_df)
)

width = 0.38

fig, ax = plt.subplots(
    figsize=(10, 6)
)

ax.bar(
    x - width / 2,
    plot_df["FDR10"],
    width,
    label="FDR < 0.10",
)

ax.bar(
    x + width / 2,
    plot_df["FDR05"],
    width,
    label="FDR < 0.05",
)

ax.set_xticks(x)
ax.set_xticklabels(
    plot_df["layer"]
)

ax.set_ylabel(
    "Number of ENVI-imputed DE genes"
)

ax.set_xlabel(
    "Spatial domain"
)

ax.set_title(
    "SCZ-associated ENVI-imputed genes by spatial domain"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUT / "08A_domain_DEG_counts.png",
    dpi=300,
)

fig.savefig(
    OUT / "08A_domain_DEG_counts.pdf",
)

plt.close(fig)


# ======================================================================================
# FIGURE 2: SHARED WITH GLOBAL VS DOMAIN-ONLY
# ======================================================================================

fig, ax = plt.subplots(
    figsize=(10, 6)
)

ax.bar(
    x,
    plot_df["shared_with_global_FDR10"],
    label="Also global FDR < 0.10",
)

ax.bar(
    x,
    plot_df["domain_FDR10_not_global"],
    bottom=plot_df["shared_with_global_FDR10"],
    label="Domain FDR < 0.10 only",
)

ax.set_xticks(x)

ax.set_xticklabels(
    plot_df["layer"]
)

ax.set_ylabel(
    "Number of FDR < 0.10 genes"
)

ax.set_xlabel(
    "Spatial domain"
)

ax.set_title(
    "Global versus domain-specific ENVI-imputed SCZ signals"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUT / "08B_global_vs_domain_overlap.png",
    dpi=300,
)

fig.savefig(
    OUT / "08B_global_vs_domain_overlap.pdf",
)

plt.close(fig)


# ======================================================================================
# FIGURE 3: PAIRWISE JACCARD HEATMAP
# ======================================================================================

fig, ax = plt.subplots(
    figsize=(8, 7)
)

mat = jaccard.to_numpy(
    dtype=float
)

im = ax.imshow(
    mat,
    aspect="auto",
)

labels = [
    DOMAINS[d]
    for d in domains
]

ax.set_xticks(
    np.arange(len(domains))
)

ax.set_yticks(
    np.arange(len(domains))
)

ax.set_xticklabels(
    labels,
    rotation=45,
    ha="right",
)

ax.set_yticklabels(
    labels
)

ax.set_title(
    "Pairwise overlap of domain-specific DE genes\nJaccard index, FDR < 0.10"
)

for i in range(len(domains)):
    for j in range(len(domains)):

        value = mat[i, j]

        if np.isfinite(value):

            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
            )

fig.colorbar(
    im,
    ax=ax,
    label="Jaccard index",
)

fig.tight_layout()

fig.savefig(
    OUT / "08C_domain_Jaccard_heatmap.png",
    dpi=300,
)

fig.savefig(
    OUT / "08C_domain_Jaccard_heatmap.pdf",
)

plt.close(fig)


# ======================================================================================
# FINAL
# ======================================================================================

section("FINAL SUMMARY")

print("Imputed genes tested:", len(global_df))
print("Global FDR < 0.10:", len(global10))
print("Global FDR < 0.05:", len(global05))

print("\nDomain summary:")
print(
    summary_df[
        [
            "domain",
            "layer",
            "FDR10",
            "FDR05",
            "shared_with_global_FDR10",
            "domain_FDR10_not_global",
            "unique_to_this_domain_FDR10",
            "shared_global_direction_match_fraction",
        ]
    ].to_string(
        index=False
    )
)

print(
    "\nGenes significant in >=1 domain:",
    int((master["n_domains_FDR10"] >= 1).sum())
)

print(
    "Genes significant in >=2 domains:",
    int((master["n_domains_FDR10"] >= 2).sum())
)

print(
    "Genes significant in >=3 domains:",
    int((master["n_domains_FDR10"] >= 3).sum())
)

print(
    "Genes significant in >=4 domains:",
    int((master["n_domains_FDR10"] >= 4).sum())
)

print("\nFINAL STATUS: DOMAIN CHARACTERIZATION COMPLETE")
