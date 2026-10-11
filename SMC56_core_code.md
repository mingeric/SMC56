# SMC56 core analysis code and data sources

Core computation code for the SMC5/6 pan-cancer figure set, one section per figure.
All data files are CSV/TSV tables read from the delivery folder 00_source_data (145
files); file names in the Data lines refer to that folder. The full builder scripts,
including all plotting and export code, are kept next to the data in the submission
package.

## Shared style and gene-set definitions (smc56_style.py)

```python
# Style/palette base imported by every figure builder (rcParams, panel helpers omitted).
# BASE: local project root (absolute path on the build machine)
DATA = os.path.join(BASE, "SMC56_FINAL_DELIVERY", "SOURCE_DATA_AND_CODE", "00_source_data")
OUT = os.path.join(BASE, "MODULE_A_REBUILD_20261011", "figures")

# SMC5/6 complex: 10-gene panel and 4 composite modules
GENES10 = ["SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
           "EID3", "SLF1", "SLF2", "SIMC1"]
GROUPS4 = ["C1_Core_SMC56", "C2_Canonical", "C3_SLF_repair", "C4_SIMC"]

# subunit -> composite module membership (S4 contract; SLF2 in C3 and C4; EID3 unassigned)
GENE2GROUPS = {
    "SMC5": ["C1_Core_SMC56", "C2_Canonical"], "SMC6": ["C1_Core_SMC56", "C2_Canonical"],
    "NSMCE1": ["C2_Canonical"], "NSMCE2": ["C2_Canonical"],
    "NSMCE3": ["C2_Canonical"], "NSMCE4A": ["C2_Canonical"],
    "EID3": [],
    "SLF1": ["C3_SLF_repair"], "SLF2": ["C3_SLF_repair", "C4_SIMC"],
    "SIMC1": ["C4_SIMC"],
}

# diverging cmap: teal (low) -> ivory (mid) -> red (high)
CMAP_DIV = LinearSegmentedColormap.from_list(
    "smc56_div", ["#5E9C96", "#A7D8D2", "#EFE8DF", "#D9A89A", "#8F3B3B", "#4B1F20"], N=256)
# sequential cmap: cream -> peach -> red -> wine
CMAP_SEQ = LinearSegmentedColormap.from_list(
    "smc56_seq", ["#FFF7E3", "#F7D7A8", "#D9A89A", "#8F3B3B", "#4B1F20"], N=256)
# correlation-matrix cmap: slate-gray -> ivory -> red (symmetric about 0)
CMAP_CORR = LinearSegmentedColormap.from_list(
    "smc56_corr", ["#6E9EA3", "#A7D8D2", "#EFE8DF", "#D9A89A", "#8F3B3B"], N=256)

def load(fname):
    return pd.read_csv(os.path.join(DATA, fname))

def save(fig, name, dpi=300):
    """Export PNG (display) + PDF (vector archive). Returns the PNG path."""
    os.makedirs(OUT, exist_ok=True)
    png = os.path.join(OUT, f"{name}.png")
    pdf = os.path.join(OUT, f"{name}.pdf")
    fig.savefig(png, dpi=dpi)
    fig.savefig(pdf)
    print(f"[saved] {png}")
    return png
```

## Module A - Pan-cancer expression baseline

### FigF0 Quantitative formula framework (build_A_figF0_formulas.py)
Data: none (methodological overview figure; quantitative definitions per 0511SMC56_formula.docx, no CSV input).

```python
# (a) Step 1 - standardisation and stoichiometric imbalance
# z-score of gene g in sample i (x = raw expression; mu, sigma = pan-cancer cohort gene mean/sd):
r"$z_{i,g}=\dfrac{x_{i,g}-\mu_g}{\sigma_g}$"
# stoichiometric imbalance, variance layer over n=9 core subunits (zbar_i = within-sample mean;
# larger SI -> poorer subunit coordination):
r"$\mathrm{SI}_i=\dfrac{1}{n-1}\sum_{g=1}^{n}\left(z_{i,g}-\bar{z}_i\right)^{2}$"

# (b) Step 2 - 9 core subunits decomposed into 3 basal sub-complexes
rows = [
    ("Core scaffold complex", "SMC5 · SMC6 · NSMCE2(NSE2)", "ATPase activity and structural basis · panel covered"),
    ("Regulatory complex", "NSMCE1(NSE1) · NSMCE3(NSE3) · NSMCE4A(NSE4)", "SUMO/ubiquitin modification centre · covered (EID3 grouping tentative)"),
    ("Basal localisation complex", "NSMCE5(NSE5) · NSMCE6(NSE6)", "localisation module / SIF anchor · not in panel (evidence pending)"),
]

# (c) Step 3 - SIF2 dynamic recruitment efficiency (barrel effect: complex output is capped
#     by the less-expressed partner). Repair complex R = SIF1(SLF1)+SIF2(SLF2), DNA damage
#     repair; silencing complex S = MC1(SIMC1)+SIF2(SLF2), antiviral/epigenetic silencing.
r"recruitment efficiency $E_{\mathrm{R/S}}=\min\left(\bar{z}_{\mathrm{partner1}},\ \bar{z}_{\mathrm{SIF2}}\right)$"

# (d) Step 4 - competitive balance score (eps: small constant avoiding division by zero;
#     R/S = repair/silencing complex activity; CBS>0 repair-dominated -> DDRi sensitivity,
#     CBS<0 silencing-dominated -> ICB resistance hint)
r"$\mathrm{CBS}=\dfrac{R-S}{R+S+\varepsilon}$"

# (e) Panel link: C1-C4 group scores = mean of member z-scores ("expression" layer, S4
#     definition); 6/9 core subunits covered (SMC5/6, NSE1-4) plus the SIF1/SIF2/MC1
#     recruitment layer; NSMCE5/6 not included. All metrics are expression-level proxies,
#     not evidence of complex assembly.
```

### FigA1 Cohort flow and gene-group map (build_A_figA1.py)
Data: tcga_patient_gene_log2tpm.csv, F1_tumour_vs_normal_12cohorts.csv, S4_sc_group_score_contract.csv — 9,523 patients across 31 cancer types (all with 10-gene log2TPM); 12 paired cohorts with n_tumour=5,846 and n_normal=607.

```python
pat = load("tcga_patient_gene_log2tpm.csv")
n_total = len(pat)
cc = pat["cancer"].value_counts()

f1 = load("F1_tumour_vs_normal_12cohorts.csv")
coh = f1.groupby("cohort")[["n_tumour", "n_normal"]].first()
coh = coh.loc[[c for c in coh.index if c in cc.index]]           # keep the 12 paired cohorts
coh = coh.sort_values("n_tumour", ascending=False)
n_cohort_patient = int(sum(cc.get(c, 0) for c in coh.index))
n_rest = n_total - n_cohort_patient
n_rest_types = int(cc.index.size - coh.index.size)

# sankey stages: 0 all patients -> 1 by cancer type -> 2 transcriptome available -> 3 tumour/normal
s0 = [(n_total, PAL["wine"], f"TCGA pan-cancer patients  n = {n_total:,}")]
s1 = [(int(cc[c]), PAL["red"], c) for c in coh.index]
s1.append((n_rest, PAL["warmgray"], f"other {n_rest_types} cancer types"))
s2 = [(n_total, PAL["slate"],
       f"transcriptome available  n = {n_total:,} (10-gene log2TPM complete)")]
s3 = [(int(r.n_tumour + r.n_normal), None, c) for c, r in coh.iterrows()]

# node stacking: bar heights proportional to counts within a fixed vertical span
Y0, Y1 = 4.0, 95.0
H = Y1 - Y0

def stack(items, gap, total_h=None):
    """items=[(val,color,label)]; returns [(val,color,label,y0,y1)]"""
    tot = sum(v for v, _, _ in items)
    avail = H - gap * (len(items) - 1)
    sc = avail / tot
    out, y = [], Y0
    for v, c, l in items:
        h = v * sc
        out.append((v, c, l, y, y + h))
        y += h + gap
    return out

# panel b: gene-group membership dictionary read back from the S4 contract
s4 = load("S4_sc_group_score_contract.csv")
mem = {g: set(str(s4.loc[s4["group"] == g, "members"].iloc[0]).split(" + ")) for g in GROUPS4}
# scoring rule: group score = mean of member-gene log2TPM z-scores (S4 contract)
```

### FigA2 GTEx tissue baseline (build_A_figA2.py)
Data: GTEx_10gene_median_tpm.csv — median TPM for 10 genes across 54 normal tissues.

```python
# gene -> group mapping (by membership)
GENE2GROUP = {}
for gene in GENES10:
    if gene in ("SMC5", "SMC6"):
        GENE2GROUP[gene] = "C1_Core_SMC56"
    elif gene.startswith("NSMCE"):
        GENE2GROUP[gene] = "C2_Canonical"
    elif gene in ("SLF1", "SLF2"):
        GENE2GROUP[gene] = "C3_SLF_repair"
    else:
        GENE2GROUP[gene] = "C4_SIMC"

gt = load("GTEx_10gene_median_tpm.csv").set_index("Description").drop(columns=["Name"])
tissues = list(gt.columns)
M = np.log2(gt.loc[GENES10, tissues].astype(float).to_numpy().T + 1.0)   # 54 tissues x 10 genes
Mdf = pd.DataFrame(M, index=tissues, columns=GENES10)

# average-linkage clustering of rows (tissues) and columns (genes)
rl = linkage(Mdf.to_numpy(), method="average", metric="euclidean")
cl = linkage(Mdf.to_numpy().T, method="average", metric="euclidean")
rord = dendrogram(rl, no_plot=True)["leaves"]
cord = dendrogram(cl, no_plot=True)["leaves"]
Hc = Mdf.iloc[rord, cord]

# Tau tissue-specificity score and peak tissue per gene
tau, peak = {}, {}
for gene in GENES10:
    v = Mdf[gene].to_numpy()
    tau[gene] = float(np.sum(1 - v / v.max()) / (len(v) - 1))
    peak[gene] = short(Mdf.index[int(np.argmax(v))])
```

### FigA3 TCGA expression panorama and tumour-vs-normal effect sizes (build_A_figA3.py)
Data: tcga_patient_gene_log2tpm.csv, F1_tumour_vs_normal_12cohorts.csv — 9,523 patients, 31 cancer types (per-cancer medians); 12 paired cohorts with log2FC and BH-corrected q.

```python
pat = load("tcga_patient_gene_log2tpm.csv")
genes = GENES10
med = pat.groupby("cancer")[genes].median()          # cancer type x gene medians
n_per = pat["cancer"].value_counts()

# A6: cluster cancer-type columns (rows of med) by median expression profile
cl = linkage(med.to_numpy(), method="ward")
cord = dendrogram(cl, no_plot=True)["leaves"]
med = med.iloc[cord].T                               # (10 genes) x (31 cancer types)

# A7: tumour-vs-normal log2FC and q, genes x cohorts
f1 = load("F1_tumour_vs_normal_12cohorts.csv")
piv = f1.pivot(index="feature", columns="cohort", values="log2FC").loc[genes]
qpiv = f1.pivot(index="feature", columns="cohort", values="q").loc[genes]
# SKCM (n_normal=1) is all-NaN in the source file -> cluster valid cohorts only, SKCM greyed to last column
valid = piv.columns[piv.notna().all()]
skcm = [c for c in piv.columns if c not in valid]
sub = piv[valid]
rl = linkage(sub.to_numpy(), method="average")
cl = linkage(sub.to_numpy().T, method="average")
rord = dendrogram(rl, no_plot=True)["leaves"]
cord = dendrogram(cl, no_plot=True)["leaves"]
sub, qsub = sub.iloc[rord, cord], qpiv[valid].iloc[rord, cord]
piv = pd.concat([sub, piv[skcm].iloc[rord]], axis=1)
qpiv = pd.concat([qsub, qpiv[skcm].iloc[rord]], axis=1)
vlim = np.abs(sub.to_numpy()).max()
```

### FigA4 Volcano gallery and inverse-variance meta-analysis (build_A_figA4.py)
Data: F1_tumour_vs_normal_12cohorts.csv — 12 paired cohorts; the 6 largest by tumour n (BRCA/KIRC/HNSC/LUAD/THCA/LUSC) shown as volcanoes; SKCM (n_normal=1) excluded from the meta-analysis.

```python
f1 = load("F1_tumour_vs_normal_12cohorts.csv")
sizes = f1.groupby("cohort")["n_tumour"].first().sort_values(ascending=False)
top6 = list(sizes.index[:6])

# A8 volcano panels: one cohort per panel; points (log2FC, -log10(p)) coloured by gene;
# gene labelled when |log2FC| > 1 and p < 0.001 (reference lines at |log2FC|=1, p=0.001).

# A9 inverse-variance meta-analysis: SE back-derived from the 95% CI
se = (f1["ci_hi"] - f1["ci_lo"]) / 3.92
f1 = f1.assign(SE=se)
merged = {}
for g in GENES10:
    d = f1[(f1["feature"] == g) & f1["log2FC"].notna()]
    w = 1 / d["SE"] ** 2
    est = float(np.sum(d["log2FC"] * w) / np.sum(w))
    se_m = float(np.sqrt(1 / np.sum(w)))
    merged[g] = (est, est - 1.96 * se_m, est + 1.96 * se_m)
```

### FigA5 Group-score heatmap and per-cancer ridges (build_A_figA5.py)
Data: tcga_patient_gene_groupscore.csv — n=9,523 patients; group score = mean of member-gene log2TPM z-scores (S4 contract).

```python
gscore = load("tcga_patient_gene_groupscore.csv")
GROUP_COLS = GROUPS4                       # C1_Core_SMC56 / C2_Canonical / C3_SLF_repair / C4_SIMC
med = gscore.groupby("cancer")[list(GROUP_COLS)].median()   # cancer type x 4 groups
n_per = gscore["cancer"].value_counts()

# cluster cancer-type columns by score profile; symmetric colour limit
cl = linkage(med.to_numpy(), method="ward")
cord = dendrogram(cl, no_plot=True)["leaves"]
med = med.iloc[cord]
vlim = float(np.abs(med.to_numpy()).max())

# A11 ridges: top-9 cancer types by sample size, one KDE per group score
top9 = list(n_per.index[:9])
for gi, g in enumerate(GROUP_COLS):
    sub = gscore[gscore["cancer"].isin(top9)]
    lo = sub[g].min() - 0.3; hi = sub[g].max() + 0.3
    grid = np.linspace(lo, hi, 220)
    gap = 0.78
    for i, c in enumerate(reversed(top9)):            # larger cohorts plotted on top
        v = sub.loc[sub["cancer"] == c, g].to_numpy()
        if len(v) < 5:
            continue
        d = gaussian_kde(v)(grid)
        d = d / d.max()
```

### FigA6 Gene-gene correlation and score contribution (build_A_figA6.py)
Data: F1_gene_correlation_matrix.csv, F1_gene_correlation_pvalues.csv, F1_gene_score_contribution.csv, S4_sc_group_score_contract.csv — Pearson correlations on n=9,523 patients (log2TPM).

```python
R = load("F1_gene_correlation_matrix.csv").set_index("gene").loc[GENES10, GENES10]
P = load("F1_gene_correlation_pvalues.csv").set_index("gene").loc[GENES10, GENES10]
con = load("F1_gene_score_contribution.csv")
s4 = load("S4_sc_group_score_contract.csv")
mem = {g: set(str(s4.loc[s4["group"] == g, "members"].iloc[0]).split(" + ")) for g in GROUPS4}

# A12: lower-triangle Pearson correlation matrix (diagonal kept, significance stars from P)
n = len(GENES10)
M = R.to_numpy().astype(float)
Mshow = M.copy()
Mshow[np.triu_indices(n, 1)] = np.nan          # mask upper triangle

# A13: gene-to-group score contribution (Pearson r of gene expression vs its group score)
pv = con.pivot(index="gene", columns="group", values="pearson_r")
# member genes: dot size/colour = 90 + 620 * |r| on the correlation colormap; non-members blank
```

### FigA7 GTEx+TCGA context clustering and integrated panorama (build_A_figA7.py)
Data: GTEx_10gene_median_tpm.csv, tcga_patient_gene_log2tpm.csv, F1_tumour_vs_normal_12cohorts.csv, S1_dataset_eligibility_registry.csv — 85 contexts (54 GTEx tissues + 31 TCGA cancer types); registry covers 14 data sources x 10 genes.

```python
# A14: GTEx 54 tissues + TCGA 31 cancer types, median log2TPM
gt = load("GTEx_10gene_median_tpm.csv").set_index("Description").drop(columns=["Name"])
G = np.log2(gt.loc[genes, list(gt.columns)].astype(float) + 1.0).T     # 54 x 10
pat = load("tcga_patient_gene_log2tpm.csv")
T = pat.groupby("cancer")[genes].median()                              # 31 x 10
combo = pd.concat([G, T], axis=0)
labels = list(G.index) + list(T.index)
types = ["GTEx"] * len(G) + ["TCGA"] * len(T)

# ward clustering on correlation distance, rows reordered
D = pdist(combo.to_numpy(), metric="correlation")
Z = linkage(D, method="ward")
ordr = dendrogram(Z, no_plot=True)["leaves"]
combo = combo.iloc[ordr]
labels = [labels[i] for i in ordr]
types = [types[i] for i in ordr]

# A15-1 reproducibility: counts of cohorts with q<0.05, up/down per gene
f1 = load("F1_tumour_vs_normal_12cohorts.csv")
cohorts = list(pd.unique(f1["cohort"]))
sig_up, sig_dn, sign_m = {}, {}, {}
for g in genes:
    d = f1[f1["feature"] == g].set_index("cohort")
    sig_up[g] = int(((d["q"] < 0.05) & (d["log2FC"] > 0)).sum())
    sig_dn[g] = int(((d["q"] < 0.05) & (d["log2FC"] < 0)).sum())
    sign_m[g] = d.loc[cohorts, "log2FC"]

# A15-2 direction-consistency matrix: +1/-1 when q<0.05, +/-0.35 when not, NaN when n=1
M = np.zeros((10, len(cohorts)))
for i, g in enumerate(genes):
    for j, c in enumerate(cohorts):
        v = sign_m[g][c]
        q = f1[(f1["feature"] == g) & (f1["cohort"] == c)]["q"].iloc[0]
        if pd.isna(v):
            M[i, j] = np.nan
        elif q < 0.05:
            M[i, j] = 1 if v > 0 else -1
        else:
            M[i, j] = 0.35 if v > 0 else -0.35

# A15-3 data-availability barcode from the S1 eligibility registry
reg = load("S1_dataset_eligibility_registry.csv")
reg = reg[reg["gene"].isin(genes)].set_index("gene").loc[genes]
src_cols = list(reg.columns[:-1])                     # 14 data sources (all columns but status)
avail = reg[src_cols].apply(lambda c: c.astype(str).str.startswith("AVAILABLE"))
```

## Module B - Mutation and survival

### FigB1 Mutation and survival (build_B_figB1_prognosis.py)
Data: B1_mutation_prognosis_pancancer.csv, B3_residue_prognosis_pancancer.csv, B4_km_input_mut_vs_wt.csv, B6_biomarker_composite_ranking.csv — stratified logrank HR, n=598 with 198 events across 15 cancer strata.

```python
from scipy.stats import chi2 as chi2_dist

b1 = load("B1_mutation_prognosis_pancancer.csv")
b3 = load("B3_residue_prognosis_pancancer.csv")
b4 = load("B4_km_input_mut_vs_wt.csv")
b6 = load("B6_biomarker_composite_ranking.csv")

d = b1.dropna(subset=["HR"]).sort_values("HR", ascending=True).reset_index(drop=True)   # (a) forest
d3 = b3.dropna(subset=["HR"]).sort_values("HR", ascending=True).reset_index(drop=True)  # (c) residues

# (b) hand-written Kaplan-Meier estimator (t = time, e = event indicator)
def km_curve(t, e):
    t = np.asarray(t, float); e = np.asarray(e, int)
    uniq = np.unique(t[e == 1])
    s = 1.0; xs, ys = [0.0], [1.0]
    for u in uniq:
        at = (t >= u).sum(); dcnt = ((t == u) & (e == 1)).sum()
        s *= (1 - dcnt / at)
        xs += [u, u]; ys += [s, s]  # step: horizontal to u, vertical drop
    return np.array(xs), np.array(ys)

# (b) hand-written logrank test (mutated vs wild-type, chi2 with 1 df)
def logrank(t1, e1, t0, e0):
    t1 = np.asarray(t1, float); e1 = np.asarray(e1, int)
    t0 = np.asarray(t0, float); e0 = np.asarray(e0, int)
    times = np.unique(np.concatenate([t1[e1 == 1], t0[e0 == 1]]))
    O1 = E1 = V = 0.0
    for u in times:
        n1 = (t1 >= u).sum(); n0 = (t0 >= u).sum(); n = n1 + n0
        if n == 0: continue
        d1 = ((t1 == u) & (e1 == 1)).sum(); d0 = ((t0 == u) & (e0 == 1)).sum()
        dd = d1 + d0
        O1 += d1; E1 += dd * n1 / n
        if n > 1:
            V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    stat = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return chi2_dist.sf(stat, 1)

# representative KM combos: n_mut>=8 and >=4 events in the mutant arm, top 4 by event count
cand = []
for (g, c), dd in b4.groupby(["gene", "cancer"]):
    mut = dd[dd.mut_status == "mutated"]; wt = dd[dd.mut_status == "wild-type"]
    if len(mut) >= 8 and mut.event.sum() >= 4:
        cand.append((int(mut.event.sum()), g, c, mut, wt))
cand.sort(key=lambda x: -x[0])
kms = cand[:4]

# (d) composite biomarker ranking: ascending sort, bar colour keyed by tier (B/C/D)
d6 = b6.sort_values("composite_score", ascending=True)
```

## Module M - Mutation landscape

### FigM1 Pan-cancer mutation load (build_M_figM1_panorama.py)
Data: M7_gene_cancer_mut_burden.csv, M1_mutations_long.csv

```python
m7 = load("M7_gene_cancer_mut_burden.csv")
m1 = load("M1_mutations_long.csv")

# (a) gene x cancer somatic mutation burden heatmap
piv = m7.pivot_table(index="hugo_symbol", columns="cancer",
                     values="n_mutated_samples", aggfunc="sum").fillna(0)
order_g = piv.sum(1).sort_values(ascending=False).index
order_c = piv.sum(0).sort_values(ascending=False).index
piv = piv.loc[order_g, order_c]
# drop too-sparse cancers: keep leading cancers covering <= 0.99 of cumulative burden
cum = piv.sum(0).cumsum() / piv.sum(0).sum()
keep = cum[cum <= 0.99].index
piv = piv[keep]

# (b) per-gene mutation class composition, rebuilt from M1
# (M2_mutation_class_matrix.csv header columns lost the cancer prefix; n=921 mutation records)
sub = m1[m1.hugo_symbol.isin(GENES10)].copy()
sub["aa_event"] = sub["aa_event"].fillna("other")
cls_of = sub.groupby(["hugo_symbol", "aa_event"]).size().unstack(fill_value=0)
CLS = ["missense", "frameshift", "nonsense", "splice", "inframe", "other"]
cls_of = cls_of.reindex(columns=[c for c in CLS if c in cls_of.columns],
                        index=cls_of.sum(1).sort_values(ascending=False).index)
```

### FigM2 Mutation classes and lollipop (build_M_figM2_lollipop.py)
Data: M3_aa_lollipop_<gene>.csv (SMC6, SLF2, SMC5, SLF1), M4_aa_domain_annotation.csv

```python
GENES = ["SMC6", "SLF2", "SMC5", "SLF1"]  # top-4 genes by M7 mutation burden
loll = {g: load(f"M3_aa_lollipop_{g}.csv") for g in GENES}
dom = load("M4_aa_domain_annotation.csv")

for gene in GENES:
    d = loll[gene].copy()
    d["aa_pos"] = d["aa_pos"].astype(int)
    agg = d.groupby("aa_pos").agg(n=("n_mutations", "sum"),
                                  cls=("aa_event", "first")).reset_index()
    plen = int(dom.loc[dom.hugo_symbol == gene, "protein_length"].max())
    # classes can be mixed: colour each site by its majority class
    domcls = d.groupby(["aa_pos", "aa_event"]).size().reset_index(name="c") \
              .sort_values("c").drop_duplicates("aa_pos").set_index("aa_pos")["aa_event"]
    agg["cls"] = agg["aa_pos"].map(domcls)
    # pack annotated domains (feature_type != "Chain") into non-overlapping band rows
    dd = dom[(dom.hugo_symbol == gene) & (dom.feature_type != "Chain")].copy()
    dd = dd.sort_values(["start", "end"])
    rows, row_end = [], []
    for _, r in dd.iterrows():
        placed = False
        for i in range(len(rows)):
            if r.start > row_end[i]:
                rows[i].append(r); row_end[i] = r.end; placed = True; break
        if not placed:
            rows.append([r]); row_end.append(r.end)
    hot = agg.loc[agg.n.idxmax()]  # dominant hotspot: class, position, count
```

### FigM3 Domain enrichment and recurrent hotspots (build_M_figM3_domain_hotspot.py)
Data: M3b_domain_enrichment.csv, M3_aa_lollipop_all.csv

```python
m3b = load("M3b_domain_enrichment.csv")
loll = load("M3_aa_lollipop_all.csv")

# (a) domain-level enrichment of mutated residues (Fisher exact per UniProt feature,
#     odds_ratio precomputed upstream; BH q in source file, none reached q < 0.05)
d = m3b[m3b.description.notna()].copy()
d = d[d.n_mut_residues_in_domain >= 2]
d["label"] = d.hugo_symbol + " · " + d.description.str.slice(0, 26) + \
             " (" + d.start.astype(int).astype(str) + "-" + d.end.astype(int).astype(str) + ")"
d = d.sort_values("odds_ratio", ascending=True)
d = d.drop_duplicates("label")

# (b) cross-cancer recurrent hotspot residues
g = loll.groupby(["hugo_symbol", "aa_pos"]).agg(
    n=("patient_id", "size"), n_cancers=("cancer", "nunique")).reset_index()
mode_cls = loll.groupby(["hugo_symbol", "aa_pos"])["aa_event"] \
               .agg(lambda s: s.mode().iat[0]).rename("cls").reset_index()
g = g.merge(mode_cls)
g = g.sort_values(["n", "n_cancers"], ascending=False).head(14).iloc[::-1]
```

### FigM4 Hotspot panorama and amino-acid exchange (build_M_figM4_aa_hotspot_pano.py)
Data: M3_aa_lollipop_all.csv, M4_aa_domain_annotation.csv, G1_hotspot_prognosis.csv

```python
loll = load("M3_aa_lollipop_all.csv")
dom = load("M4_aa_domain_annotation.csv")
hot = load("G1_hotspot_prognosis.csv")

# amino-acid exchange parsed from protein_change (e.g. R123K -> alt "K"; stop = "*")
def parse_alt(pc):
    s = str(pc)
    m = re.match(r'^[A-Z](\d+)([A-Z\*])$', s)
    return m.group(2) if m else None

loll["alt_aa"] = loll.protein_change.map(parse_alt)
n_parsed = int(loll.alt_aa.notna().sum())
ex = loll.dropna(subset=["alt_aa"]).groupby(["ref_aa", "alt_aa"]).size() \
         .sort_values(ascending=False)
ex_top = ex.head(12).iloc[::-1]
ref_top = loll.ref_aa.value_counts()  # most-mutated reference amino acids

# cross-cancer recurrent hotspots (top 12)
hs = loll.groupby(["hugo_symbol", "aa_pos"]).agg(
    n=("patient_id", "size"), n_cancers=("cancer", "nunique")).reset_index()
hs = hs.sort_values(["n", "n_cancers"], ascending=False).head(12)
per = loll[loll.set_index(["hugo_symbol", "aa_pos"]).index.isin(
    pd.MultiIndex.from_frame(hs[["hugo_symbol", "aa_pos"]]))] \
    .groupby(["hugo_symbol", "aa_pos", "cancer"]).size().reset_index(name="cnt")
cls_map = loll.groupby(["hugo_symbol", "aa_pos"])["aa_event"] \
              .agg(lambda s: s.mode().iat[0]).to_dict()

# hotspot OS Cox summary (computed upstream); carriers < 5 are marked underpowered
hot_ok = hot[hot.status == "OK"].copy()
n_hot_all = len(hot)
n_hot_up = int((hot.status != "OK").sum())
hot_ok = hot_ok.sort_values(["HR"], ascending=True).reset_index(drop=True)
```

## Module C - Copy-number variation

### FigC1 CNA distribution (build_CNA_figC1.py)
Data: M5_cna_matrix.csv, M6_cna_summary.csv

```python
mat = pd.read_csv(os.path.join(DATA, "M5_cna_matrix.csv"), index_col=0)
m6 = pd.read_csv(os.path.join(DATA, "M6_cna_summary.csv"))

# gene rows grouped by complex: core pair -> canonical -> SLF -> SIMC/EID3
gene_order = ["SMC5", "SMC6", "SMC1A", "SMC3", "SMC2", "SMC4",
              "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
              "SLF1", "SLF2", "SIMC1", "EID3"]
mat = mat.reindex(gene_order)

# cancer columns sorted by mean GISTIC CNA (descending) to highlight amplification bias
col_order = mat.mean(axis=0).sort_values(ascending=False).index
mat = mat[col_order]

# (b) pan-cancer CNA category fractions, n=10,845 TCGA samples
m6 = m6.set_index("hugo_symbol").reindex(gene_order[::-1])  # bottom-to-top = gene_order
CN_COLS = ["pct_amplification", "pct_gain", "pct_shallow_deletion", "pct_deep_deletion"]
pct = m6[CN_COLS] * 100  # percent of pan-cancer samples
```

### FigC2 CNA-expression dose response (build_CNA_figC2.py)
Data: M9_cna_expression_paired_long.csv, M9_cna_expression_dose_response.csv

```python
GENES6 = ["SMC5", "SMC6", "NSMCE2", "NSMCE3", "SLF1", "SLF2"]
CATS = [(-2.0, "Deep del", PAL["slate"]),
        (-1.0, "Het loss", PAL["aqua"]),
        (0.0,  "Neutral",  PAL["tan"]),
        (1.0,  "Gain",     PAL["peach"]),
        (2.0,  "Amp",      PAL["red"])]

m9 = pd.read_csv(os.path.join(DATA, "M9_cna_expression_paired_long.csv"))
m9 = m9[m9["hugo_symbol"].isin(GENES6)]
dose = pd.read_csv(os.path.join(DATA, "M9_cna_expression_dose_response.csv")).set_index("hugo_symbol")

# (a) CN category x expression distributions with embedded quartile lines (categories with n > 5)
for g in GENES6:
    sub = m9[m9["hugo_symbol"] == g]
    data = [sub.loc[sub["cna"] == c, "toil_log2_rsem"].values for c, _, _ in CATS]
    for i, d in enumerate(data):
        if len(d) > 5:
            q1, med, q3 = np.percentile(d, [25, 50, 75])
    rho, q = dose.loc[g, "rho"], dose.loc[g, "q"]  # Spearman rho + BH q from M9

# (b) within-sample CNA-expression OLS fit, log2(TPM) ~ CNA (derived statistic)
rng = np.random.default_rng(56)
for g in GENES6:
    sub = m9[m9["hugo_symbol"] == g]
    x = sub["cna"].values + rng.uniform(-0.16, 0.16, len(sub))  # jitter for display only
    y = sub["toil_log2_rsem"].values
    b1, b0 = np.polyfit(sub["cna"].values, y, 1)
```

### FigC3 CNA prognosis (build_CNA_figC3.py)
Data: G1_cna_prognosis_cox.csv, G1_cna_prognosis_stratified.csv

```python
cox = pd.read_csv(os.path.join(DATA, "G1_cna_prognosis_cox.csv"))
st = pd.read_csv(os.path.join(DATA, "G1_cna_prognosis_stratified.csv"))

# (a) per gene x cancer Cox (CNA high vs low): powered fits only, BH q < 0.05
sig = cox[(cox["status"] == "OK") & (cox["q"] < 0.05)].copy()
# group rows by GENES10 order; within group sort by HR ascending
gorder = [g for g in GENES10 if g in sig["hugo_symbol"].unique()]
sig["g"] = sig["hugo_symbol"].map({g: i for i, g in enumerate(gorder)})
sig = sig.sort_values(["g", "HR"]).reset_index(drop=True)

# (b) pan-cancer stratified Cox (strata = cancer), HR per unit GISTIC, n=10,622, 3,446 events
st = st.set_index("hugo_symbol").reindex(GENES10).reset_index()
```

## Module CL - Clinical prognosis

### FigCL1 Univariate vs multivariable Cox (build_CL1_clinical.py)
Data: F8_pancancer_univariate_cox.csv, F8_pancancer_multivariable_cox.csv, F8_clinical_molecular_incremental_cindex.csv

```python
ux = load("F8_pancancer_univariate_cox.csv")
mx = load("F8_pancancer_multivariable_cox.csv")
ci = load("F8_clinical_molecular_incremental_cindex.csv")

sig = ux[ux.q < 0.05].copy()
sig["pair"] = sig.cancer + " · " + sig.feature
# cluster rows by feature (GROUPS4 first), within cluster by HR descending
forder = GROUPS4 + ["SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
                    "EID3", "SLF1", "SLF2", "SIMC1"]
sig["_f"] = sig.feature.map({f: i for i, f in enumerate(forder)})
sig = sig.sort_values(["_f", "HR"], ascending=[True, False]).reset_index(drop=True)

# (b) match the same significant pairs to the ridge-adjusted multivariable Cox
mkey = mx.dropna(subset=["HR"]).set_index(["cancer", "feature"])
rows_b, miss = [], []
for _, r in sig.iterrows():
    if (r.cancer, r.feature) in mkey.index:
        m = mkey.loc[(r.cancer, r.feature)]
        r2 = r.copy()
        r2["HR"], r2["lo"], r2["hi"] = m["HR"], m["lo"], m["hi"]
        r2["q"], r2["direction"] = m["q"], m["direction"]
        rows_b.append(r2)
    else:
        miss.append(f"{norm_cancer(r.cancer)}·{r.feature}")  # 7 uncovered pairs (ACC x6, PRAD x1)
rows_b = pd.DataFrame(rows_b)
rows_b["_f"] = rows_b.feature.map({f: i for i, f in enumerate(forder)})
rows_b = rows_b.sort_values(["_f", "HR"], ascending=[True, False])

# (c) molecular-only Harrell C-index (21 powered cancers; clinical-only arms underpowered upstream)
ok = ci[(ci.status == "OK") & (ci.model == "molecular_only")].sort_values("c_index")
```

### FigCL2 Gene-set scores and patient survival (build_CL2_group_patient_prognosis.py)
Data: tcga_patient_gene_groupscore.csv, F8_pancancer_univariate_cox.csv, F8_pancancer_multivariable_cox.csv, B4_km_input_mut_vs_wt.csv

```python
G4 = ["C1_Core_SMC56", "C2_Canonical", "C3_SLF_repair", "C4_SIMC"]

gs = load("tcga_patient_gene_groupscore.csv")
uni = load("F8_pancancer_univariate_cox.csv")
mvar = load("F8_pancancer_multivariable_cox.csv")
km_in = load("B4_km_input_mut_vs_wt.csv")

# (a) patient-level cohort overview
cnt = gs.cancer.value_counts()
n_total = int(cnt.sum())
os_pat = km_in.groupby("patient").agg(cancer=("cancer", "first"),
                                      time=("time", "first"), event=("event", "first")).reset_index()
os_cancers = set(os_pat.cancer.str.upper())
matched = gs[gs.patient.isin(set(os_pat.patient))]
n_os = len(os_pat)
n_matched = len(matched)
n_ev = int(os_pat.event.sum())
n_mv = len(mvar[mvar.feature.isin(G4) & (mvar.status == "OK")])
n_mv_sig = int((mvar[mvar.feature.isin(G4)].q < 0.05).sum())

# (b) univariate Cox rows for the four complex scores (HR from published Cox summaries)
ug = uni[uni.feature.isin(G4) & (uni.status == "OK")].copy()
n_uni_ok = len(ug)
sig = ug[ug.p < 0.05].sort_values("p", ascending=True).reset_index(drop=True)

# (c) KM: patient-level OS merged with the four complex scores
os_pat["cancer_u"] = os_pat.cancer.str.upper()
mm = os_pat.merge(gs[["patient", "cancer"] + G4], on="patient", suffixes=("", "_gs"))
# take the 5 largest cancers with events (UCEC dropped: 0 events)
ev_by_c = mm.groupby("cancer_u").event.sum().sort_values(ascending=False)
use_cancers = [c for c in ev_by_c[ev_by_c > 0].index][:5]
km_sub = mm[mm.cancer_u.isin(use_cancers)].copy()

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

for feat in G4:
    # high/low = median split of the complex score within each cancer, then pooled
    rows = []
    for c in use_cancers:
        sub = km_sub[km_sub.cancer_u == c]
        med = sub[feat].median()
        rows.append((sub[sub[feat] > med], sub[sub[feat] <= med]))
    hi_df = pd.concat([r[0] for r in rows])
    lo_df = pd.concat([r[1] for r in rows])
    for df_, lab in [(lo_df, f"low half (n={len(lo_df)}, events={int(lo_df.event.sum())})"),
                     (hi_df, f"high half (n={len(hi_df)}, events={int(hi_df.event.sum())})")]:
        kmf = KaplanMeierFitter()
        kmf.fit(df_.time, df_.event, label=lab)
    lr = logrank_test(hi_df.time, lo_df.time, hi_df.event, lo_df.event)
    p = lr.p_value
```

## Module T - Tumour microenvironment and pathways

### FigT1 TME signature associations (build_T1_tme.py)
Data: F7_tme_signature_matrix_wide.csv, F7_tme_feature_signature_correlation.csv, F7_tme_checkpoint_association.csv

```python
wide = load("F7_tme_signature_matrix_wide.csv")
corr = load("F7_tme_feature_signature_correlation.csv")
cp = load("F7_tme_checkpoint_association.csv").drop_duplicates(["group", "checkpoint"])

SIGS = ["CD8_T", "Cytotoxicity", "IFNG_response", "MHC_I", "T_exhaustion"]
CKPT = ["CD274", "CTLA4", "HAVCR2", "IDO1", "LAG3", "PDCD1", "TIGIT"]
GORDER = GROUPS4 + ["SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
                    "EID3", "SLF1", "SLF2", "SIMC1"]

# (a) cancer x signature medians of patient-level ssGSEA scores (n=6,260)
med = wide.groupby("cancer")[SIGS].median()
med = med.sort_values("CD8_T", ascending=False)

# (b) feature x immune-checkpoint Spearman rho (BH q for stars; duplicate PDCD1 rows deduplicated)
pv = cp.pivot_table(index="group", columns="checkpoint", values="rho").reindex(index=GORDER, columns=CKPT)
qs = cp.pivot_table(index="group", columns="checkpoint", values="q").reindex(index=GORDER, columns=CKPT)

# (c) per-cancer group x signature Spearman rho (20 cancers, status OK)
cancers20 = sorted(corr.cancer.unique())
for g in GROUPS4:
    sub = corr[(corr.group == g)
               & (corr.signature.isin(SIGS))
               & (corr.status == "OK")]
    pvt = sub.pivot_table(index="cancer", columns="signature", values="rho").reindex(index=cancers20, columns=SIGS)
    qst = sub.pivot_table(index="cancer", columns="signature", values="q").reindex(index=cancers20, columns=SIGS)
```

### FigG GSEA after virtual knockout (build_VP_figG.py)
Data: P6_gsea_hallmark_matrix.csv, P6_gsea_all.csv, P6_gsea_group.csv, P6_gsea_summary.csv, V12_vp_gsea_landscape.csv, s6a_circularity.tsv

```python
# sys.path also points at the sibling MODULE_A_REBUILD_20261011/scripts dir (provides smc56_style, DATA)
hm   = load("P6_gsea_hallmark_matrix.csv")
allg = load("P6_gsea_all.csv")
grpg = load("P6_gsea_group.csv")
summ = load("P6_gsea_summary.csv")
v12  = load("V12_vp_gsea_landscape.csv")
cir  = pd.read_csv(os.path.join(DATA, "s6a_circularity.tsv"), sep="\t")

# (a) Hallmark NES matrix (13 features x 40 terms); FDR looked up by (feature_key, term)
pathways = [c for c in hm.columns if c != "feature_key"]
feats = hm["feature_key"].tolist()
nes = hm.set_index("feature_key")[pathways].astype(float)
fdr_lookup = {}
for df in (allg, grpg):
    for _, r in df.iterrows():
        fdr_lookup[(r["feature_key"], r["term"])] = r["FDR"]

# (c) post-virtual-KO GSEA landscape (9 pathways x 14 features), ordered by max |NES|
pw_order = (v12.groupby("pathway")["NES"].apply(lambda s: s.abs().max())
            .sort_values(ascending=False).index.tolist())
ft_order = (v12.assign(a=v12["NES"].abs()).groupby("feature")["a"].max()
            .sort_values(ascending=False).index.tolist())
m12 = v12.pivot(index="feature", columns="pathway", values="NES")\
         .reindex(index=ft_order, columns=pw_order)
f12 = v12.pivot(index="feature", columns="pathway", values="fdr")\
         .reindex(index=ft_order, columns=pw_order)

# (d) circularity audit: score-expression association vs cell-composition confounding
c = cir.reset_index(drop=True)
raw = c[c["analysis"] == "raw_spearman"].reset_index(drop=True)  # Spearman rho, BH-corrected FDR
kw  = c[c["analysis"] == "celltype_kw"]                          # between-group Kruskal-Wallis
mvn = c[c["analysis"].str.startswith("mal_vs_nonmal")]           # malignant vs non-malignant
mvn_max_f = mvn["fdr"].max()  # all p > 0.05: no circularity after deconvolution
```

### FigF1 Activity and SI scores (build_F1_activity_si.py)
Data: F1_activity_si_patient_level.csv, F1_cancer_feature_summary.csv

```python
f1 = load("F1_activity_si_patient_level.csv")
summ = load("F1_cancer_feature_summary.csv")

cancers = (f1.groupby("cancer").size().sort_values(ascending=False).index.tolist())
order = (f1.groupby("cancer").Activity.median().sort_values(ascending=False).index.tolist())

# (a) per-cancer ridgeline densities: Activity z and log10 SI (SI clipped at 1e-3)
n = len(order)
dy = 1.0
scale = 0.82  # ridge peak height (overlapping ridges)
act_vals = [f1.loc[f1.cancer == c, "Activity"] for c in order]
si_vals = [np.log10(f1.loc[f1.cancer == c, "SI"].clip(lower=1e-3)) for c in order]

def ridge(ax, values_by_cancer, color, edge, xlab, xticks):
    lo = min(v.min() for v in values_by_cancer)
    hi = max(v.max() for v in values_by_cancer)
    pad = 0.06 * (hi - lo)
    xs = np.linspace(lo - pad, hi + pad, 200)
    for i, c in enumerate(order):
        v = values_by_cancer[i].dropna().values
        kde = gaussian_kde(v)
        d = kde(xs)
        d = d / d.max() * scale
        y0 = (n - 1 - i) * dy
        # per-cancer KDE rows drawn as overlapping filled ridges (plot calls omitted)

# (b) feature means z-scored across cancers (rows = 10 genes + 4 complex scores)
feat_order = ["SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3",
              "SLF1", "SLF2", "SIMC1", "C1_Core_SMC56", "C2_Canonical",
              "C3_SLF_repair", "C4_SIMC"]
mat = summ.pivot_table(index="feature", columns="cancer", values="mean").reindex(feat_order)
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1, ddof=1), axis=0)
col_order = z.loc["C1_Core_SMC56"].sort_values(ascending=False).index.tolist()
z = z[col_order]

# (c) Activity vs SI patient-level scatter with global Spearman correlation
other = f1[~f1.cancer.isin(["BRCA", "LUAD", "COADREAD", "HNSC", "KIRC", "LGG"])]
hl = ["BRCA", "COADREAD", "HNSC", "KIRC", "LUAD", "LGG"]
rho, pv = spearmanr(f1.Activity, f1.SI)  # n=9523 patients, 31 cancers
```

## Module ML - Machine learning and dependency

### FigML1 Model benchmark (build_FigML_model_perf.py)
Data: P2_coadread_nestedcv_benchmark.csv, P2_coadread_calibration.csv, P2_coadread_dca.csv, P2_coadread_learning_curve.csv

```python
bench = load("P2_coadread_nestedcv_benchmark.csv")
calib = load("P2_coadread_calibration.csv")
dca   = load("P2_coadread_dca.csv")
lc    = load("P2_coadread_learning_curve.csv")
# (a) nested 5x4 patient-level CV benchmark (seed 42, locked splits), OOF AUC; n=227, 84 events
b = bench.sort_values("auc", ascending=True).reset_index(drop=True)
best_i = int(b["auc"].idxmax())
# (b) calibration: overall OOF bins only (horizon NA), 10 bins per model
cal = calib[calib["horizon"].isna()]
show = [("ElasticNet-logit", "clinical_only"), ("ElasticNet-logit", "molecular_only"),
        ("CoxBoost-GB", "clinical_plus_molecular"), ("Survival-SVM", "molecular_only")]
for m, fs in show:
    d = cal[(cal.model == m) & (cal.feature_set == fs)].sort_values("pred_mean")
# (c) DCA: net_benefit vs threshold; net_benefit_all = treat-all reference curve
d0 = dca[(dca.model == "ElasticNet-logit") & (dca.feature_set == "clinical_only")]
# (d) learning curve: aggregate the 5 repeats per training fraction
g = lc.groupby("frac").agg(n_train=("n_train", "first"),
                           tr_m=("train_auc", "mean"), tr_s=("train_auc", "std"),
                           te_m=("test_auc", "mean"), te_s=("test_auc", "std")).reset_index()
# training AUC identically 1.0 at every fraction (overfitting flag); validation AUC flat
```

### FigML1b Calibration, DCA, learning curves (build_FigML1b_calibration_dca.py)
Data: P2_coadread_calibration.csv, P2_coadread_dca.csv, P2_coadread_learning_curve.csv

```python
calib = load("P2_coadread_calibration.csv")
dca   = load("P2_coadread_dca.csv")
lc    = load("P2_coadread_learning_curve.csv")
# upstream (P2_manuscript 4.3): glmnet elastic-net/LASSO-Cox, CoxBoost, survival RF/XGBoost, SuperLearner;
# preprocessing frozen on training folds; calibration/DCA/learning curves on held-out patients only.
# (a) calibration: overall OOF bins (horizon NA); molecular models conservative at high risk
cal = calib[calib["horizon"].isna()]
for m, fs in show:
    d = cal[(cal.model == m) & (cal.feature_set == fs)].sort_values("pred_mean")
d0 = dca[(dca.model == "ElasticNet-logit") & (dca.feature_set == "clinical_only")]  # (b) DCA, n=227 / 84 events
# (c) learning curve: validation AUC 0.517-0.571 regardless of n (same groupby agg as FigML1d)
```

### FigML2 Feature importance and stability (build_FigML_feature_robust.py)
Data: P2_coadread_feature_importance.csv, P2_coadread_bootstrap_selection_frequency.csv, P2_coadread_seed_stability.csv

```python
fi   = load("P2_coadread_feature_importance.csv")
boot = load("P2_coadread_bootstrap_selection_frequency.csv")
stab = load("P2_coadread_seed_stability.csv")
GROUP_FEATS = {"C1_Core_SMC56", "C2_Canonical", "C3_SLF_repair", "C4_SIMC"}  # complex-group feats
# (a) permutation importance: mean AUC drop (+/-SD, 30 permutations), CoxBoost-GB molecular model
d = fi[fi.model == "CoxBoost-GB"].sort_values("mean_drop_auc", ascending=True)  # NSMCE2/NSMCE3/SLF2 dominate
# (b) bootstrap selection frequency: % selected over 300 resamples; C4_SIMC / SIMC1 most stable
d = boot.sort_values("selection_freq", ascending=True)
# (c) cross-seed stability: 20 seeds per feature, per-run delta-AUC jittered around the mean
order = (stab.groupby("feature")["drop_auc"].mean()
         .sort_values(ascending=True).index.tolist())
rng = np.random.default_rng(7)
for i, f in enumerate(order):
    v = stab[stab.feature == f]["drop_auc"].values
    # scatter(v, i + rng.uniform(-0.16, 0.16, len(v))); vline at np.mean(v); negative = noise
```

### FigML4 DepMap dependency panorama (build_FigML4_depmap_panorama.py)
Data: P2_depmap_lineage_dependency.csv, P2_depmap_dependency_summary.csv, P2_depmap_homolog_comparator.csv

```python
lin  = load("P2_depmap_lineage_dependency.csv")
summ = load("P2_depmap_dependency_summary.csv")
homo = load("P2_depmap_homolog_comparator.csv")
# DepMap 24Q4 Chronos (<0 = fitness loss on inhibition), 1,178 cell lines, 25 lineages;
# homo flags SMC1A/SMC3/SMC2/SMC4 as not_in_DepMap24Q4_release -> honest missing, no fabricated values.
# lineage x gene matrices; pct_common_dependent = % of lines with chronos < -0.5 (from source file)
piv_m = lin.pivot_table(index="lineage", columns="gene", values="mean_chr", aggfunc="mean")
piv_p = lin.pivot_table(index="lineage", columns="gene", values="pct_common_dependent", aggfunc="mean")
genes = [g for g in GENES10 if g in piv_m.columns]
piv_m = piv_m[genes].dropna()
piv_p = piv_p[genes].reindex(piv_m.index)
order = piv_m["SMC6"].sort_values().index.tolist()  # strongest SMC6 dependency first
piv_m = piv_m.loc[order]; piv_p = piv_p.loc[order]
n_lines = lin.groupby("lineage")["n_lines"].max().reindex(order)
# (c) member ranking: 1,178-line mean bar + per-lineage mean scatter; pct_below_minus05 = % < -0.5
g = summ.sort_values("mean_chr", ascending=False).reset_index(drop=True)  # SMC6: -0.85, 93.5% lines < -0.5
```

### FigML5 Immuno-oncology external validation (build_FigML5_immuno.py)
Data: P2_imvigor210_response_auc.csv, P2_imvigor210_tmb_cox.csv

```python
resp = load("P2_imvigor210_response_auc.csv")
cox  = load("P2_imvigor210_tmb_cox.csv")
# IMvigor210 (atezolizumab): response n=298 (68 resp / 230 non-resp); OS TMB-adjusted Cox n=242 (159 events); precomputed scores.
# (a) response AUC, 95% CI as given in the source table (CI method not stated there)
r = resp.sort_values("AUC", ascending=True).reset_index(drop=True)
# all AUC 0.508-0.580, every CI crosses 0.5 (negative result, reported as-is)
# (b) TMB-adjusted Cox forest per module score (Core-SMC / Canonical-NSE / SLF-repair / SIMC-targeting)
c = cox.iloc[::-1].reset_index(drop=True)
for i, rr in enumerate(c.itertuples()):
    star = p_stars(rr.q)   # q = BH-FDR; all q >= 0.34, HR 1.01-1.18
```

### FigML3 DepMap + immuno-oncology combined (build_FigML_depmap_immuno.py)
Data: P2_depmap_dependency_with_lineage.csv, P2_depmap_dependency_summary.csv, P2_imvigor210_response_auc.csv, P2_imvigor210_tmb_cox.csv

```python
dep  = load("P2_depmap_dependency_with_lineage.csv")
summ = load("P2_depmap_dependency_summary.csv")
resp = load("P2_imvigor210_response_auc.csv")
cox  = load("P2_imvigor210_tmb_cox.csv")
# (a) gene x lineage mean chronos over 30 lineages, rows sorted by SMC6 dependency strength
piv = dep.pivot_table(index="OncotreeLineage", columns="gene", values="chronos", aggfunc="mean")
genes = [g for g in GENES10 if g in piv.columns]
piv = piv[genes].dropna()
order = piv["SMC6"].sort_values().index.tolist()  # sort by SMC6 dependency
piv = piv.loc[order]
n_cell = dep.groupby("OncotreeLineage")["ModelID"].nunique()
# (b) member comparison: 1,178-line mean bar + per-lineage mean scatter
g = summ.sort_values("mean_chr", ascending=False)
lineage_mean = dep.pivot_table(index="OncotreeLineage", columns="gene", values="chronos", aggfunc="mean")
# comparators SMC1A/SMC3/SMC2/SMC4 absent from DepMap 24Q4 -> intra-complex comparison only
```

## Module SC - Single-cell atlas

### SC helpers (sc_common.py)
Data: F3_sc_cell_level.csv

```python
# CLUSTER_DEF: 27 display groups as (cluster_label, celltype, [Cell_subtype members]) over 6 major
# cell types (T 7 / Epithelial 6 / B 3 / Myeloids 4 / Stromal 6 / Mast 1); e.g. ("T cells · CD8+ cell", "T cells", ["CD8+ T cells"]).
# Other shared helpers: build_palette, kde_underlay (gaussian_kde underlay), ondata_labels, legend_n.
def assign_clusters(df):
    lab = pd.Series(index=df.index, dtype=object)
    sub = df["Cell_subtype"].astype(str)
    ct = df["celltype"].astype(str)
    for label, ctype, subs in CLUSTER_DEF:
        m = (ct == ctype) & sub.isin(subs)
        lab[m] = label
    # Unknown cells: KMeans(k=2) blocks in UMAP space, then 30-NN majority vote from known subtypes
    for ctype in ("T cells", "Myeloids", "B cells"):
        u_mask = (ct == ctype) & (sub == "Unknown")
        if not u_mask.any():
            continue
        k_mask = (ct == ctype) & lab.notna()
        known = df.loc[k_mask]
        unk = df.loc[u_mask]
        km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(unk[["umap_1", "umap_2"]])
        nn = NearestNeighbors(n_neighbors=30).fit(known[["umap_1", "umap_2"]].values)
        for blk in (0, 1):
            b_idx = unk.index[km.labels_ == blk]
            if len(b_idx) == 0:
                continue
            _, ind = nn.kneighbors(df.loc[b_idx, ["umap_1", "umap_2"]].values)
            votes = known["Cell_subtype"].values[ind]          # (n_blk, 30)
            for i, cidx in enumerate(b_idx):
                winner = pd.Series(votes[i]).mode()[0]
                for label, ctype2, subs in CLUSTER_DEF:
                    if ctype2 == ctype and winner in subs:
                        lab[cidx] = label
                        break
    return lab
def cluster_counts(ctab):
    return ctab["cluster"].value_counts().reindex(CLUSTER_ORDER)
```

### FigSC1 Clustering and annotation (build_sc_fig1_full.py)
Data: F3_sc_cell_level.csv

```python
df = load("F3_sc_cell_level.csv")
assert len(df) == 61024, "cell count mismatches delivery contract"  # n=61,024 cells, 36 subtypes
cnt = df["Cell_subtype"].value_counts()
fam_of = df.groupby("Cell_subtype")["celltype"].first().to_dict()
groups = list(cnt.index)  # 36 Cell_subtypes, descending abundance; no merging, Unknown (n=2,181) kept
n_tumor = int((df["Class"] == "Tumor").sum()); n_normal = int((df["Class"] == "Normal").sum())
n_pat = df["patient"].nunique()
# family color = blend_white(FAM_BASE[fam], 0.75 * i / max(k - 1, 1)) by descending abundance rank
# (blend_white mixes the hex color toward white by t; COLOR["Unknown"] = "#A7A29E")
# density cloud recipe: histogram2d -> gaussian_filter(sigma=2.6) -> normalized density field
XMIN, XMAX = df["umap_1"].min() - 1.5, df["umap_1"].max() + 1.5
YMIN, YMAX = df["umap_2"].min() - 1.5, df["umap_2"].max() + 1.5
GX, GY = 420, 300
XX_RANGE = [XMIN, XMAX]; YY_RANGE = [YMIN, YMAX]; SIGMA = 2.6
def density_field(sub):
    z, _, _ = np.histogram2d(sub["umap_1"], sub["umap_2"], bins=[GX, GY],
                             range=[XX_RANGE, YY_RANGE])
    z = gaussian_filter(z.astype(float), SIGMA)
    z = z.T  # histogram2d axis0=x -> transpose (rows=y, cols=x)
    z /= (z.max() + 1e-9)
    return z
# draw_cloud(ax, sub, color, alpha_peak=0.62, floor=0.10): alpha = clip((z - floor) / (1 - floor), 0, 1) ** 0.75
# * alpha_peak, zero below floor; imshow(rgba, extent=[XMIN, XMAX, YMIN, YMAX], origin="lower")
big = [s for s in groups if cnt[s] >= 1000]         # groups labeled on the cloud
draw_order = sorted(groups, key=lambda s: -cnt[s])  # large groups underneath, small on top
```

### FigSC2 Cell composition (build_sc_fig2_composition.py)
Data: F3_sc_cell_level.csv (via sc_common.get_cluster_table cache: sc_cell_clusters.csv)

```python
ctab = get_cluster_table()
ct_order = ["Epithelial cells", "T cells", "B cells", "Myeloids", "Stromal cells", "Mast cells"]
both = [p for p in ctab.patient.unique() if (ctab[ctab.patient == p].Class == "Normal").any()]
tum_only = [p for p in ctab.patient.unique() if p not in both]
def epi_frac(p):
    g = ctab[ctab.patient == p]
    return (g.celltype == "Epithelial cells").mean()
both.sort(key=epi_frac); tum_only.sort(key=epi_frac)
pat_order = both + tum_only
# (a) per-patient major-type composition (%)
comp = (ctab.groupby(["patient", "celltype"]).size().unstack(fill_value=0)
        .reindex(pat_order)[ct_order])
comp_pct = comp.div(comp.sum(1), axis=0) * 100
# (b) display-cluster x patient abundance, log10(n+1); 27 clusters x 23 patients
mat = (ctab.groupby(["cluster", "patient"]).size().unstack(fill_value=0)
       .reindex(index=CLUSTER_ORDER, columns=pat_order))
mat_log = np.log10(mat + 1)
ph_order = ["G1", "S", "G2M"]  # (c) cell-cycle phase share per major type
ph = (ctab.groupby(["celltype", "phase"]).size().unstack(fill_value=0)
      .reindex(index=ct_order, columns=ph_order))
ph_pct = ph.div(ph.sum(1), axis=0) * 100
```

### FigSC3 Module scores (build_sc_fig3_modules.py)
Data: sc_cell_clusters.csv (cache built from F3_sc_cell_level.csv)

```python
df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "_cache", "sc_cell_clusters.csv"))
ct_order = ["Epithelial cells", "T cells", "B cells", "Myeloids", "Stromal cells", "Mast cells"]
SCORES = [("Core", "C1_Core_SMC56", "C1 Core · SMC5/6 core pair"),
          ("Canonical", "C2_Canonical", "C2 Canonical · canonical six subunits"),
          ("SLF-repair", "C3_SLF_repair", "C3 SLF-repair · repair complex"),
          ("SIMC", "C4_SIMC", "C4 SIMC · SIMC1 interactome")]
# (a-d) feature plots: symmetric color scale, vmax = 95th percentile of |z| (not 99),
# gamma=0.45 two-sided power stretch (near-zero contrast gain ~3x)
GAMMA = 0.45
vlims = []
for pi, (feat, gkey, ttl) in enumerate(SCORES):
    vmax = float(np.percentile(np.abs(df[feat]), 95))   # 95th pct, symmetric about 0
    vlims.append(vmax)
    # soft_feature_image + feature_dots on umap_1/umap_2, CMAP_DIV, [-vmax, vmax], gamma=GAMMA
# (e) violins per major type x module (width 0.21); overlay q1, med, q3 = np.percentile(d, [25, 50, 75])
# module score floors (source-data plateau): Core -0.31 / Canonical -0.36 / SLF-repair -0.18 / SIMC -0.16
```

### FigSC4 Pseudotime and trajectories (build_sc_fig4_pseudotime.py)
Data: sc_cell_clusters.csv (cache built from F3_sc_cell_level.csv)

```python
df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "_cache", "sc_cell_clusters.csv"))
ct_order = ["Epithelial cells", "T cells", "B cells", "Myeloids", "Stromal cells", "Mast cells"]
SCORES = [("Core", "C1_Core_SMC56", "C1 Core"), ("Canonical", "C2_Canonical", "C2 Canonical"),
          ("SLF-repair", "C3_SLF_repair", "C3 SLF-repair"), ("SIMC", "C4_SIMC", "C4 SIMC")]
vmax_pt = 0.40  # (a) pseudotime UMAP color clip at the 99th pct (~1% terminal mast cells above)
# (b) ridgeline per major type: xs = np.linspace(0, 1, 240);
#     kde = gaussian_kde(v, bw_method=0.12)(xs), fill to y0 + kde / kde.max() * 0.82
# (c) module score vs pseudotime: 40 equal-width bins, bin means + rolling mean (window 7, shrinking edges)
bins = np.linspace(0, 1, 41)
mid = (bins[:-1] + bins[1:]) / 2
for name, key, short in SCORES:
    cut = pd.cut(df["pseudotime"], bins, labels=False)
    m = np.array([df.loc[cut == i, name].mean() for i in range(len(mid))])
    sm = np.copy(m)
    for i in range(len(m)):
        lo, hi = max(0, i - 3), min(len(m), i + 4)
        seg = m[lo:hi]
        sm[i] = seg[~np.isnan(seg)].mean() if (~np.isnan(seg)).any() else np.nan
```

### FigSC5 Subtype-level analysis (build_sc_fig5_analysis.py)
Data: F3_sc_cell_level.csv

```python
df = load("F3_sc_cell_level.csv")
submap = {...}  # display names for 10 key subtypes, e.g. "CD8+ T cells": "T · CD8+", "SPP1+": "Myeloid · SPP1+"
df["disp"] = df["Cell_subtype"].map(lambda s: submap.get(s, s))
cnt = df["disp"].value_counts()
fam_of = df.groupby("Cell_subtype")["celltype"].first().to_dict()
fam_short = {"T cells": "T", "Epithelial cells": "Epi", "B cells": "B",
             "Myeloids": "Myeloid", "Stromal cells": "Stromal", "Mast cells": "Mast"}
rare = [k for k, v in cnt.items() if v < 600 and k != "Unknown"]
rare_map = {k: f"{fam_short.get(fam_of.get(k, '?'), '?')} · other" for k in rare}
df["plot_group"] = df["disp"].mask(df["disp"].isin(rare), df.loc[df["disp"].isin(rare), "Cell_subtype"].map(rare_map))
SCORES = [("Core", "C1_Core_SMC56"), ("Canonical", "C2_Canonical"),
          ("SLF-repair", "C3_SLF_repair"), ("SIMC", "C4_SIMC")]
# detection floor per module = mode of the raw score (source-data plateau)
FLOOR = {c: float(df[c].value_counts().idxmax()) for c, _ in SCORES}
for c, _ in SCORES:
    df["act_" + c] = df[c] > FLOOR[c] + 1e-6
# per-group aggregation: counts, module means, per-gene means, % cells with gene z > 0
agg = df.groupby("plot_group").agg(
    n=("barcode", "size"),
    **{c: (c, "mean") for c, _ in SCORES},
    **{g: (g, "mean") for g in GENES10},
)
for g in GENES10:
    agg[g + "_pos"] = df.assign(pos=df[g] > 0).groupby("plot_group")["pos"].mean()
# tumor/normal abundance log2 ratio with pseudocount 0.1
piv_t = df[df["Class"] == "Tumor"]["plot_group"].value_counts()
piv_n = df[df["Class"] == "Normal"]["plot_group"].value_counts()
agg["log2ratio"] = np.log2((piv_t.reindex(agg.index).fillna(0) + 0.1) /
                           (piv_n.reindex(agg.index).fillna(0) + 0.1))
# module detection T vs N: Fisher exact (two-sided) on 2x2 counts, BH within each module
nT = (df["Class"] == "Tumor").groupby(df.loc[df["Class"] == "Tumor", "plot_group"]).size().reindex(agg.index).fillna(0)
nN = (df["Class"] == "Normal").groupby(df.loc[df["Class"] == "Normal", "plot_group"]).size().reindex(agg.index).fillna(0)
for c, _ in SCORES:
    aT = df[df["Class"] == "Tumor"].groupby("plot_group")["act_" + c].mean().reindex(agg.index)
    aN = df[df["Class"] == "Normal"].groupby("plot_group")["act_" + c].mean().reindex(agg.index)
    ps = []
    for k in agg.index:
        vt, vn = aT.get(k, np.nan), aN.get(k, np.nan)
        if not (np.isfinite(vt) and np.isfinite(vn)) or nT[k] < 20 or nN[k] < 20:
            ps.append(np.nan); continue
        a = int(round(vt * nT[k])); b = int(nT[k]) - a
        cc = int(round(vn * nN[k])); d = int(nN[k]) - cc
        ps.append(fisher_exact([[a, b], [cc, d]], alternative="two-sided")[1])
    pv = pd.Series(ps, index=agg.index)
    valid = pv.notna().values
    adj = np.full(len(pv), np.nan)
    mvs = pv.values[valid]
    oi = np.argsort(mvs)
    run = np.minimum.accumulate((mvs[oi] * len(mvs) / (np.arange(len(mvs)) + 1))[::-1])[::-1]
    tmp = np.empty(len(mvs)); tmp[oi] = np.minimum(run, 1)
    adj[valid] = tmp
    agg["q_" + c] = adj
    agg["actT_" + c] = aT
    agg["actN_" + c] = aN
agg = agg.sort_values("Core", ascending=False)  # shared row order for all panels
```

### FigSC6 Ten-subunit expression panorama (build_sc_fig6_genes.py)
Data: sc_cell_clusters.csv, F3_sc_cell_level.csv

```python
df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "_cache", "sc_cell_clusters.csv"))
if "SMC5" not in df.columns:  # cache lacks gene columns: merge from source by barcode
    src = load("F3_sc_cell_level.csv")
    df = df.merge(src[["barcode"] + GENES10], on="barcode", how="left")
    assert df["SMC5"].notna().all(), "missing gene values after barcode merge"
# module membership badges (S4 contract): C1=SMC5+SMC6; C2=SMC5/6+NSMCE1-4A; C3=SLF1+SLF2; C4=SIMC1+SLF2
TAG = {"C1_Core_SMC56": "C1", "C2_Canonical": "C2",
       "C3_SLF_repair": "C3", "C4_SIMC": "C4"}
det = {}
for gi, g in enumerate(GENES10):
    vmax = float(np.percentile(df[g], 99)) or 1.0   # symmetric color clip at the 99th pct
    # soft_feature_image + feature_dots on umap_1/umap_2, CMAP_DIV, [-vmax, vmax]
    det[g] = float((df[g] > 0).mean())              # positive rate (z > 0)
```

## Module V - Virtual perturbation and pharmacology

### FigV1 Knockout downstream and rewiring (build_VP_figV1.py)
Data: V1_vp_ko_downstream.csv, V2_vp_rewiring_summary.csv, V2b_vp_edgeloss_long.csv

```python
v1 = load("V1_vp_ko_downstream.csv")
v2 = load("V2_vp_rewiring_summary.csv")
v2b = load("V2b_vp_edgeloss_long.csv")
# (a) top18 downstream genes by max |Z| across subunits, self edges excluded (n=1,556 genes x 10 KO)
sub_order = list(v2["subunit"])  # same subunit order as V2
nons = v1[~v1.is_self]
top_genes = (nons.groupby("gene")["Z"].max().sort_values(ascending=False)
             .head(18).index.tolist())
mat = (v1[v1.gene.isin(top_genes)]
       .pivot(index="gene", columns="subunit", values="log2FC")
       .reindex(index=top_genes, columns=sub_order))
logfc_c = mat.clip(-5, 5)   # log2FC (KO vs control) clipped to +/-5 for display
# (b) rewiring summary: Spearman rho of wiring strength vs perturbation-response Z, per subunit
v2s = v2.sort_values("rho_wiring_vs_Z")
# n_signif = significantly rewired genes per subunit (p_adj < 0.05)
# (c) edge loss vs gain per subunit: delta = edge_out - edge_in, summed by sign
v2b["delta"] = v2b["edge_out"] - v2b["edge_in"]
agg = (v2b.groupby("subunit")
       .apply(lambda g: pd.Series({
           "lost": -g.loc[g["delta"] < 0, "delta"].sum(),
           "gained": g.loc[g["delta"] > 0, "delta"].sum()}), include_groups=False)
       .reindex(sub_order))   # NSMCE3/SLF2: no rewiring (raw values ~ 0)
```

### FigV2 Fate shifts and TF sensitivity (build_VP_figV2.py)
Data: V6_vp_fate_shift.csv, V5_vp_tf_sensitivity.csv, V6b_vp_lineage_vulnerability.csv, V4_vp_vectorfield_qa.csv

```python
v6  = load("V6_vp_fate_shift.csv")
v5  = load("V5_vp_tf_sensitivity.csv")
v6b = load("V6b_vp_lineage_vulnerability.csv")
v4  = load("V4_vp_vectorfield_qa.csv")
# (a) fate shift: mean projection onto the CMS transition field, 9 subunits x 3 CMS subtypes
# (Epithelial lineage; n=2,487-4,366 cells per panel); frac = fraction shifting toward CMS2
sub_order = ["SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
             "SLF1", "SLF2", "SIMC1"]
cms_order = ["CMS1", "CMS2", "CMS3"]
mat = (v6.pivot(index="subunit", columns="subtype", values="mean_proj_cms")
       .reindex(index=sub_order, columns=cms_order))
frac = (v6.pivot(index="subunit", columns="subtype", values="frac_toward_CMS2")
        .reindex(index=sub_order, columns=cms_order))
vlim = np.nanmax(np.abs(mat.values))
# (b) TF sensitivity: strongest-|Z| KO condition per TF, top 12 of 62 TFs; negative Z enriched
v5["absZ"] = v5["Z"].abs()
top_tf = (v5.sort_values("absZ", ascending=False)
          .drop_duplicates("tf").head(12).sort_values("Z"))
# (c) lineage vulnerability: mean relative perturbation magnitude, 9 subunits x 6 lineages
lin_order = ["B cells", "Epithelial cells", "Mast cells",
             "Myeloids", "Stromal cells", "T cells"]
mat2 = (v6b.pivot(index="subunit", columns="lineage", values="mean_rel_mag")
        .reindex(index=sub_order, columns=lin_order))
vmax = np.nanmax(mat2.values)   # n_cells per grid cell 44-4,366; B cells and epithelium most sensitive
v4s = v4.sort_values("ratio")  # (d) QA: observed vs permutation null magnitude, ratio > 1 for all
```

### FigV3 Pharmacology and synthetic lethality (build_VP_figV3.py)
Data: V7_vp_ddr_drug.csv, V7b_vp_gdsc_full.csv, V8_vp_synthetic_lethal.csv, V8b_vp_sl_null_record.csv

```python
v7  = load("V7_vp_ddr_drug.csv")
v7b = load("V7b_vp_gdsc_full.csv")
v8  = load("V8_vp_synthetic_lethal.csv")
v8b = load("V8b_vp_sl_null_record.csv")
# (a) DDR drug response: Spearman rho of expression vs drug IC50 (rho < 0 = higher expression,
# greater sensitivity); top 20 drugs x 13 features; V7 = FDR-filtered subset of the 574
# GDSC associations in V7b (all q < 0.05 here; cell labels = p_stars(q))
top_drugs = (v7.groupby("drug")["rho"].apply(lambda s: s.abs().max())
             .sort_values(ascending=False).head(20).index.tolist())
feat_order = (v7.groupby("feature")["rho"].apply(lambda s: s.abs().max())
              .sort_values(ascending=False).index.tolist())
mat = (v7[v7.drug.isin(top_drugs)]
       .pivot(index="drug", columns="feature", values="rho")
       .reindex(index=top_drugs, columns=feat_order))
qv = (v7[v7.drug.isin(top_drugs)]
      .pivot(index="drug", columns="feature", values="q")
      .reindex(index=top_drugs, columns=feat_order))
vlim = np.nanmax(np.abs(mat.values))
# (b) synthetic lethal candidates: top 24 of 403 pairs by |depmap_r| (DepMap co-dependency r);
# null check: gene-level DepMap x KO-rewiring intersection = 0 -> no candidate passes the null
v8s = v8.sort_values("abs_depmap_r", ascending=False).head(24).iloc[::-1]
null_inter = int(v8b.loc[v8b["key"] == "total_intersecting_genes", "value"].iloc[0])
# (c) combined top hits: V8 top6 co-dependency partners + V7 top6 drugs by |rho|
v7t = v7.reindex(v7["rho"].abs().sort_values(ascending=False).index).drop_duplicates("drug").head(6)
v8t = v8.sort_values("abs_depmap_r", ascending=False).drop_duplicates("candidate").head(6)
```

### FigV4 Multi-method KO consensus (build_VP_figV4_consensus.py)
Data: V10_vp_layer_overlap.csv, V1_vp_ko_downstream.csv, V2b_vp_edgeloss_long.csv, V5_vp_tf_sensitivity.csv, P6_gsea_all.csv, GTEx_10gene_median_tpm.csv, tcga_patient_gene_log2tpm.csv

```python
v10 = load("V10_vp_layer_overlap.csv")
v1 = load("V1_vp_ko_downstream.csv")
v2b = load("V2b_vp_edgeloss_long.csv")
v5 = load("V5_vp_tf_sensitivity.csv")
p6 = load("P6_gsea_all.csv")
gtex = load("GTEx_10gene_median_tpm.csv")
tcga = load("tcga_patient_gene_log2tpm.csv")
# (a) cross-method overlap: DepMap top40 co-dependency partners x KO rewired genes
# (hypergeometric chance expectation); observed overlap = 0 for all 10 subunits
d = v10.sort_values("expected_by_chance", ascending=True)
med_exp = float(d["expected_by_chance"].median())   # median chance expectation ~14.9, fold = 0
# (b) within-KO consensus: V1 significant downstream (p_adj < 0.05, self excluded)
#     intersected with V2b lost edges (edge_out < edge_in)
v2b["delta"] = v2b["edge_out"] - v2b["edge_in"]
rows = []
for sub, g in v1.groupby("subunit"):
    sig = g[(~g.is_self) & (g.p_adj < 0.05)]
    lost = set(v2b[(v2b.subunit == sub) & (v2b.delta < 0)]["gene"])
    inter = sig[sig.gene.isin(lost)].sort_values("Z", ascending=False)
    rows.append(dict(subunit=sub, n_v1=len(sig), n_cons=len(inter),
                     genes=";".join(inter["gene"])))
cons = pd.DataFrame(rows).sort_values(["n_cons", "n_v1"], ascending=True).reset_index(drop=True)
# 5 unique consensus genes: BST2/ANXA5/IFI6/TIMP1/C19orf33
# (c) upstream TF -> subunit -> downstream triples (TF = V5 top3 by rank; downstream = V1 top3 by Z)
sub_order = GENES10
for i, sub in enumerate(sub_order):
    tf = v5[v5.subunit == sub].sort_values("rank").head(3)["tf"].tolist()
    dn = v1[(v1.subunit == sub) & (~v1.is_self)].sort_values("Z", ascending=False).head(3)
    # TF ordering is directional only (all p_adj > 0.05, does not pass FDR)
# (e) reference tissue baseline: GTEx (54 normal tissues, median TPM) vs TCGA (9,523 tumors, log2TPM)
tiss = [c for c in gtex.columns if c not in ("Name", "Description")]
gm = gtex.set_index("Description")[tiss].astype(float).T          # tissue x gene
gz = (gm - gm.mean()) / gm.std()
gprof = gz.median()                                               # cross-tissue median z per gene
tc = tcga[GENES10]
tz = (tc - tc.mean()) / tc.std()
tprof = tz.median()
rho, rho_p = spearmanr(gprof[GENES10], tprof[GENES10])
cmp = pd.DataFrame({"gtex": gprof, "tcga": tprof})
cmp["diff"] = cmp["tcga"] - cmp["gtex"]
cmp = cmp.sort_values("diff", ascending=True)
# z-normalized within the 10-gene panel (platform scales differ; no absolute comparison)
# (d) GSEA top terms: FDR < 0.05 ranked by |NES| (GO_BP top10 / KEGG top8 of 2065/366 rows)
VMAX = 5.5   # symmetric NES color scale; NES < 0 = downregulated after KO
for lib, ntop in (("GO_BP", 10), ("KEGG", 8)):
    d = p6[(p6.library_canon == lib) & (p6.FDR < 0.05)].copy()
    d["absN"] = d["NES"].abs()
    d = d.sort_values("absN", ascending=False).head(ntop)  # labels: term stripped of "(GO:\d+)", p_stars(FDR)
```
