import networkit as nk
import numpy as np
import pandas as pd
import time
import os
import powerlaw
from networkit import vizbridges
import argparse
from scipy.stats import spearmanr
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

valid_aa = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def log(msg):
    """
    Logging function. Log all the print outputs
    """
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(str(msg) + "\n")


def load_junction_aa_mapping(MAP_FILE):

    log("Reading the map file.")
    start = time.time()

    id_to_cdr3 = {}
    cdr3_to_id = {}

    with open(MAP_FILE) as f:
        for line in f:
            node_id_str, cdr3 = line.rstrip("\n").split("\t")
            node_id = int(node_id_str)

            id_to_cdr3[node_id] = cdr3
            cdr3_to_id[cdr3] = node_id

    log(f"Loaded {len(cdr3_to_id):,} nodes")
    log(f"Mapping read time: {time.time() - start:.2f} sec\n")

    return cdr3_to_id, id_to_cdr3


def load_graph(GRAPH_FILE):
    """
    Load graph from binary format
    """
    start = time.time()
    g = nk.readGraph(GRAPH_FILE, nk.Format.NetworkitBinary)
    log(f"Time to read the binary graph file: {time.time() - start:.2f} sec\n")
    return g

def load_epitope_info(EPITOPE_INFO_FILE):
    """
    Load Epitope information
    """
    df = pd.read_parquet(EPITOPE_INFO_FILE)
    df = df[df.akc_measurement_category != 'Negative']
    df = df[df["junction_aa"].str.match(valid_aa) & df["akc_epitope_seq_aa"].str.match(valid_aa)]
    df = df.drop_duplicates()
    log(f"Epitop dataframe shape: {df.shape}")
    return df

def create_bipartite_mapping(epitope_df):
    
    epitope_to_cdr3 = {}
    cdr3_to_epitope = {}

    for _, row in epitope_df.iterrows():
        cdr3 = row['junction_aa']
        epitope = row['akc_epitope_seq_aa']
        
        if epitope not in epitope_to_cdr3:
            epitope_to_cdr3[epitope] = []
        if cdr3 not in cdr3_to_epitope:
            cdr3_to_epitope[cdr3] = []
        
        epitope_to_cdr3[epitope].append(cdr3)
        cdr3_to_epitope[cdr3].append(epitope)
    
    return epitope_to_cdr3, cdr3_to_epitope


def get_receptor_count_per_epitope(G, epitope_to_cdr3, cdr3_to_id, id_to_cdr3):

    epitope_receptor_count = {}
    # Iterate through each epitope
    for epitope, cdr3_nodes in epitope_to_cdr3.items():
        # Set to store unique receptors (including distance=1 neighbors)
        receptors = set()
        # Add CDR3 nodes themselves (distance 0)
        receptors.update(cdr3_nodes)
        # For each CDR3 node, get its distance-1 neighbors (connected nodes)
        for cdr3 in cdr3_nodes:
            node = cdr3_to_id[cdr3]  # Map the CDR3 sequence to node ID
            for neighbor in G.iterNeighbors(node):
                neighbor_cdr3 = id_to_cdr3[neighbor]  # convert to CDR3 string
                receptors.add(neighbor_cdr3)
        # Store the number of unique receptors (including neighbors) for this epitope
        epitope_receptor_count[epitope] = len(receptors)

    return epitope_receptor_count

def get_epitope_distribution_per_receptor(G, cdr3_to_epitope, cdr3_to_id, id_to_cdr3):

    receptor_epitope_distribution = {}

    for cdr3, epitopes in cdr3_to_epitope.items():
        node = cdr3_to_id[cdr3]
        receptor_epitopes = set(epitopes)
        
        for neighbor in G.iterNeighbors(node):
            neighbor_cdr3 = id_to_cdr3.get(neighbor)
            
            if neighbor_cdr3 is None:
                continue
            neighbor_epitopes = cdr3_to_epitope.get(neighbor_cdr3)
            if neighbor_epitopes:
                receptor_epitopes.update(neighbor_epitopes)
        receptor_epitope_distribution[cdr3] = len(receptor_epitopes)

    return receptor_epitope_distribution

# def plot_inferred_specificity(epitope_specificity_df, receptor_specificity_df):
#     # sns.set(style="whitegrid")
#     sns.set(style="whitegrid")
#     fig, axes = plt.subplots(2, 2, figsize=(16, 12))
#     # -------------------------
#     # Epitope propagation
#     # -------------------------
#     sns.scatterplot(data=epitope_specificity_df, x="real_receptors", y="inferred_receptors", hue="new_receptors", palette="viridis", alpha=0.7, ax=axes[0,0])

#     max_val = epitope_specificity_df["inferred_receptors"].max()
#     axes[0,0].plot([0, max_val], [0, max_val], linestyle="--", color="gray")
#     axes[0,0].set_title("Epitope Specificity Propagation")
#     axes[0,0].set_xlabel("Real receptors")
#     axes[0,0].set_ylabel("Inferred(total) receptors")
#     axes[0,0].set_xscale("log")
#     axes[0,0].set_yscale("log")
#     # -------------------------
#     # Receptor propagation
#     # -------------------------
#     sns.scatterplot(data=receptor_specificity_df,x="real_epitopes", y="inferred_epitopes",hue="new_epitopes",palette="viridis",alpha=0.7,ax=axes[0,1],legend=True)
#     max_val = receptor_specificity_df["inferred_epitopes"].max()
#     axes[0,1].plot([0, max_val], [0, max_val], linestyle="--", color="gray")
#     axes[0,1].set_title("Receptor Specificity Propagation")
#     axes[0,1].set_xlabel("Real epitopes")
#     axes[0,1].set_ylabel("Inferred(total) epitopes")
#     # -------------------------
#     # Distribution of new receptors
#     # -------------------------
#     sns.histplot(epitope_specificity_df["new_receptors"],bins=70,ax=axes[1,0],color="steelblue")
#     axes[1,0].set_title("Distribution of New Receptors per Epitope")
#     axes[1,0].set_xlabel("New receptors inferred")
#     axes[1,0].set_ylabel("Count")
#     axes[1,0].set_xscale("log")
#     axes[1,0].set_yscale("log")
#     # -------------------------
#     # Distribution of new epitopes
#     # -------------------------
#     sns.histplot(receptor_specificity_df["new_epitopes"],bins=70,ax=axes[1,1],color="darkorange")
#     axes[1,1].set_title("Distribution of New Epitopes per Receptor")
#     axes[1,1].set_xlabel("New epitopes inferred")
#     axes[1,1].set_ylabel("Count")
#     axes[1,1].set_yscale("log")

#     plt.tight_layout()
#     plt.savefig(f"{FIG_DIR}/log_real_vs_inferred_specificity_distribution.png", dpi=600)
#     plt.close()
    
def plot_inferred_specificity(epitope_specificity_df, receptor_specificity_df):
    # sns.set(style="whitegrid")
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    # -------------------------
    # Epitope propagation
    # -------------------------
    sns.scatterplot(data=epitope_specificity_df, x="real_receptors", y="inferred_receptors", hue="new_receptors", palette="viridis", alpha=0.7, ax=axes)

    max_val = epitope_specificity_df["inferred_receptors"].max()
    axes.plot([0, max_val], [0, max_val], linestyle="--", color="gray")
    axes.set_xlabel("Observed CDR3")
    axes.set_ylabel("Inferred(total) CDR3")
    axes.set_xscale("log")
    axes.set_yscale("log")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/receptor_real_vs_inferred_specificity_distribution.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    # -------------------------
    # Receptor propagation
    # -------------------------
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    sns.scatterplot(data=receptor_specificity_df,x="real_epitopes", y="inferred_epitopes",hue="new_epitopes",palette="viridis",alpha=0.7,ax=axes,legend=True)
    max_val = receptor_specificity_df["inferred_epitopes"].max()
    axes.plot([0, max_val], [0, max_val], linestyle="--", color="gray")
    axes.set_xlabel("Observed epitopes")
    axes.set_ylabel("Inferred(total) epitopes")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/epitope_real_vs_inferred_specificity_distribution.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    
    
def epitope_df_figures(df):
    
    log(f"Total unique junction_aa: {df.junction_aa.nunique()}")
    log(f"Total unique epitope: {df.akc_epitope_seq_aa.nunique()}")
    
    # ----------------------------
    log(f"\nDegree per CDR3: ")
    # ----------------------------
    cdr3_degree = df.groupby("junction_aa")["akc_epitope_seq_aa"].nunique()
    log(cdr3_degree.describe())
    
    # ----------------------------
    log(f"\nDegree per Epitope: ")
    # ----------------------------
    
    epitope_degree = df.groupby("akc_epitope_seq_aa")["junction_aa"].nunique()
    log(epitope_degree.describe())
    
    log(f"\nCross-Reactivity Analysis")
    cross_reactive = cdr3_degree[cdr3_degree > 1]
    log(f"Total cross reactive CDR3 that has degree > 1: {len(cross_reactive)}")
    
    # ----------------------------
    # Plot CDR3 degree distribution
    # ----------------------------
    
    plt.figure(figsize=(5, 4))
    plt.hist(cdr3_degree, bins=50)
    plt.xlabel("Number of Unique Epitopes per CDR3")
    plt.ylabel("Number of CDR3")
    # plt.title("Epitope Specificity Distribution per Receptor")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/junction_aa_epitope_distribution.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    
    
    plt.figure(figsize=(5, 4))
    plt.hist(cdr3_degree, bins=50)
    plt.yscale("log")
    plt.xlabel("Number of Unique Epitopes per CDR3")
    plt.ylabel("Number of CDR3 (log scale)")
    # plt.title("Epitope Specificity Distribution ")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/junction_aa_epitope_distribution_log.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    
    # ----------------------------
    # Plot epitope degree distribution
    # ----------------------------
    plt.figure(figsize=(5, 4))
    plt.hist(epitope_degree, bins=50)
    plt.xlabel("Number of Unique CDR3 per Epitopes")
    plt.ylabel("Number of Epitopes")
    # plt.title("Number of Receptors Recognizing Each Epitope")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/epitope_receptor_distribution.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    
    
    plt.figure(figsize=(5, 4))
    plt.hist(epitope_degree, bins=50)
    plt.yscale("log")
    plt.xlabel("Number of Unique CDR3 per Epitope")
    plt.ylabel("Number of Epitopes (log scale)")
    # plt.title("Number of Receptors Recognizing Each Epitope")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/epitope_receptor_distribution_log.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    
    # ----------------------------
    # Select top nodes
    # ----------------------------
    top_epitopes = epitope_degree.sort_values(ascending=False).head(8).index
    top_cdr3 = cdr3_degree.sort_values(ascending=False).head(30).index
    
    # ----------------------------
    # Plot top CDR3 length and their binding information
    # ----------------------------
    sub_df = df[df["akc_epitope_seq_aa"].isin(top_epitopes)]
    # ----------------------------
    # Plot top epitopes CDR3 length distribution.
    # ----------------------------
    top6 = top_epitopes[:2]
    sub_df = df[df["akc_epitope_seq_aa"].isin(top_epitopes)].copy()
    sub_df = sub_df.drop_duplicates(["akc_epitope_seq_aa", "junction_aa"])
    sub_df["junction_aa_length"] = sub_df["junction_aa"].str.len()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    axes = axes.flatten()
    for i, epitope in enumerate(top6):
        ep_df = sub_df[sub_df["akc_epitope_seq_aa"] == epitope]
        median_len = ep_df["junction_aa_length"].median()
        sns.histplot(
            ep_df["junction_aa_length"],
            # bins=range(sub_df["junction_aa_length"].min(),sub_df["junction_aa_length"].max() + 2),
            ax=axes[i],
            # kde=True
        )
        axes[i].axvline(median_len, color="red", label=f"Median = {median_len}", linestyle="--")
        axes[i].set_title(f"{epitope} ({len(ep_df)})")
        axes[i].set_xlabel("CDR3 Length")
        axes[i].set_ylabel("Count")
        axes[i].legend()

    # plt.suptitle("CDR3 Length Distribution for Top 5 Epitopes", y=1.02)
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/top_n_epitope_junction_aa_length_distribution.png", bbox_inches = 'tight', dpi=300)
    plt.close()
    # ----------------------------
    # Plot top 10 CDR3's epitope length distribution.
    # ----------------------------
    top_n_cdr3 = top_cdr3[:2]
    sub_df = df[df["junction_aa"].isin(top_n_cdr3)].copy()
    sub_df = sub_df.drop_duplicates(["akc_epitope_seq_aa", "junction_aa"])
    sub_df["epitope_length"] = sub_df["akc_epitope_seq_aa"].str.len()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4),)
    axes = axes.flatten()
    for i, junction_aa in enumerate(top_n_cdr3):
        ep_df = sub_df[sub_df["junction_aa"] == junction_aa]
        median_len = ep_df["epitope_length"].median()
        sns.histplot(
            ep_df["epitope_length"],
            # bins=range(sub_df["epitope_length"].min(), sub_df["epitope_length"].max() + 2),
            ax=axes[i],
            # kde=True
        )
        axes[i].axvline(median_len, color="red", label=f"Median = {median_len}", linestyle="--")
        axes[i].set_title(f"{junction_aa} ({len(ep_df)})")
        axes[i].set_xlabel("Epitope Length")
        axes[i].set_ylabel("Count")
        axes[i].legend()

    # plt.suptitle("CDR3 Length Distribution for Top 10 CDR3's", y=1.02)
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/top_n_cdr3s_epitope_length_distribution.png", bbox_inches = 'tight', dpi=300)
    plt.close()
    

    n = 30
    log(f"\nTop {n} CDR3 with most epitopes: ")
    log(f"{cdr3_degree.sort_values(ascending=False).head(n)}")
    
    log(f"\nTop {n} Epitopes with most receptors: ")
    log(f"{epitope_degree.sort_values(ascending=False).head(n)}")
    
    sub_df = df[(df["akc_epitope_seq_aa"].isin(top_epitopes)) &(df["junction_aa"].isin(top_cdr3))]
    # sub_df = sub_df.drop_duplicates(["akc_epitope_seq_aa", "junction_aa"])
    log(sub_df.head())
    # ----------------------------
    # Create binary matrix
    # ----------------------------
    matrix = pd.crosstab(sub_df["akc_epitope_seq_aa"],sub_df["junction_aa"])
    # Ensure consistent ordering
    matrix = matrix.reindex(index=top_epitopes, columns=top_cdr3, fill_value=0)
    matrix = (matrix > 0).astype(int)
     
    matrix = matrix.loc[(matrix != 0).any(axis=1), (matrix != 0).any(axis=0)]

    # Plot heatmap
    plt.figure(figsize=(5, 4))
    sns.heatmap(matrix, cmap="viridis",
                cbar=False, 
                xticklabels=True, 
                yticklabels=True,
                linewidths=0.5, 
                linecolor='white')

    plt.ylabel("Epitope")
    plt.xlabel("CDR3")
    plt.tight_layout()
    # plt.title("Epitope–CDR3 Interaction (Filtered for Common Values)")
    plt.savefig(f"{FIG_DIR}/epitope_cdr3_heatmap.png", bbox_inches = 'tight', dpi=600)
    plt.close()
    # # linkage methodmetric="euclidean",    # distance metriccbar_kws={"label": "Interaction Count"}
    # g = sns.clustermap(matrix,cmap="viridis",figsize=(16, 12),xticklabels=True,yticklabels=True,method="average", )

    # # g.fig.suptitle("Epitope–CDR3 Interaction Clustermap\nTop 50 Epitopes vs Top 200 CDR3", y=1.02)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

    # plt.savefig(f"{FIG_DIR}/epitope_cdr3_cluster_heatmap.png", dpi=600)
    # plt.close() 
    
    # =========================
    # TOP CDR3 ANALYSIS
    # =========================

    top_cdr3 = ( cdr3_degree.sort_values(ascending=False).index)
    cdr3_df = df[df["junction_aa"].isin(top_cdr3)].copy()

    # Remove duplicate interactions
    cdr3_df = cdr3_df.drop_duplicates(["junction_aa", "akc_epitope_seq_aa"])
    # Compute CDR3 length
    cdr3_df["junction_aa_length"] = cdr3_df["junction_aa"].str.len()
    # Count number of epitopes per CDR3
    cdr3_summary = (cdr3_df.groupby(["junction_aa", "junction_aa_length"])["akc_epitope_seq_aa"].nunique().reset_index(name="epitope_count"))

    # =========================
    # TOP EPITOPE ANALYSIS
    # =========================

    epitope_degree = (df.drop_duplicates(["junction_aa", "akc_epitope_seq_aa"]).groupby("akc_epitope_seq_aa")["junction_aa"].nunique())

    top_epitopes = (epitope_degree.sort_values(ascending=False).index)

    epitope_df = df[df["akc_epitope_seq_aa"].isin(top_epitopes)].copy()
    epitope_df = epitope_df.drop_duplicates(["junction_aa", "akc_epitope_seq_aa"])

    # Compute epitope length
    epitope_df["epitope_length"] = epitope_df["akc_epitope_seq_aa"].str.len()

    # Count number of CDR3 per epitope
    epitope_summary = (epitope_df.groupby(["akc_epitope_seq_aa", "epitope_length"])["junction_aa"].nunique().reset_index(name="cdr3_count"))
    log(epitope_summary.head())
    # =========================
    #  Correlation Stats
    # =========================

    rho1, p1 = spearmanr(cdr3_summary["junction_aa_length"],cdr3_summary["epitope_count"])

    rho2, p2 = spearmanr(epitope_summary["epitope_length"],epitope_summary["cdr3_count"])
    
    # =========================
    # PLOTTING
    # =========================

    # sns.set_style("whitegrid")

    fig, axes = plt.subplots(1, 1, figsize=(5, 4))

    # ---- CDR3 perspective ----
    sns.scatterplot(data=cdr3_summary, x="junction_aa_length", y="epitope_count", ax=axes)

    axes.set_xlabel("CDR3 Length")
    axes.set_ylabel("Number of Epitopes Bound")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/junction_aa_length_and_binding_distribution.png", bbox_inches= 'tight', dpi=600)
    plt.close()
    # ---- Epitope perspective ----
    
    
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    sns.scatterplot(data=epitope_summary,x="epitope_length",y="cdr3_count",ax=axes)

    axes.set_xlabel("Epitope Length (aa)")
    axes.set_ylabel("Number of Unique CDR3")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/epitope_length_and_binding_distribution.png", bbox_inches= 'tight', dpi=600)
    plt.close()

    log(f"CDR3 perspective Spearman cdr3 length vs epitope count: {rho1} p = {p1}")
    log(f"Epitope perspective Spearman epitope length vs cdr3 count: {rho2} p = {p2}")
    
def create_directories_if_not_exist(path):
    """Create directories if they do not exist"""
    if not os.path.exists(path):
        os.makedirs(path)
        
LOG_FILE = "run.log"
FIG_DIR = "Figures/"
def main():
    parser = argparse.ArgumentParser("CDR3 Graph and Epitope Analysis Parameters.")
    parser.add_argument("--LOCUS", default='trb', help="Locus to build/analyze the graph (tra, trb, trd, trg)")
    parser.add_argument("--VERSION", default="v1", help="Version of the table name that will be put on the graph")
    parser.add_argument("--DATA_DIR", default="/ak_graph_data/", help="Version of the table name that will be put on the graph")
    parser.add_argument("--MAX_THREADS", default=64, help="Maximum number of threads to use by the graph algorithm")
    
    args = parser.parse_args()
    LOCUS = args.LOCUS
    VERSION = args.VERSION
    DATA_DIR = args.DATA_DIR
    MAX_THREADS = args.MAX_THREADS
    
    global LOG_FILE
    global FIG_DIR
    
    if LOCUS not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code

    MAP_FILE = f"{DATA_DIR}/pair_files/{LOCUS}_output_seq_map_{VERSION}.tsv"
    GRAPH_FILE = f"{DATA_DIR}/graph_files/{LOCUS}_graph_{VERSION}.nkbg003"
    LCC_FILE = f'{DATA_DIR}/graph_files/{LOCUS}_largest_connected_component_{VERSION}.edgelist.txt'
    FIG_DIR = f"{DATA_DIR}/figures/{LOCUS}_{VERSION}/epitope_analysis/"
    LOG_FILE = f"{DATA_DIR}/logs/{LOCUS}_cdr3_epitope_analysis_{VERSION}.log"
    EPITOPE_INFO_FILE = f"{DATA_DIR}/epitope_data/{LOCUS}_cdr3_epitope_info_{VERSION}.parquet"
    # local_fig_dir = f'./ak_graph_data_v2/figures/{LOCUS}_{VERSION}'
    
    # Create the necessary directories if they do not exist
    create_directories_if_not_exist(f"{FIG_DIR}") ## Change it later.
    create_directories_if_not_exist(f"{DATA_DIR}/logs")
    
    # FIG_DIR = local_fig_dir
    log("================================================")
    log("                   Parameters                   ")
    log("================================================")
    log(f"LOCUS:              {LOCUS}")
    log(f"VERSION:            {VERSION}")
    log(f"WORKDIR:            {DATA_DIR}")
    log(f"MAX_THREADS:        {MAX_THREADS}")
    log("================================================")
    
    # log(f"Maximum number of available threads: {nk.getMaxNumberOfThreads()}")
    log(f"Setting max number of threads to {MAX_THREADS}.")
    
    nk.setNumberOfThreads(MAX_THREADS)
    # load the mapping of junction_aa to node_id
    cdr3_to_id, id_to_cdr3 = load_junction_aa_mapping(MAP_FILE)
    #Load the epitope database
    epitope_df = load_epitope_info(EPITOPE_INFO_FILE)
    #load the graph
    g = load_graph(GRAPH_FILE)
    
    cdr3_to_id_keys_set = set(cdr3_to_id.keys())
    epitope_cdr3_set = set(epitope_df.junction_aa.values)
    
    start = time.time()
    common_cdr3 = cdr3_to_id_keys_set & epitope_cdr3_set
    log(f"Time to get common cdr3 : {time.time() - start:.2f} sec\n")
    # common_cdr3 = np.intersect1d(list(seq_dict.keys()), epitope_df.junction_aa.values)
    log(f"Total number of common CDR3: {len(common_cdr3)}")
    epitope_df = epitope_df[epitope_df.junction_aa.isin(common_cdr3)]
    epitope_df.to_csv('common_trb_and_epitopes_akc_v2.csv', index = False)
    log(f"Epitope DF Shape: {epitope_df.shape}")
    log("Plot revised epitope df figures.")
    epitope_df_figures(epitope_df)
    # Create bipartite mappings from epitope_df
    epitope_to_cdr3, cdr3_to_epitope = create_bipartite_mapping(epitope_df)

    # Get the number of unique receptors (including neighbors) per epitope
    epitope_receptor_count = get_receptor_count_per_epitope(g, epitope_to_cdr3, cdr3_to_id, id_to_cdr3)
    log(f"Epitope receptor counts length : {len(epitope_receptor_count)}")
    
    epitope_degree = epitope_df.groupby("akc_epitope_seq_aa")["junction_aa"].nunique()
    epitope_degree = epitope_degree.reset_index()
    
    epitope_degree.columns = ["akc_epitope_seq_aa", "real_receptors"]
    inf_df = pd.DataFrame(epitope_receptor_count.items(),columns=["akc_epitope_seq_aa", "inferred_receptors"])
    epitope_specificity_df = epitope_degree.merge(inf_df, on="akc_epitope_seq_aa")
    epitope_specificity_df["new_receptors"] = (epitope_specificity_df["inferred_receptors"] - epitope_specificity_df["real_receptors"])
    
    # Get the epitope distribution (number of unique epitopes per receptor)
    receptor_epitope_distribution = get_epitope_distribution_per_receptor(g, cdr3_to_epitope, cdr3_to_id, id_to_cdr3)
    log(f" Len of receptor_epitope_distribution: {len(receptor_epitope_distribution)}")
    
    cdr3_degree = epitope_df.groupby("junction_aa")["akc_epitope_seq_aa"].nunique()
    cdr3_degree = cdr3_degree.reset_index()
    cdr3_degree.columns = ["cdr3", "real_epitopes"]
    inf_receptor_df = pd.DataFrame(receptor_epitope_distribution.items(),columns=["cdr3", "inferred_epitopes"])
    receptor_specificity_df = cdr3_degree.merge(inf_receptor_df,on="cdr3",how="left")
    receptor_specificity_df["new_epitopes"] = ( receptor_specificity_df["inferred_epitopes"] - receptor_specificity_df["real_epitopes"])
    
    # epitope_specificity_df.to_csv("epitope_specificity_df.csv", index = False)
    # receptor_specificity_df.to_csv("receptor_specificity_df.csv", index = False)
    
    plot_inferred_specificity(epitope_specificity_df, receptor_specificity_df)
    
    
if __name__ == "__main__":
    main()
    