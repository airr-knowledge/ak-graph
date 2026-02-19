import networkit as nk
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import powerlaw
from networkit import vizbridges

# ==========================================================================================
# File paths
# ==========================================================================================
LOG_FILE = "graph_analysis_with_networkit_v3.log"
fig_dir = "figures_with_networkit_v2"
edge_file = "output_pairs_trb_v3.tsv"
map_file = "output_seq_map_trb_v3.tsv"

# WORKDIR = "graph_data/"
# EDGE_FILE = f"{WORKDIR}/edges_sorted.tsv"
# MAP_FILE = f"{WORKDIR}/cdr3_to_id_with_indels.tsv"


WORKDIR = "graph_data/"
EDGE_FILE = f"{edge_file}"
MAP_FILE = f"{map_file}"

os.makedirs(fig_dir, exist_ok=True)

# ==========================================================================================
# Logging function
# ==========================================================================================
def log(msg):
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(str(msg) + "\n")

# log(f"Maximum number of available threads: {nk.getMaxNumberOfThreads()}")
log(f"Setting max number of threads to 64.")
nk.setNumberOfThreads(64)

log(f"Reading using the binary file.")
# ==============================
# Read node mapping (cdr3 sequences)
# ==============================
start = time.time()
id_dict = {}
with open(MAP_FILE, "r") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            id_dict[int(parts[0])] = parts[1]

log(f"Loaded {len(id_dict):,} nodes")
log(f"Mapping read time: {time.time() - start:.2f} sec\n")

# ==============================
# Load graph from edge list
# ==============================
start = time.time()
g = nk.readGraph(
    EDGE_FILE,
    nk.Format.EdgeList,
    separator="\t",
    firstNode=0,
    continuous=True,  # ensures all nodes 0..max_id are included, including isolates
    directed=False
)
log(f"Graph loaded with {g.numberOfNodes():,} nodes and {g.numberOfEdges():,} edges")
log(f"Edge load time: {time.time() - start:.2f} sec\n")

# ==============================
# Save the network in a binary format
# ==============================
start = time.time()

nk.writeGraph(g, f"{WORKDIR}/trb_graph_with_indel_v3.nkbg003", nk.Format.NetworkitBinary, chunks=16)
log(f"Time to write the file in binary: {time.time() - start:.2f} sec\n")


# # ==========================================================================================
# # Load graph from binary format
# # ==========================================================================================
# start = time.time()
# g = nk.readGraph(f"{WORKDIR}/trb_graph_with_indel_v1.nkbg003", nk.Format.NetworkitBinary)
# log(f"Time to read the binary graph file: {time.time() - start:.2f} sec\n")

# ==========================================================================================
# Network overview
# ==========================================================================================
start = time.time()
log("Graph overview with all nodes.")
log(f"{nk.overview(g)}")

log(f"Graph overview time: {time.time() - start:.2f} sec\n")

# ==========================================================================================
# Degree distribution Before removing isolates
# ==========================================================================================
start = time.time()
dd = sorted(nk.centrality.DegreeCentrality(g).run().scores(), reverse=True)
degrees, numberOfNodes = np.unique(dd, return_counts=True)
plt.figure(figsize=(12,10))
plt.xscale("log")
plt.xlabel("degree")
plt.yscale("log")
plt.ylabel("number of nodes")
plt.plot(degrees, numberOfNodes)
plt.title("Degree distribution before removing isolates")
plt.savefig(f"{fig_dir}/graph_degree_distribution_with_isolates.png", dpi = 300)
plt.close()
log(f"Time for plotting degree distribution: {time.time() - start:.2f} sec\n")


# ==========================================================================================
# Removing isolates
# ==========================================================================================
start = time.time()
#Remove isolates and calculate largest connected components
isolates = [u for u in g.iterNodes() if g.degree(u) == 0]

# Remove isolates
for u in isolates:
    g.removeNode(u)

print(f"Removed {len(isolates)} isolated nodes. New node count: {g.numberOfNodes()}")
log(f"Time to remove isolates from the graph: {time.time() - start:.2f} sec\n")

# ==========================================================================================
# Network overview after removing isolates
# ==========================================================================================
start = time.time()
log(f"{nk.overview(g)}")

log(f"Graph overview time after removing isolates: {time.time() - start:.2f} sec\n")

# ==========================================================================================
# Degree distribution after removing isolates
# ==========================================================================================
start = time.time()
dd = sorted(nk.centrality.DegreeCentrality(g).run().scores(), reverse=True)
degrees, numberOfNodes = np.unique(dd, return_counts=True)
plt.figure(figsize=(12,10))
plt.xscale("log")
plt.xlabel("degree")
plt.yscale("log")
plt.ylabel("number of nodes")
plt.plot(degrees, numberOfNodes)
plt.savefig(f"{fig_dir}/graph_degree_distribution_without_isolates.png", dpi = 300)
plt.close()
log(f"Time for plotting degree distribution: {time.time() - start:.2f} sec\n")


# ==========================================================================================
# ClusteringPerDegree
# ==========================================================================================
start = time.time()
plt.figure(figsize=(16,12))
nk.plot.clusteringPerDegree(g)
plt.tight_layout()
plt.savefig(f"{fig_dir}/clustering_per_degreee_without_isolates.png", dpi=300)
plt.close()
log(f"ClusteringPerDegree plot time: {time.time() - start:.2f} sec\n")
# ==========================================================================================
# Extract Largest Connected Component and write it in a file
# ==========================================================================================
start = time.time()
cc = nk.components.ConnectedComponents(g)
cc.run()
largest_component = cc.extractLargestConnectedComponent(g, compactGraph=False)
log(f"Largest connected component finiding time: {time.time() - start:.2f} sec\n")

start = time.time()
nk.writeGraph(largest_component, f'{WORKDIR}/trb_largest_connected_component_v1.edgelist.txt', nk.Format.EdgeListTabOne)
log(f"Largest connected component write time: {time.time() - start:.2f} sec\n")


# ==========================================================================================
# Connected Component Sizes plot without isolates and biggest one at the end
# ==========================================================================================
start = time.time()
cc = nk.components.ConnectedComponents(g)
cc.run()
component_sizes = cc.getComponentSizes().values()        # list of component sizes
num_components = len(component_sizes)
print(f"Found {num_components} connected components")

sizes = list(component_sizes)
giant = max(sizes)
rest = [s for s in sizes if s != giant]
plt.figure(figsize=(10,7))
plt.hist(rest, bins=np.logspace(np.log10(1), np.log10(max(rest)), 100), log=True, color='blue', alpha=0.7)
plt.xscale('log')
plt.axvline(giant, color='red', linestyle='--', label=f"Giant component = {giant}")
plt.xlabel("Component size")
plt.ylabel("Count (log)")
plt.title("Component size distribution (excluding giant component)")
plt.legend()
plt.savefig(f"{fig_dir}/component_size_distribution_without_isolates.png", dpi=300)
plt.close()
log(f"Connected component sizes plot time: {time.time() - start:.2f} sec\n")


# ==========================================================================================
# Powerlaw fit after removing isolates
# ==========================================================================================
start = time.time()
deg = nk.centrality.DegreeCentrality(g).run().scores()
deg_nonzero = [d for d in deg if d > 0]

fit = powerlaw.Fit(deg_nonzero, discrete=True)
log("Distribution comparisons (normalized ratio):")
for dist in ['exponential', 'lognormal','stretched_exponential', 'truncated_power_law']:
    R, p = fit.distribution_compare('power_law', dist, normalized_ratio=True)
    log(f"Power-law vs {dist}: R={R:.3f}, p={p:.3f}")


fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# ----------------------------
# PDF plot (left)
# ----------------------------
ax = axes[0]
fit.plot_pdf(color='b', linewidth=2, ax=ax, label='Empirical PDF')
fit.power_law.plot_pdf(color='r', linestyle='--', ax=ax, label='Power-law')
fit.exponential.plot_pdf(color='g', linestyle=':', ax=ax, label='Exponential')
fit.lognormal.plot_pdf(color='m', linestyle='-.', ax=ax, label='Lognormal')
fit.stretched_exponential.plot_pdf(color='orange', linestyle='-', ax=ax, label='Stretched exp')

# Mark xmin for power-law
ax.axvline(fit.power_law.xmin, color='r', linestyle=':', label='Power-law xmin')
ax.set_xlabel("Degree")
ax.set_ylabel("PDF")
ax.set_title("Degree Distribution (PDF)")
ax.legend()

# ----------------------------
# CCDF plot (right)d
# ----------------------------
ax = axes[1]
fit.plot_ccdf(color='b', linewidth=2, ax=ax, label='Empirical CCDF')
fit.power_law.plot_ccdf(color='r', linestyle='--', ax=ax, label='Power-law')
fit.exponential.plot_ccdf(color='g', linestyle=':', ax=ax, label='Exponential')
fit.lognormal.plot_ccdf(color='m', linestyle='-.', ax=ax, label='Lognormal')
fit.stretched_exponential.plot_ccdf(color='orange', linestyle='-', ax=ax, label='Stretched exp')

# Mark xmin for power-law
ax.axvline(fit.power_law.xmin, color='r', linestyle=':', label='Power-law xmin')
ax.set_xlabel("Degree")
ax.set_ylabel("CCDF")
ax.set_title("Degree Distribution (CCDF)")
ax.legend()

plt.tight_layout()
plt.savefig(f"{fig_dir}/degree_distribution_pdf_ccdf.png", dpi=300)
plt.show()

log(f"Time for powerlaw fit: {time.time() - start:.2f} sec\n")

# ==========================================================================================
# Community Dectection
# ==========================================================================================

start = time.time()
communities = nk.community.detectCommunities(g)

nk.community.Modularity().getQuality(communities, g)
log(f"Time Community detection: {time.time() - start:.2f} sec\n")

# ==============================
# Core Decomposition
# A k-core decomposition of a graph is performed by successicely peeling away nodes with degree less than. The remaining nodes form the k-core of the graph.
# ==============================
start = time.time()
coreDec = nk.centrality.CoreDecomposition(g)
coreDec.run()

cores = coreDec.scores()
plt.figure(figsize=(12,10))
cores = np.array(coreDec.scores())
plt.hist(cores, bins=200)
plt.xlabel("Core Number")
plt.ylabel("# Nodes")
plt.title("k-core Distribution")
plt.savefig(f"{fig_dir}/kcore_distribution.png", dpi=300)
plt.close()
log(f"Time Core Decomposition plot: {time.time() - start:.2f} sec\n")

# # ==========================================================================================
# # Betweenness Calculation
# # ==========================================================================================

# ab = nk.centrality.ApproxBetweenness(g, epsilon=0.001)
# ab.run()
# scores = ab.scores()
# log(f"Approximate betweenness time: {time.time() - start:.2f} sec\n")

# log(f"Approximate betweenness statistics")
# scores = np.array(scores)
# log("Betweenness stats:")
# log(f"Min: {scores.min():.4f}")
# log(f"Max: {scores.max():.4f}")
# log(f"Mean: {scores.mean():.4f}")
# log(f"Median: {np.median(scores):.4f}")
# log(f"Std: {scores.std():.4f}")


# plt.figure(figsize=(12,10))
# plt.hist(scores, bins=100, log=True, color='skyblue', edgecolor='black')
# plt.xlabel("Betweenness centrality")
# plt.ylabel("Number of nodes")
# plt.title("Betweenness Centrality Distribution")
# plt.tight_layout()
# plt.savefig(f"{fig_dir}/betweenness_centrality_distribution.png", dpi=300)
# plt.show()


# log(f"Core decomposition and plot total time: {time.time() - start:.2f} sec\n")
# # ==========================================================================================
# # END
# # ==========================================================================================