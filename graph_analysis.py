import time
import igraph as ig
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG_FILE = "graph_analysis_with_indels.log"
# GRAPH_FILE = "graph_full.pkl"  # iGraph pickle file
fig_dir = "figures_with_indels"

WORKDIR="graph_data/"
EDGE_FILE=f"{WORKDIR}/edges_id_with_indels.tsv"
MAP_FILE=f"{WORKDIR}/cdr3_to_id_with_indels.tsv"

# Create directory if it doesn't exist
os.makedirs(fig_dir, exist_ok=True)

def log(msg):
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(str(msg) + "\n")
# ================================================================
# Building Graph From Integer Edge list
# ================================================================
start = time.time()
graph_time = time.time()
id_dict = {}
try:
    with open(MAP_FILE, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 2:
                print("len: ", len(parts))
                print("parts: ", parts)
                continue
            else:
                id_str, seq = parts[0], parts[1]
                id_dict[int(id_str)] = seq
except Exception as e:
    print("Exception: ", e)

print("Size of id_dict (Node) : ", len(id_dict))

log(f"Map file reading time: {time.time() - start:.2f} sec\n")

start = time.time()
g = ig.Graph.Read_Edgelist(EDGE_FILE, directed=False)
log(f"Edge reading Time: {time.time() - start:.2f} sec\n")

start = time.time()
# Map node IDs to sequences
g.vs["cdr3"] = [id_dict.get(i, None) for i in range(g.vcount())]

log(f"Graph read: {g.vcount():,} nodes, {g.ecount():,} edges")
log(f"Edge Mapping Time: {time.time() - start:.2f} sec\n")
log(f"Total time to build the graph: {time.time() - graph_time:.2f} sec\n")

# ===========================================================================
# Step 3: Simplify the graph if needed
# ===========================================================================
log("Simplifying graph to remove duplicate edges and self-loops...")
start_simplify = time.time()
g.simplify()
log(f"Graph read after simplify: {g.vcount():,} nodes, {g.ecount():,} edges")
log(f"Time elapsed for simplify: {time.time() - start_simplify:.2f} seconds")


# ================================================================
# Load the graph
# ================================================================
# start = time.time()

# g = ig.Graph.Read_Pickle(GRAPH_FILE)
# stats_time = time.time() - start
# log(f"\nGraph loaded in {stats_time:.2f} s")

# is_directed = g.is_directed()
# log(f"Is the loaded graph directed? {is_directed}")
# if g.is_directed():
#     g = g.as_undirected() 
#     log(f"Graph converted to undirected for mutational analysis.")

# ================================================================
# Degree statistics
# ================================================================
start = time.time()
log("\nComputing graph statistics...")

degrees = g.degree()
log("Degree statistics:")
# log(f"Degrees: {degrees}")
log(f"Min degree: {min(degrees)}")
log(f"Max degree: {max(degrees)}")
log(f"Mean degree: {sum(degrees)/len(degrees):.2f}")
log(f"Finding degree stats time: {time.time() - start:.2f} sec\n")

# ================================================================
# Finding clustering coefficients
# ================================================================
start = time.time()
log(f"clustering coefficient: {g.transitivity_undirected()}")
log(f"Finding clustering coefficiet time: {time.time() - start:.2f} sec\n")

# ================================================================
# Plot the degree distribution as a histogram
# ================================================================
start = time.time()
# 'density=True' normalizes the histogram to represent probabilities
plt.figure(figsize=(12,10))
plt.hist(degrees, bins=range(min(degrees), max(degrees) + 2), density=True) 
plt.xlabel("Degree (k)")
plt.ylabel("Probability P(k)")
plt.title("Degree Distribution")
plt.savefig(f"{fig_dir}/igraph_degree_distribution.png", dpi = 300)
plt.close()

# plot degreee distributions using loglog scale
degrees = np.array(degrees)  # degree sequence
hist, bin_edges = np.histogram(degrees, bins=range(max(degrees)+1))
plt.figure(figsize=(12,10))
plt.loglog(bin_edges[:-1], hist, marker='o', linestyle='None')
plt.xlabel("Degree k")
plt.ylabel("Number of nodes")
plt.title("Degree Distribution (log-log)")
plt.grid(True)
plt.savefig(f"{fig_dir}/degree_distribution_loglog.png", dpi=300)
plt.close()

log(f"Degree distribution plot time: {time.time() - start:.2f} sec\n")


# ================================================================
# Find Largest connected components
# ================================================================
start = time.time()
log("Find number of connected components.")
# Connected components
components = g.connected_components()
num_components = len(components)
lcc_index = components.sizes().index(max(components.sizes()))
lcc_size = components.sizes()[lcc_index]

log(f"Number of connected components: {num_components:,}")
log(f"Largest connected component (nodes): {lcc_size:,}")

# Extract LCC
lcc_graph = g.subgraph(components[lcc_index])
log(f"LCC edges: {lcc_graph.ecount():,}")

log(f"Largest connected components finding time: {time.time() - start:.2f} sec\n")
# ================================================================
# Plot Component size distribution
# ================================================================

start = time.time()
component_sizes = components.sizes()
plt.figure(figsize=(12,10))
plt.hist(component_sizes, bins=100, log=True, color="blue", alpha=0.7)
plt.xlabel("Component size")
plt.ylabel("Frequency (log scale)")
plt.title("Connected component size distribution")
plt.savefig(f"{fig_dir}/component_size_distribution.png", dpi=300)
plt.close()

log(f"Component size plot time: {time.time() - start:.2f} sec\n")

# ================================================================
# Plot Component size distribution without the LCC
# ================================================================

start = time.time()
sizes = np.array(components.sizes())
largest = sizes.max()
sizes_no_lcc = sizes[sizes != largest]

plt.figure(figsize=(12,10))
plt.hist(sizes_no_lcc, bins=100, log=True, color="blue", alpha=0.7)
plt.xscale("log")      # optional but recommended
plt.xlabel("Component size (log scale)")
plt.ylabel("Frequency (log scale)")
plt.title("Connected component size distribution (excluding LCC)")
plt.savefig(f"{fig_dir}/component_size_distribution_no_LCC.png", dpi=300)
plt.close()

log(f"No LLC component size plot time: {time.time() - start:.2f} sec\n")

# ================================================================
# Plot Degree vs local clustering coefficient
# ================================================================

start = time.time()
lcc_degrees = lcc_graph.degree()
plt.figure(figsize=(12,10))
clustering = lcc_graph.transitivity_local_undirected(mode="zero")
plt.scatter(lcc_degrees, clustering, alpha=0.2, s=2)
plt.xscale("log")
plt.xlabel("Degree")
plt.ylabel("Local clustering coefficient")
plt.title("Degree vs Clustering")
plt.savefig(f"{fig_dir}/degree_vs_clustering.png", dpi=300)
plt.close()

log(f"Degree vs clustering coef plot time: {time.time() - start:.2f} sec\n")

# ================================================================
# Plot a hub node
# ================================================================

# Select hub nodes with a desired threshold
start = time.time()

degree_threshold = 150
hub_nodes = [v.index for v in lcc_graph.vs if lcc_graph.degree(v) > degree_threshold]
log(f"Number of hubs (degree > {degree_threshold}): {len(hub_nodes):,}")
degree_threshold = 200
hub_nodes = [v.index for v in lcc_graph.vs if lcc_graph.degree(v) > degree_threshold]
log(f"Number of hubs (degree > {degree_threshold}): {len(hub_nodes):,}")

degree_threshold = 250
hub_nodes = [v.index for v in lcc_graph.vs if lcc_graph.degree(v) > degree_threshold]
log(f"Number of hubs (degree > {degree_threshold}): {len(hub_nodes):,}")

degree_threshold = 270
hub_nodes = [v.index for v in lcc_graph.vs if lcc_graph.degree(v) > degree_threshold]
log(f"Number of hubs (degree > {degree_threshold}): {len(hub_nodes):,}")

#extract subgraph of hubs
hub_subgraph = lcc_graph.subgraph(hub_nodes)
print(f"Hub subgraph nodes: {hub_subgraph.vcount():,}, edges: {hub_subgraph.ecount():,}")


layout = hub_subgraph.layout("large_graph")
visual_style = {}
visual_style["vertex_size"] = 20
visual_style["vertex_color"] = "gray"
visual_style["vertex_label"] = hub_subgraph.vs["cdr3"]
visual_style["vertex_label_color"] = "black"
visual_style["vertex_label_angle"] = 90
visual_style["layout"] = layout
visual_style["bbox"] = (1600, 1600)
visual_style["margin"] = 40
visual_style["shape"] = "diamond"

fig, ax = plt.subplots()
figname = f"{fig_dir}/lcc_hub_subgraph.png"
ig.plot(hub_subgraph, figname, **visual_style)
plt.close()
log(f"LCC hub node plot time: {time.time() - start:.2f} sec\n")

# # # Community detection using Louvain / multilevel
# # communities = lcc_graph.community_multilevel()
# # community_sizes = [len(c) for c in communities]

# # print(f"LCC graph number of communities: {len(communities):,}")
# # print(f"LCC graph Largest community size: {max(community_sizes):,}")

# # # Community size distribution
# # plt.figure(figsize=(7,5))
# # plt.hist(community_sizes, bins=50, log=True, color="green", alpha=0.7)
# # plt.xlabel("Community size")
# # plt.ylabel("Frequency (log scale)")
# # plt.title("Community size distribution")
# # plt.savefig(f"{fig_dir}/community_size_distribution.png", dpi=300)
# # plt.close()


# ================================================================
# Calculate diameter, betweenness and 
# assortativity for largest connected components
# ================================================================
start = time.time()
#average shortest path for LCC graph
# log(f"LCC avgerage shortest path: {lcc_graph.average_path_length()}")
#The longest shortest path in the graph for largest connected components
log(f"LCC Diameter: {lcc_graph.get_diameter()}")
log(f"LCC Diameter calculation time: {time.time() - start:.2f} sec\n")

start = time.time()
# For LCC only (computationally heavy)
bc = lcc_graph.betweenness()
log(f"Betweenness: min={min(bc)}, max={max(bc)}, mean={sum(bc)/len(bc):.2f}")
log(f"LCC betweenness calculation time: {time.time() - start:.2f} sec\n")

#Check whether high-degree sequences tend to connect to other high-degree sequences.
start = time.time()
assortativity = lcc_graph.assortativity_degree(directed=False)
log(f"LCC Graph Degree assortativity: {assortativity:.3f}")
log(f"LCC Graph Degree assortativity calculation time: {time.time() - start:.2f} sec\n")

# stats_time = time.time() - start
# log(f"\n[igraph] Stats computed in {stats_time:.2f} s")

# ================================================================
# End
# ================================================================