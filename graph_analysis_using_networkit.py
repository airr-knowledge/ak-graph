import networkit as nk
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import powerlaw
from networkit import vizbridges
import argparse

def log(msg):
    '''
    Logging function. Log all the print outputs
    '''
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(str(msg) + "\n")
    
def build_and_save_graph(EDGE_FILE, graph_output_path):
    '''
    Load graph from edge list and save the graph in networkit graph format for faster load later
    '''
    start = time.time()
    g = nk.readGraph(
        EDGE_FILE,
        nk.Format.EdgeList,
        separator="\t",
        firstNode=0,
        continuous=True,  # ensures all nodes 0...max_id are included, including isolates
        directed=False
    )
    log(f"Graph loaded with {g.numberOfNodes():,} nodes and {g.numberOfEdges():,} edges")
    log(f"Edge load time: {time.time() - start:.2f} sec\n")
    
    start = time.time()
    nk.writeGraph(g, graph_output_path, nk.Format.NetworkitBinary, chunks=16)
    log(f"Time to write the file in binary: {time.time() - start:.2f} sec\n")
    
    return g

def load_junction_aa_mapping(MAP_FILE):
    '''
    Read node mapping (junction_aa sequences)
    '''
    log(f"Reading the map file.")
    start = time.time()
    id_dict = {}
    with open(MAP_FILE, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                id_dict[int(parts[0])] = parts[1]

    log(f"Loaded {len(id_dict):,} nodes")
    log(f"Mapping read time: {time.time() - start:.2f} sec\n")
    return id_dict


def load_graph(graph_path):
    '''
    Load graph from binary format
    '''
    start = time.time()
    g = nk.readGraph(graph_path, nk.Format.NetworkitBinary)
    log(f"Time to read the binary graph file: {time.time() - start:.2f} sec\n")
    return g

def comunity_detection_stats(g, FIG_DIR):
    '''
    Community Dectection statistics
    '''
    start = time.time()
    communities = nk.community.detectCommunities(g)
    
    nk.community.Modularity().getQuality(communities, g)
    log(f"Time Community detection: {time.time() - start:.2f} sec\n")
    

def core_decomposition_stats(g, FIG_DIR):
    '''
    Core Decomposition
    A k-core decomposition of a graph is performed by successicely peeling away nodes with degree
    less than the remaining nodes form the k-core of the graph.
    '''
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
    plt.savefig(f"{FIG_DIR}/kcore_distribution.png", dpi=300)
    plt.close()
    log(f"Time Core Decomposition plot: {time.time() - start:.2f} sec\n")
    
def analyze_power_law(g, FIG_DIR):
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
    # PDF plot (left)

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

    # CCDF plot (right)d
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
    plt.savefig(f"{FIG_DIR}/degree_distribution_pdf_ccdf.png", dpi=300)
    plt.show()

    log(f"Time for powerlaw fit: {time.time() - start:.2f} sec\n")


def analyze_betweenness(g, FIG_DIR):
    '''
    Betweenness Calculation
    '''
    ab = nk.centrality.ApproxBetweenness(g, epsilon=0.001)
    ab.run()
    scores = ab.scores()
    log(f"Approximate betweenness time: {time.time() - start:.2f} sec\n")

    log(f"Approximate betweenness statistics")
    scores = np.array(scores)
    log("Betweenness stats:")
    log(f"Min: {scores.min():.4f}")
    log(f"Max: {scores.max():.4f}")
    log(f"Mean: {scores.mean():.4f}")
    log(f"Median: {np.median(scores):.4f}")
    log(f"Std: {scores.std():.4f}")


    plt.figure(figsize=(12,10))
    plt.hist(scores, bins=100, log=True, color='skyblue', edgecolor='black')
    plt.xlabel("Betweenness centrality")
    plt.ylabel("Number of nodes")
    plt.title("Betweenness Centrality Distribution")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/betweenness_centrality_distribution.png", dpi=300)
    plt.show()


    log(f"Core decomposition and plot total time: {time.time() - start:.2f} sec\n")


def get_general_graph_statistics(g, FIG_DIR):
    print("We are here.")
    # Network overview
    start = time.time()
    log("Graph overview with all nodes.")
    log(f"{nk.overview(g)}")

    log(f"Graph overview time: {time.time() - start:.2f} sec\n")

    
    # Degree distribution Before removing isolates
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
    plt.savefig(f"{FIG_DIR}/graph_degree_distribution_with_isolates.png", dpi = 300)
    plt.close()
    log(f"Time for plotting degree distribution: {time.time() - start:.2f} sec\n")
    
    # Removing isolates
    start = time.time()
    #Remove isolates and calculate largest connected components
    isolates = [u for u in g.iterNodes() if g.degree(u) == 0]
    # Remove isolates
    for u in isolates:
        g.removeNode(u)
    print(f"Removed {len(isolates)} isolated nodes. New node count: {g.numberOfNodes()}")
    log(f"Time to remove isolates from the graph: {time.time() - start:.2f} sec\n")


    # Network overview after removing isolates
  
    start = time.time()
    log(f"{nk.overview(g)}")
    log(f"Graph overview time after removing isolates: {time.time() - start:.2f} sec\n")

    # Degree distribution after removing isolates
    start = time.time()
    dd = sorted(nk.centrality.DegreeCentrality(g).run().scores(), reverse=True)
    degrees, numberOfNodes = np.unique(dd, return_counts=True)
    plt.figure(figsize=(12,10))
    plt.xscale("log")
    plt.xlabel("degree")
    plt.yscale("log")
    plt.ylabel("number of nodes")
    plt.plot(degrees, numberOfNodes)
    plt.savefig(f"{FIG_DIR}/graph_degree_distribution_without_isolates.png", dpi = 300)
    plt.close()
    log(f"Time for plotting degree distribution: {time.time() - start:.2f} sec\n")

    # ClusteringPerDegree
    start = time.time()
    plt.figure(figsize=(16,12))
    nk.plot.clusteringPerDegree(g)
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/clustering_per_degreee_without_isolates.png", dpi=300)
    plt.close()
    log(f"ClusteringPerDegree plot time: {time.time() - start:.2f} sec\n")
    
    # Extract Largest Connected Component and write it in a file
    start = time.time()
    cc = nk.components.ConnectedComponents(g)
    cc.run()
    largest_component = cc.extractLargestConnectedComponent(g, compactGraph=False)
    log(f"Largest connected component finiding time: {time.time() - start:.2f} sec\n")

    start = time.time()
    nk.writeGraph(largest_component, f'{WORKDIR}/trb_largest_connected_component_v1.edgelist.txt', nk.Format.EdgeListTabOne)
    log(f"Largest connected component write time: {time.time() - start:.2f} sec\n")


    # Connected Component Sizes plot without isolates and biggest one at the end
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
    plt.savefig(f"{FIG_DIR}/component_size_distribution_without_isolates.png", dpi=300)
    plt.close()
    log(f"Connected component sizes plot time: {time.time() - start:.2f} sec\n")


WORKDIR = "./ak_graph_data/"

LOG_FILE = None
def main():
    parser = argparse.ArgumentParser("Generate graph using networkit.")
    parser.add_argument("analysis_type", help="Either build_graph or graph_stat")
    parser.add_argument("locus", help="Locus to build/analyze the graph (tra, trb, trd, trg)")
    parser.add_argument("--version", default="v1", help="Version of the table name that will be put on the graph")
    parser.add_argument("--max_threads", default=64, help="Maximum number of threads to use by the graph algorithm")
    
    args = parser.parse_args()
    analysis_type = args.analysis_type
    locus = args.locus
    version = args.version
    max_threads = args.max_threads
    
    print("Parameters: ")
    print("analysis_type: ", analysis_type)
    print("locus: ", locus)
    print("version: ", version)
    print("max_threads: ", max_threads)
    global LOG_FILE
    EDGE_FILE = f"{WORKDIR}/{locus}_output_pairs_{version}.tsv"
    GRAPH_OUTPUT_PATH = f"{WORKDIR}/{locus}_graph_{version}.nkbg003"
    MAP_FILE = f"{WORKDIR}/{locus}_output_seq_map_{version}.tsv"
    FIG_DIR = f"{WORKDIR}/figures/{locus}_{version}/"
    LOG_FILE = f"{WORKDIR}/logs/{locus}_analysis_{version}.log"
    #create figure directory if not exist
    os.makedirs(FIG_DIR, exist_ok=True)
    
    # log(f"Maximum number of available threads: {nk.getMaxNumberOfThreads()}")
    log(f"Setting max number of threads to 64.")
    
    nk.setNumberOfThreads(64)

    if analysis_type == 'build_graph':
        build_and_save_graph(EDGE_FILE, GRAPH_OUTPUT_PATH)
    elif analysis_type == 'graph_stat':
        g = load_graph(GRAPH_OUTPUT_PATH)
        get_general_graph_statistics(g, FIG_DIR)
    else:
        print("Wrong analysis type. Please put build_graph or graph_stat")

if __name__ == "__main__":
    main()
    
        