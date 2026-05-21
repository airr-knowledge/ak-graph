import pandas as pd
import numpy as np

import sys
import re
import argparse
from dotenv import dotenv_values
import os
import time

import json
import sqlite3
import csv
import gc

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import seaborn as sns
import logomaker
from matplotlib.ticker import MaxNLocator

from collections import defaultdict
from collections import Counter

# SQLITE_DB_PATH = "/ak_graph_data/airrkb_v2_tilde.db"
SQLITE_DB_PATH = "./airrkb_v2_tilde.db"
# GRAPH_FILE = f"{DATA_DIR}/graph_files/{LOCUS}_graph_{VERSION}.nkbg003"
# MAP_FILE = f"{DATA_DIR}/pair_files/{LOCUS}_output_seq_map_{VERSION}.tsv"

# --------------------------------------------------------------------------------- 
# Load repertoire sequences and find locus and species
# --------------------------------------------------------------------------------- 
def load_airr_file(filepath):
    """Load AIRR TSV and extract junction_aa."""
    df = pd.read_csv(filepath, sep="\t", low_memory=False)
    df = df[(df.productive == 'T') | (df.productive == True)]
    if "junction_aa" not in df.columns:
        raise ValueError("AIRR file must contain 'junction_aa' column")
    
    df = df[['sequence_id', 'junction_aa', 'duplicate_count', 'locus']]
    print(df.head())

    return df

def get_species(rep_id, airr_metadata):
    with open(airr_metadata, 'r') as file:
        rep_data = json.load(file)

    for rep in rep_data['Repertoire']:
        repertoire_id = rep.get("repertoire_id", None)
        if repertoire_id == rep_id:
            species = rep.get("subject", {}).get('species', {}).get('id', None)
            return species
    return None


# ---------------------------------------------------------------------------------   
# Query Builders
# ---------------------------------------------------------------------------------         

def get_query_for_locus_prejoined(locus, species, chunk_size):
    placeholders = ', '.join(['?'] * chunk_size)
    species_renamed = species.replace(':', '_')
    table_name = f"{locus}_{species_renamed}_tilde".lower()
    return f"""
    SELECT
        akc_assay_akc_id,
        akc_complex_akc_id,
        akc_epitope_akc_id,
        akc_epitope_seq_aa,
        akc_source_protein,
        akc_source_organism,
        junction_aa,
        akc_species,
        akc_v_call,
        akc_j_call
    FROM {table_name} jtd
    WHERE jtd.junction_aa IN ({placeholders})
    """
    
def get_query_for_chain_prejoined(locus, species, chunk_size):
    placeholders = ",".join(["?"] * chunk_size)
    species_renamed = species.replace(':', '_')
    table_name = f"{locus}_{species_renamed}_tilde".lower()
    query = f"""
    SELECT junction_aa FROM {table_name} jtd WHERE jtd.junction_aa IN ({placeholders})
    """
    return query

def get_query_for_assay_object(chunk_size):
    placeholders = ",".join(["?"] * chunk_size)
    query = f"""SELECT * FROM "QueryAssay" qa WHERE qa.akc_id IN ({placeholders})"""
    return query

def get_connection(SQLITE_DB_PATH):
    
    db_path = f"file:{SQLITE_DB_PATH}?mode=ro&immutable=1"
    conn = sqlite3.connect(db_path, uri=True)
    # TACC optimizations
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode = OFF;")
    conn.execute("PRAGMA synchronous = OFF;")
    # conn.execute("PRAGMA temp_store = MEMORY;")
    # conn.execute("PRAGMA cache_size = -2000000;")  # ~2GB cache
    
    conn.row_factory = sqlite3.Row
    return conn


def chunk_list(data, chunk_size):
    """Helper function to chunk the data into smaller chunks."""
    for i in range(0, len(data), chunk_size):
        yield data[i:i + chunk_size]


def query_database_stream(parameter, query_type, SQLITE_DB_PATH, locus, species, chunk_size=1):
    # Open a fresh connection for this specific chunk
    for chunk in chunk_list(parameter, chunk_size):
        conn = get_connection(SQLITE_DB_PATH) 
        cur = conn.cursor()
        try:
            if query_type == "junction_aa":
                query = get_query_for_locus_prejoined(locus, species, len(chunk))
                cur.execute(query, tuple(chunk))
                # Yield rows one by one to keep memory low
                for row in cur:
                    yield row
            
            elif query_type == "assay":
                query = get_query_for_assay_object(len(chunk))
                cur.execute(query, tuple(chunk))
                for row in cur:
                    yield row
                    
            elif query_type == "chain":
                query = get_query_for_chain_prejoined(locus, species, len(chunk))
                cur.execute(query, tuple(chunk))
                for row in cur:
                    yield row
        finally:
            conn.close()


def process_query_results_stream(rows_generator):
    """
    Processes assay rows from a generator to keep memory usage low.
    Returns a DataFrame of metadata and the dictionary of assay objects.
    """
    processed_data = []
    all_assay_dict = {}

    # Iterate through the generator (one row at a time)
    for row in rows_generator:
        # Postgres might return 'assay_object' as a dict or a JSON string
        assay_raw = row["assay_object"]

        if not assay_raw:
            continue
        # Handle JSON parsing only if necessary
        if isinstance(assay_raw, str):
            assay_dict = json.loads(assay_raw)
        elif isinstance(assay_raw, dict):
            assay_dict = assay_raw
        else:
            continue

        assay_id = row["akc_id"]
        # Store the original object for the final JSON output
        all_assay_dict[assay_id] = assay_dict
        
        # Extract nested metadata
        specimen = assay_dict.get('specimen', {})
        participant = assay_dict.get('participant', {})
        investigation = assay_dict.get('investigation', {})

        # Build the flat metadata list
        processed_data.append({
            'akc_id': assay_id,
            'data_type': assay_dict.get('type'),
            'assay_type': assay_dict.get('assay_type'),
            'specimen_tissue': specimen.get('tissue'),
            'participant_species': participant.get('species'),
            'investigation_name': investigation.get('name'),
            'investigation_description': investigation.get('description')
        })

    # Convert the collected metadata to a DataFrame
    assay_df = pd.DataFrame(processed_data)
    # Instead of returning a massive string here, we return the dict.
    return assay_df, all_assay_dict

def plot_cdr3_vs_epitope_stats(summary_df, output_file_base, output_dir = './', n = 5):
    # plot top n junction_aa vs number of epitopes
    sns.set_theme()
    temp = summary_df.head(15)
    fig, axes = plt.subplots(1, 1, figsize = (7, 6))
    sns.barplot(data = temp, y = 'query_cdr3', x = 'n_unique_epitope_seq', ax = axes)
    axes.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("Number of Unique Epitope Sequence")
    plt.ylabel("CDR3")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{output_file_base}.tilde.top_n_cdr3_vs_epiope_distribution_figure.png", bbox_inches = 'tight', dpi=300)
    plt.close()
    
    # Calculate stats for cross reactivity plot
    cdr3_to_epitope_counts = defaultdict(lambda: defaultdict(int))
    
    for idx, row in summary_df.iterrows():
        cdr3 = row['query_cdr3']
        for ep in row['unique_epitope_seq'].split(','):
            cdr3_to_epitope_counts[cdr3][ep] = 1
    
    #Calculation for per epitope logo
    epitope_to_cdr3s = defaultdict(list)

    for _, row in summary_df.iterrows():
        for ep in row['unique_epitope_seq'].split(','):
            epitope_to_cdr3s[ep].append(row['query_cdr3'])
    # print(epitope_to_cdr3s)
            
    # Cross-reactivity / Sparsity Plots
    x = summary_df['n_unique_epitope_seq']

    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(x, bins=range(0, x.max() + 2), color='skyblue', edgecolor='black', align='left')  # bins are integers
    plt.xlabel("Number of unique epitopes per CDR3", fontsize=12)
    plt.yscale('log')
    plt.ylabel("Count of CDR3s(log)", fontsize=12)
    plt.title("CDR3 Cross-reactivity Distribution", fontsize=14)

    # Make x-axis integer only
    plt.xticks(range(0, x.max() + 1, max(1, x.max() // 10)))  # step intelligently based on range
    plt.yticks(fontsize=10)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{output_file_base}.tilde.cross_reactivity_distribution_plot.png", bbox_inches = 'tight', dpi=300)
    plt.close()
        
    # ----------------------------------------------------------------------------------------
    # Sequence logos for top 5 epitopes CDR3
    # ----------------------------------------------------------------------------------------
    # top_epitopes from your calculation
    top_epitopes = sorted(epitope_to_cdr3s.items(), key=lambda x: len(x[1]), reverse=True)[:n]
    top_epitopes = [ep for ep, cdr3s in top_epitopes]  # just the epitope strings
    
    for rank, ep in enumerate(top_epitopes, start=1):
        cdr3_list = epitope_to_cdr3s[ep]
        if len(cdr3_list) == 0:
            continue
        # Count lengths
        lengths_count = Counter([len(s) for s in cdr3_list])
        # Take top 3 most common lengths
        top_lengths = [l for l, _ in lengths_count.most_common(3)]
        plt.figure(figsize=(max(12, max(top_lengths)), 4 * len(top_lengths)))

        for i, l in enumerate(top_lengths):
            sequences = [s for s in cdr3_list if len(s) == l]
            if not sequences:
                continue

            counts_df = logomaker.alignment_to_matrix(sequences, to_type='counts')
            prob_df = counts_df.div(counts_df.sum(axis=1), axis=0)
            prob_df = prob_df.loc[:, (prob_df.sum(axis=0) > 0)]
            
            ax = plt.subplot(len(top_lengths), 1, i + 1)
            logomaker.Logo(
                prob_df,
                ax=ax,
                color_scheme='chemistry',
                shade_below=0.5,
                fade_below=0.5,
                vpad=0.05
            )
            ax.set_title(f"Length {l} | n={len(sequences)}")
            ax.set_xlabel("Position in CDR3")
            ax.set_ylabel("Probability")

        
        # save_path = f"{output_dir}/{output_file_base}_logo_{ep.replace('/', '_')}_multi_length.png"
        save_path = f"{output_dir}/{output_file_base}.tilde.top_epitope_rank_{rank}.png"
        plt.suptitle(f"Top Epitope {rank}: {ep}")
        plt.tight_layout()
        plt.savefig(save_path, bbox_inches = 'tight', dpi=300)
        plt.close()
        print(f"Saved multi-length logo for epitope {ep} -> {save_path}")         


def create_directories_if_not_exist(path):
    """Create directories if they do not exist"""
    if not os.path.exists(path):
        os.makedirs(path)

def main(SQLITE_DB_PATH): 
    parser = argparse.ArgumentParser("Please provide the parameters for tilde analysis ")
    parser.add_argument("--data_dir", default="./local_ak_graph_data_v2/", help="Version of the table name that will be put on the graph")
    parser.add_argument("--input_file", default= "test_input_blood.airr.tsv", help="Name of the input file")
    parser.add_argument("--output_file_base", default="test_input_blood", help="Output file base name")
    # parser.add_argument("--input_file", default= "test_input.airr.tsv", help="Name of the input file")
    # parser.add_argument("--output_file_base", default="test_input", help="Output file base name")
    parser.add_argument("--top_n_epitopes", default=5, type=int, help="Top N epitopes to plot")
    parser.add_argument("--AIRRMetadata",default='repertoires.airr.json', help="Airr Metadata File")
    
    args = parser.parse_args()
    input_file = args.input_file
    output_file_base = args.output_file_base
    top_n_epitopes = args.top_n_epitopes
    airr_metadata = args.AIRRMetadata
    data_dir = args.data_dir

    output_dir = f"{data_dir}/sqlite/tilda_output"
    
    #create figure directory if not exist
    create_directories_if_not_exist(output_dir)
    
    #read airr file for locus information and junction_aa list
    airr_df = load_airr_file(input_file)
    locus = airr_df.locus.unique()
    print("Unique Locus and Size: ", len(locus))
    locus = locus[0].lower()
    #repertoire id should be output_file_base. If not then change here.
    rep_id = input_file.strip().split('.')[0]
    species = get_species(rep_id, airr_metadata)
    if not species:
        print("Species not found in repertoire metadata file. Setting species as NCBITAXON:9606 aka human for as default.")
        species = "NCBITAXON:9606" #human
        # species = 'NCBITAXON:10090' #mouse
    
    if locus not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code
        
    print("=======================================================================================")
    print("                                    Parameters                                         ")
    print("=======================================================================================")
    print(f"\t\tlocus:              {locus}")
    print(f"\t\tspecies:            {species}")
    print(f"\t\tdata_dir:           {data_dir}")
    print(f"\t\tinput_file:         {input_file}")
    print(f"\t\toutput_file_base:   {output_file_base}")
    print("=======================================================================================")
    
    print("=======================================================================================")
    print("                                 Loading AIRR file.                                    ")
    print("=======================================================================================")
    
    print("Extracting junction_aa sequences...")
    junction_aa_list = airr_df["junction_aa"].dropna().tolist()
    all_unique_junction_aa = list(set(junction_aa_list))

    print("Pre-filtering junction_aa not in AKC DB...")
    # get available junction_aa in the database
    unique_junction_aa = set()
    
    for row in query_database_stream(all_unique_junction_aa, "chain", SQLITE_DB_PATH, locus, species):
        j_aa = row['junction_aa']
        unique_junction_aa.add(j_aa)
        
    unique_junction_aa = list(unique_junction_aa)
    
    print("Pre-calculating duplicate counts...")
    dup_counts = airr_df.groupby("junction_aa")["duplicate_count"].sum().to_dict()
    j_aa_freqs = airr_df["junction_aa"].value_counts().to_dict()

    print(f"Total productive sequences: {len(airr_df)}")
    print(f"Total Unique junction_aa in the airr file: {len(all_unique_junction_aa)}")
    print(f"Total Unique junction_aa in AKC being queried: {len(unique_junction_aa)}")

    #deleting the airr file to save space
    del airr_df
    gc.collect()
    
    print("=======================================================================================")
    print("                Querying Database for junction_aa and Epitope match.                   ")
    print("=======================================================================================")
    
    
    matched_columns = ['akc_assay_akc_id', 'akc_complex_akc_id', 'akc_epitope_akc_id', 'akc_epitope_seq_aa', 'akc_source_protein', 
                       'akc_source_organism', 'junction_aa', 'akc_species', 'akc_v_call', 'akc_j_call']
    
    detailed_tsv = f"{output_dir}/{output_file_base}.tilde.detail.tsv"
    summary_data = {}
    unique_assay_ids = set()
    processed_rows = 0
    last_printed_progress = 0
    # print(f"Writing detailed output to a tsv file")
    with open(detailed_tsv, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(matched_columns)

        total_junctions = len(unique_junction_aa)
        print("Total rows: ", total_junctions)
        # Counter for how many rows have been processed
        start_time = time.time()
        # Iterate through the generator
        for row in query_database_stream(unique_junction_aa, "junction_aa", SQLITE_DB_PATH, locus, species):
            # Write to detailed file immediately
            writer.writerow(row)

            # Track unique assays for the second query
            if row['akc_assay_akc_id']:
                unique_assay_ids.add(row['akc_assay_akc_id'])

            # Accumulate summary statistics in a memory-efficient dict
            j_aa = row['junction_aa']
            if j_aa not in summary_data:
                summary_data[j_aa] = {
                    "n_row_matches_akc_db": 0,
                    "unique_epitope_id": set(),
                    "unique_epitope_seq": set(),
                    "unique_orgs": set(),
                    "unique_proteins": set(),
                    "junction_repeat_count": j_aa_freqs.get(j_aa, 0),
                    "junction_total_dup_count": dup_counts.get(j_aa, 0)
                }
            
            s = summary_data[j_aa]
            s["n_row_matches_akc_db"] += 1
            if row['akc_epitope_akc_id']:s["unique_epitope_id"].add(row['akc_epitope_akc_id'])
            if row['akc_epitope_seq_aa']:s["unique_epitope_seq"].add(row['akc_epitope_seq_aa'])
            if row['akc_source_organism']:s["unique_orgs"].add(row['akc_source_organism'])
            if row['akc_source_protein']:s["unique_proteins"].add(row['akc_source_protein'])
            
            processed_junctions = len(summary_data)
            progress_percentage = (processed_junctions / total_junctions) * 100
            if progress_percentage >= last_printed_progress + 10 or progress_percentage == 100.0:
                elapsed = (time.time() - start_time) / 60
                sys.stdout.write(
                    f"\rProgress: {progress_percentage:.2f}% | "
                    f"{processed_junctions}/{total_junctions} junctions | Elapsed: {elapsed:.2f} minutes"
                )
                sys.stdout.flush()
                last_printed_progress = int(progress_percentage // 10) * 10
    if len(summary_data) < 1:
        print("--------------------------------------------------------------------------------------------")
        print(f"No matching records for:\ninput: {input_file}\nLocus: {locus}\nSpecies: {species}")
        print("--------------------------------------------------------------------------------------------")
        return
    print(f"\nDone writing detailes to {detailed_tsv} file")

    # Finalize Summary Dataframe
    summary_rows = []
    for j_aa, s in summary_data.items():
        summary_rows.append({
            "query_cdr3": j_aa,
            "n_row_matches_akc_db": s["n_row_matches_akc_db"],
            "n_unique_epitope_id": len(s["unique_epitope_id"]),
            "n_unique_epitope_seq": len(s["unique_epitope_seq"]),
            "unique_epitope_seq": ",".join(sorted(s["unique_epitope_seq"])),
            "unique_source_organism": ",".join(sorted(s["unique_orgs"])),
            "unique_source_protein": ",".join(sorted(s["unique_proteins"])),
            "junction_aa_repeat_count": s["junction_repeat_count"],
            "junction_total_dup_count": s["junction_total_dup_count"]
        })
    summary_df = pd.DataFrame(summary_rows).sort_values(by='n_unique_epitope_seq', ascending=False)
    summary_df = summary_df[summary_df.n_unique_epitope_id>0].reset_index(drop = True)

    print("\nSummary Dataframe \n")
    print(f"Total Unique Junction_aa match/Summary dataframe Shape: {summary_df.shape}")
    print(summary_df.head())
    print("\nWriting the junction_aa summary into a tsv file...")
    summary_df.to_csv(f"{output_dir}/{output_file_base}.tilde.summary.tsv", sep="\t", index=False)
    print("\nCreating plots for top junction_aa's\n")
    plot_cdr3_vs_epitope_stats(summary_df, output_file_base, output_dir, n = top_n_epitopes)
    
    print("=======================================================================================")
    print("                                   Querying assay objects...                         \n")
    print("=======================================================================================")
    
    # Query returns a generator
    assay_rows_gen = query_database_stream(list(unique_assay_ids), "assay", SQLITE_DB_PATH, locus, species)
    # Process the generator
    assay_df, all_assay_dict = process_query_results_stream(assay_rows_gen)
    print("\nAssay Dataframe:\n")
    print(f"Total unique number of Assay/Dataframe shape: {assay_df.shape}")
    print(f"\nType of Assay and their count: {assay_df.data_type.value_counts()}")
    print("\nWriting the assay objects into a json file...\n")
    json_filename = f"{output_file_base}.tilde.assay.json"
    # Serialize to file only at the very end
    with open(f'{output_dir}/{json_filename}', 'w') as f:
        json.dump(all_assay_dict, f, indent=4)
    print("=======================================================================================")
    print("=======================================================================================")
    print("                                     Analysis Complete!                                ")
    print("=======================================================================================")
        
if __name__ == "__main__":
    # find_index_on_a_table_sqlite(SQLITE_DB_PATH_1)
    main(SQLITE_DB_PATH)
    
#TILDE: TCR/Ig Linkage via CDR3 similarity for Discovery of Epitopes
#TILDE (TCR/Ig Linkage via similarity metrics for Discovery of Epitopes)