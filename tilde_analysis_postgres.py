import psycopg
import sys
import re
import argparse
from dotenv import dotenv_values
import os
import time
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import json
import csv
matplotlib.use("Agg")

config = dotenv_values(".env")

DB_CONFIG = {
    "host": config['POSTGRES_HOST'],
    "dbname": config['POSTGRES_DB'],
    "user": config['POSTGRES_USER'],
    "password": config['POSTGRES_PASSWORD'],
}

def get_query_for_locus(locus):
    query = f"""
    SELECT
        atc.assay_akc_id AS akc_assay_akc_id,
        c.akc_id AS akc_complex_akc_id,
        c.epitope AS akc_epitope_akc_id,
        e.sequence_aa AS akc_epitope_seq_aa,
        e.source_protein AS akc_source_protein,
        e.source_organism AS akc_source_organism,
        ch.junction_aa AS junction_aa,
        ch.species AS akc_species,
        ch.v_call AS akc_v_call,
        ch.j_call AS akc_j_call
    FROM "TCRpMHCComplex" c
    JOIN "Assay_tcr_complexes" atc
        ON c.akc_id = atc.tcr_complexes_akc_id
    JOIN "TCellReceptor" t
        ON c.tcr = t.akc_id
    JOIN "Chain" ch
        ON t.{locus}_chain = ch.akc_id
    LEFT OUTER JOIN "Epitope" e
        ON c.epitope = e.akc_id
    WHERE ch.junction_aa = ANY(%s)
    """
    return query

def get_query_for_assay_object():

    query = f"""
    SELECT 
        *
    FROM "QueryAssay" qa
    WHERE qa.akc_id = ANY(%s)
    """
    return query

def load_airr_file(filepath):
    """Load AIRR TSV and extract junction_aa."""
    df = pd.read_csv(filepath, sep="\t", low_memory=False)
    df = df[(df.productive == 'T') | (df.productive == True)]
    if "junction_aa" not in df.columns:
        raise ValueError("AIRR file must contain 'junction_aa' column")
    
    df = df[['sequence_id', 'junction_aa', 'duplicate_count']]
    print(df.head())

    return df

# def query_database(parameter, query_type, locus=None):
#     with psycopg.connect(**DB_CONFIG) as conn:
#         with conn.cursor(row_factory=psycopg.rows.dict_row) as cur:
#             start_time = time.time()
            
#             # Choose query based on query_type
#             if query_type == "junction_aa":
#                 if locus is None:
#                     raise ValueError("Locus must be provided for CDR3 query.")
#                 query = get_query_for_locus(locus)
#                 cur.execute(query, (parameter,))
#             elif query_type == "assay":
#                 query = get_query_for_assay_object()
#                 cur.execute(query, (parameter,))
#             else:
#                 raise ValueError(f"Unknown query_type: {query_type}")
            
#             rows = cur.fetchall()
#             print(f"Total time for the {query_type} query: {time.time() - start_time:.2f} sec\n")
#             print(f"Total {query_type} rows returned: {len(rows)}")
            
#     return rows

def query_database_stream(parameter, query_type, locus=None):
    """Yields rows one by one to save memory."""
    with psycopg.connect(**DB_CONFIG) as conn:
        # dict_row makes the cursor return dictionaries automatically
        with conn.cursor(row_factory=psycopg.rows.dict_row) as cur:
            if query_type == "junction_aa":
                if locus is None:
                    raise ValueError("Locus must be provided.")
                query = get_query_for_locus(locus)
                cur.execute(query, (parameter,))
            elif query_type == "assay":
                query = get_query_for_assay_object()
                cur.execute(query, (parameter,))
            
            # Streaming: Yield each row as it comes from the DB
            for row in cur:
                yield row
                
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
        assay_raw = row.get("assay_object")

        if not assay_raw:
            continue

        # Handle JSON parsing only if necessary
        if isinstance(assay_raw, str):
            assay_dict = json.loads(assay_raw)
        elif isinstance(assay_raw, dict):
            assay_dict = assay_raw
        else:
            continue

        assay_id = row.get('akc_id')
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
    # We will serialize it to JSON in the main block.
    return assay_df, all_assay_dict
    

def plot_top_junction_aa(summary_df, output_dir, n = 15):
    temp = summary_df.head(n)
    fig, axes = plt.subplots(1, 2, figsize = (10, 6), gridspec_kw={'width_ratios': [2, 3]})
    sns.barplot(data = temp, y = 'query_cdr3', x = 'n_unique_epitope_seq', ax = axes[0])
    sns.scatterplot(data = summary_df, x = 'n_unique_epitope_id', y = 'n_unique_epitope_seq', ax = axes[1])
    sns.regplot(
        data=summary_df,
        x='n_unique_epitope_id',
        y='n_unique_epitope_seq',
        scatter=False,
        ax=axes[1],
        color='red'
    )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/postgres_summary_distribution_top_junction_aa.png", dpi=300)
    plt.close()

def create_directories_if_not_exist(path):
    """Create directories if they do not exist"""
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    parser = argparse.ArgumentParser("Create and populate junction table.")
    
    parser.add_argument("--LOCUS", default="trb", help="Locus [tra, trb, trd, trg, igh, igk, igl]")
    parser.add_argument("--SPECIES", default="NCBITAXON:9606", help="Species CURIE (default: NCBITAXON:9606)")
    parser.add_argument("--VERSION", default="v1", help="Version of the table name that will be put on the table")
    parser.add_argument("--DATA_DIR", default="./local_ak_graph_data/", help="Version of the table name that will be put on the graph")
    parser.add_argument("--INPUT_FILE", default= "test_input.airr.tsv", help="Version of the table name that will be put on the graph")
    parser.add_argument("--OUTPUT_FILE", default="test_output.tsv", help="Version of the table name that will be put on the graph")
    
    args = parser.parse_args()
    LOCUS = args.LOCUS.lower()
    SPECIES = args.SPECIES
    VERSION = args.VERSION
    DATA_DIR = args.DATA_DIR
    INPUT_FILE = args.INPUT_FILE
    OUTPUT_FILE = args.OUTPUT_FILE

    
    OUTPUT_DIR = f"{DATA_DIR}/postgres/tilda_output"
    
    #create figure directory if not exist
    create_directories_if_not_exist(OUTPUT_DIR)
    
    if LOCUS not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code

    
    print("=======================================================================================")
    print("                                        Parameters                                     ")
    print("=======================================================================================")
    print(f"\t\tLOCUS:          {LOCUS}")
    print(f"\t\tSPECIES:        {SPECIES}")
    print(f"\t\tVERSION:        {VERSION}")
    print(f"\t\tDATA_DIR:       {DATA_DIR}")
    print(f"\t\tOUTPUT_DIR:     {OUTPUT_DIR}")
    print(f"\t\tINPUT_FILE:     {INPUT_FILE}")
    print(f"\t\tOUTPUT_FILE:    {OUTPUT_FILE}")
    print("=======================================================================================")
    
    
    output_basename = OUTPUT_FILE.strip().split('.')[0]
    
    print("                                       Loading AIRR file...                            ")
    print("=======================================================================================")
    
    airr_df = load_airr_file(INPUT_FILE)

    print("Extracting junction_aa sequences...")
    
    junction_aa_list = airr_df["junction_aa"].dropna().tolist()
    unique_junction_aa = list(set(junction_aa_list))
    
    print("Pre-calculating duplicate counts...")

    dup_counts = airr_df.groupby("junction_aa")["duplicate_count"].sum().to_dict()

    j_aa_freqs = airr_df["junction_aa"].value_counts().to_dict()
    
    print(f"Total productive sequences: {len(airr_df)}")
    print(f"Total Unique junction_aa being queried: {len(unique_junction_aa)}")
    
    print("=======================================================================================")
    
    print("                              Querying Database for junction_aa and Epitope match...\n")
    print("=======================================================================================")
    
    
    matched_columns = ['akc_assay_akc_id', 'akc_complex_akc_id', 'akc_epitope_akc_id', 'akc_epitope_seq_aa', 'akc_source_protein', 
                       'akc_source_organism', 'junction_aa', 'akc_species', 'akc_v_call', 'akc_j_call']
    
    detailed_tsv = f"{OUTPUT_DIR}/sqlite_{output_basename}_junction_aa_matches_detailed.tsv"
    summary_data = {}
    unique_assay_ids = set()
    
    print(f"Writing detailed output to a tsv file")
    with open(detailed_tsv, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=matched_columns, delimiter='\t')
        writer.writeheader()

        # Iterate through the generator
        for row in query_database_stream(unique_junction_aa, "junction_aa", LOCUS):
            # Write to detailed file immediately
            writer.writerow(row)

            # Track unique assays for the second query
            if row['akc_assay_akc_id']:
                unique_assay_ids.add(row['akc_assay_akc_id'])

            # Accumulate summary statistics in a memory-efficient dict
            j_aa = row['junction_aa']
            if j_aa not in summary_data:
                summary_data[j_aa] = {
                    "n_row_matches": 0,
                    "unique_epitope_id": set(),
                    "unique_epitope_seq": set(),
                    "unique_orgs": set(),
                    "unique_proteins": set(),
                    "junction_repeat_count": j_aa_freqs.get(j_aa, 0),
                    "junction_total_dup_count": dup_counts.get(j_aa, 0)
                }
            
            s = summary_data[j_aa]
            s["n_row_matches"] += 1
            if row['akc_epitope_akc_id']:s["unique_epitope_id"].add(row['akc_epitope_akc_id'])
            if row['akc_epitope_seq_aa']:s["unique_epitope_seq"].add(row['akc_epitope_seq_aa'])
            if row['akc_source_organism']:s["unique_orgs"].add(row['akc_source_organism'])
            if row['akc_source_protein']:s["unique_proteins"].add(row['akc_source_protein'])
    print(f"Done writing detailes to {detailed_tsv} file")

    # Finalize Summary Dataframe
    summary_rows = []
    for j_aa, s in summary_data.items():
        summary_rows.append({
            "query_cdr3": j_aa,
            "n_row_matches": s["n_row_matches"],
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
    summary_df.to_csv(f"{OUTPUT_DIR}/sqlite_{output_basename}_junction_aa_matches_summary.tsv", sep="\t", index=False)
    print("\nCreating plots for top junction_aa's\n")
    plot_top_junction_aa(summary_df, OUTPUT_DIR)
    print("=======================================================================================")
    
    print("                                   Querying assay objects...                         \n")
    print("=======================================================================================")
    
    # Query returns a generator
    assay_rows_gen = query_database_stream(list(unique_assay_ids), "assay")

    # Process the generator
    assay_df, all_assay_dict = process_query_results_stream(assay_rows_gen)
    

    print("\nAssay Dataframe:\n")
    print(f"Total unique number of Assay/Dataframe shape: {assay_df.shape}")
    print(f"\nType of Assay and their count: {assay_df.data_type.value_counts()}")
    print("\nWriting the assay objects into a json file...\n")
    json_filename = f"{output_basename}_all_assay_dict.json"
    # Serialize to file only at the very end
    with open(f'{OUTPUT_DIR}/sqlite_{json_filename}', 'w') as f:
        json.dump(all_assay_dict, f, indent=4)
    print("=======================================================================================")
    
    print("=======================================================================================")
    print("                                     Analysis Complete!                                ")
    print("=======================================================================================")
    
if __name__ == "__main__":
    main()