import psycopg
import argparse
from dotenv import dotenv_values
import pandas as pd
import os

config = dotenv_values(".env")

DB_CONFIG = {
    "host": config["POSTGRES_HOST"],
    "dbname": config["POSTGRES_DB"],
    "user": config["POSTGRES_USER"],
    "password": config["POSTGRES_PASSWORD"],
}

def get_participant_counts(species=None):
    base_query = """
        SELECT sex, COUNT(*) FROM "Participant"
    """
    params = []

    if species:
        base_query += " WHERE species = %s"
        params.append(species)

    base_query += " GROUP BY sex"

    return base_query, params

def get_cdr3_and_epitopes(locus):
    params = []
    if locus == 'both':
        base_query = """
        SELECT
            e.sequence_aa,
            e.source_protein,
            e.source_organism,
            ch.junction_aa,
            ch.locus,
            ch.species,
            a.measurement_category
        FROM "TCRpMHCComplex" c
        JOIN "TCellReceptor" t
            ON c.tcr = t.akc_id
        JOIN "Chain" ch
            ON t.tra_chain = ch.akc_id
        JOIN "Epitope" e
            ON c.epitope = e.akc_id
        JOIN "Assay" a
            ON a.epitope = e.akc_id

        UNION ALL

        SELECT
            e.sequence_aa,
            e.source_protein,
            e.source_organism,
            ch.junction_aa,
            ch.locus,
            ch.species,
            a.measurement_category
        FROM "TCRpMHCComplex" c
        JOIN "TCellReceptor" t
            ON c.tcr = t.akc_id
        JOIN "Chain" ch
            ON t.trb_chain = ch.akc_id
        JOIN "Epitope" e
            ON c.epitope = e.akc_id
        JOIN "Assay" a
            ON a.epitope = e.akc_id
        """
    else: 
        base_query = f"""
        SELECT
            e.sequence_aa,
            e.source_protein,
            e.source_organism,
            ch.junction_aa,
            ch.locus,
            ch.species,
            a.measurement_category
        FROM "TCRpMHCComplex" c
        JOIN "TCellReceptor" t
            ON c.tcr = t.akc_id
        JOIN "Chain" ch
            ON t.{locus}_chain = ch.akc_id
        JOIN "Epitope" e
            ON c.epitope = e.akc_id
        JOIN "Assay" a
            ON a.epitope = e.akc_id
        """
    return base_query, params


def get_epitope_counts(species=None):
    base_query = """
        SELECT sequence_aa, count(*) FROM "Epitope"
    """
    params = []
    base_query += " GROUP BY sequence_aa"
    return base_query, params

def create_directories_if_not_exist(path):
    """Create directories if they do not exist"""
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("LOCUS", default = 'tra', help="Locus to get the epitope data (tra, trb, both)")
    parser.add_argument("--DATA_DIR", default="/ak_graph_data/", help="Version of the table name that will be put on the graph")
    parser.add_argument("--VERSION", default="v1", help="Version of the table name that will be put on the graph")
    
    args = parser.parse_args()
    LOCUS = args.LOCUS
    DATA_DIR = args.DATA_DIR
    VERSION = args.VERSION
    
    if LOCUS not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code

    query, params = get_cdr3_and_epitopes(LOCUS)
    FILE_PATH = f"{DATA_DIR}/epitope_data/{LOCUS}_cdr3_epitope_info_{VERSION}.parquet"
    # local_fig_dir = f'./ak_graph_data_v2/epitope_data'
    
    # Create the necessary directories if they do not exist
    create_directories_if_not_exist(f"{DATA_DIR}/epitope_data")
    # create_directories_if_not_exist(f"{local_fig_dir}/epitope_data")
    
    print("================================================")
    print("                   Parameters                   ")
    print("================================================")
    print(f"LOCUS:              {LOCUS}")
    print(f"VERSION:            {VERSION}")
    print(f"DATA_DIR:           {DATA_DIR}")
    print("================================================")
    
    with psycopg.connect(**DB_CONFIG) as conn:
        with conn.cursor() as cur:
            print("Executing SQL:\n", query)
            cur.execute(query, params)
            rows = cur.fetchall()
            print("Total number of rows returned from the query: ", len(rows))
            columns = ['akc_epitope_seq_aa', 'akc_source_protein', 'akc_source_organism', 'junction_aa', 'locus', 'akc_species', 'akc_measurement_category']
            df = pd.DataFrame(rows, columns=columns)
            df.to_parquet(FILE_PATH, index=False)
            print(df.head())
        conn.commit()

if __name__ == "__main__":
    main()

