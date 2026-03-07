import psycopg
import argparse
from dotenv import dotenv_values
import pandas as pd

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
        c.cdr3_aa,
        c.locus,
        e.sequence_aa,
        a.measurement_category
        FROM "Assay" a
        JOIN "Epitope" e
            ON a.epitope = e.akc_id
        JOIN "Assay_tcell_receptors" atr
            ON a.akc_id = atr.assay_akc_id
        JOIN "TCellReceptor" tcr
            ON atr.tcell_receptors_akc_id = tcr.akc_id
        JOIN "Chain" c
            ON c.akc_id = tcr.tra_chain

        UNION ALL

        SELECT
            c.cdr3_aa,
            c.locus,
            e.sequence_aa,
            a.measurement_category
        FROM "Assay" a
        JOIN "Epitope" e
            ON a.epitope = e.akc_id
        JOIN "Assay_tcell_receptors" atr
            ON a.akc_id = atr.assay_akc_id
        JOIN "TCellReceptor" tcr
            ON atr.tcell_receptors_akc_id = tcr.akc_id
        JOIN "Chain" c
            ON c.akc_id = tcr.trb_chain;
        """
    else: 
        base_query = f"""
        SELECT
            c.cdr3_aa,
            c.locus,
            e.sequence_aa,
            a.measurement_category
        FROM "Assay" a
        JOIN "Epitope" e
            ON a.epitope = e.akc_id
        JOIN "Assay_tcell_receptors" atr
            ON a.akc_id = atr.assay_akc_id
        JOIN "TCellReceptor" tcr
            ON atr.tcell_receptors_akc_id = tcr.akc_id
        JOIN "Chain" c
            ON tcr.{locus}_chain = c.akc_id
        """
    return base_query, params


def get_epitope_counts(species=None):
    base_query = """
        SELECT sequence_aa, count(*) FROM "Epitope"
    """
    params = []
    base_query += " GROUP BY sequence_aa"
    return base_query, params


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--LOCUS", default = 'tra', help="Locus to get the epitope data (tra, trb, both)")
    parser.add_argument("--WORKDIR", default="./ak_graph_data/", help="Version of the table name that will be put on the graph")
    parser.add_argument("--VERSION", default="v1", help="Version of the table name that will be put on the graph")
    
    args = parser.parse_args()
    LOCUS = args.LOCUS
    WORKDIR = args.WORKDIR
    VERSION = args.VERSION
    
    if LOCUS not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code

    query, params = get_cdr3_and_epitopes(LOCUS)
    FILE_PATH = f"{WORKDIR}/{LOCUS}_cdr3_epitope_info_{VERSION}.parquet"
    
    with psycopg.connect(**DB_CONFIG) as conn:
        with conn.cursor() as cur:
            print("Executing SQL:\n", query)
            cur.execute(query, params)
            rows = cur.fetchall()
            print("Total number of rows returned from the query: ", len(rows))
            df = pd.DataFrame(rows, columns=["cdr3_aa", "locus", "epitope", "measurement_category"])
            df.to_parquet(FILE_PATH, index=False)
        conn.commit()

if __name__ == "__main__":
    main()