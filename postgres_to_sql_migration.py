import pandas as pd
import json
import logging
import time
from sqlalchemy import create_engine, inspect, Table, MetaData, text, JSON
from sqlalchemy.dialects.postgresql import JSONB
from dotenv import dotenv_values
import psycopg
import sqlite3

# Setup Logging & Config
config = dotenv_values(".env")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

DB_CONFIG = {
    "host": config["POSTGRES_HOST"],
    "dbname": config["POSTGRES_DB"],
    "user": config["POSTGRES_USER"],
    "password": config["POSTGRES_PASSWORD"],
}

def convert_to_sqlite(TARGET_TABLES, db_name='airrkb_test.db'):
    # Engines
    pg_url = f"postgresql+psycopg2://{DB_CONFIG['user']}:{DB_CONFIG['password']}@{DB_CONFIG['host']}/{DB_CONFIG['dbname']}"
    pg_engine = create_engine(pg_url)
    sqlite_engine = create_engine(f"sqlite:///{db_name}")

    # TARGET_TABLES = ["joined_tcr_data", "QueryAssay"]
    TARGET_TABLES = TARGET_TABLES + ["QueryAssay"]
    print("TARGET_TABLES: ", TARGET_TABLES)
    
    # --- PERFORMANCE: Enable SQLite Turbo Mode ---
    with sqlite_engine.connect() as conn:
        conn.execute(text("PRAGMA synchronous = OFF;"))
        conn.execute(text("PRAGMA journal_mode = OFF;"))
        conn.execute(text("PRAGMA cache_size = 100000;"))
        conn.execute(text("PRAGMA temp_store = MEMORY;"))

    inspector = inspect(pg_engine)
    pg_meta = MetaData()
    # tables = inspector.get_table_names()
    # print("Tables: ", tables)
    tables = [t for t in inspector.get_table_names() if t in TARGET_TABLES]
    print("Table Names: ", tables)

    # Store index definitions to run them ALL at the very end
    all_index_queries = []

    for table_name in tables:
        logger.info(f"--- Processing Table: {table_name} ---")
        try:
            table_schema = Table(table_name, pg_meta, autoload_with=pg_engine)

            # Identify JSON columns ONCE per table
            json_cols = []
            for column in table_schema.columns:
                if isinstance(column.type, JSONB):
                    column.type = JSON() # Prepare for SQLite
                    json_cols.append(column.name)
                    
            # Collect index queries for later
            indexes = inspector.get_indexes(table_name)
            for idx in indexes:
                # Skip PKs because they are created with the table automatically
                if not idx.get('primary_key'):
                    cols = ", ".join([f'"{c}"' for c in idx['column_names']])
                    unique = "UNIQUE" if idx['unique'] else ""
                    # This builds the string and appends it to our master list
                    idx_sql = f"CREATE {unique} INDEX IF NOT EXISTS {idx['name']} ON {table_name} ({cols})"
                    all_index_queries.append(idx_sql)

            table_schema.drop(sqlite_engine, checkfirst=True)
            table_schema.create(sqlite_engine)

            query = f'SELECT * FROM "{table_name}"'
            for i, chunk in enumerate(pd.read_sql(query, pg_engine, chunksize=200000)):
                
                # Only run this if we actually found JSON columns in this specific table
                if json_cols:
                    for col in json_cols:
                        chunk[col] = chunk[col].apply(
                            lambda x: json.dumps(x) if x is not None else None
                        )

                chunk.to_sql(table_name, sqlite_engine, if_exists='append', index=False)
                logger.info(f"Table {table_name}: Chunk {i+1} transferred.")

        except Exception as e:
            logger.error(f"Failed to migrate {table_name}: {str(e)}")

    # Bulk Index Creation ---
    logger.info(f"Starting index creation for {len(all_index_queries)} indexes...")
    start_idx_time = time.time()
    
    with sqlite_engine.begin() as conn:
        for idx_sql in all_index_queries:
            try:
                conn.execute(text(idx_sql))
                logger.info(f"Executed: {idx_sql[:50]}...")
            except Exception as e:
                logger.warning(f"Could not create index: {str(e)}")

    logger.info(f"Index creation took {time.time() - start_idx_time:.2f} sec")
    logger.info("Migration Fully Complete.")



def prejoin_table(locus, species, table_name):
    with psycopg.connect(**DB_CONFIG) as conn:
        with conn.cursor(row_factory=psycopg.rows.dict_row) as cur:
            start_time = time.time()
            cur.execute(f"DROP TABLE IF EXISTS {table_name};")
            query = f"""
            CREATE TABLE {table_name} AS
            WITH chain_filtered AS (
                SELECT *
                FROM "Chain"
                WHERE species = '{species}'
            )
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
            JOIN chain_filtered ch
                ON t.{locus}_chain = ch.akc_id
            LEFT OUTER JOIN "Epitope" e
                ON c.epitope = e.akc_id
            """
            cur.execute(query)
            #Create index here
            cur.execute(f"CREATE INDEX idx_{table_name}_junction ON {table_name} (junction_aa);")
            cur.execute(f"SELECT COUNT(DISTINCT junction_aa) AS unique_cdr3 FROM {table_name};")
            row = cur.fetchone()
            print("Unique junction_aa:", row["unique_cdr3"])
            conn.commit()

            print(f"Total time for the query: {time.time() - start_time:.2f} sec\n")
            
    return  


def prejoin_tables():
    table_names = []
    # locus_list = ['tra']
    # species_list = ['NCBITAXON:10090']
    locus_list = ['tra', 'trb', 'trd', 'trg']
    species_list = ['NCBITAXON:9606', 'NCBITAXON:10090']

    for locus in locus_list:
        for species in species_list:
            species_renamed = species.replace(':', '_')
            table_name = f"{locus}_{species_renamed}_tilde".lower()
            
            prejoin_table(locus, species, table_name)
            print(f"Done joining table {table_name}")
            table_names.append(table_name)
    return table_names
        
# --------------------------------------------------------------------------------- 
# Find tables that has index on it
# --------------------------------------------------------------------------------- 

def find_index_on_a_table_sqlite(SQLITE_DB_PATH):
    # Connect to the SQLite database
    db_path = f"file:{SQLITE_DB_PATH}?mode=ro&immutable=1"
    conn = sqlite3.connect(db_path, uri=True)
    # conn = sqlite3.connect(SQLITE_DB_PATH)
    cur = conn.cursor()
    # table_names = ['TCRpMHCComplex', 'Assay_tcr_complexes', 'TCellReceptor', 'Chain', 'Epitope', 'QueryAssay']
    table_names = []
    locus_list = ['tra', 'trb', 'trd', 'trg']
    species_list = ['NCBITAXON:9606', 'NCBITAXON:10090']
    for locus in locus_list:
        for species in species_list:
            species_renamed = species.replace(':', '_')
            table_name = f"{locus}_{species_renamed}_tilde".lower()
            table_names.append(table_name)
    table_names_tuple = tuple(table_names)
    query = """
    SELECT name, tbl_name, sql
    FROM sqlite_master
    WHERE type = 'index' AND tbl_name IN ({})
    """.format(', '.join('?' for _ in table_names_tuple))  # Use parameterized queries

    # Execute the query with the list of tables as arguments
    cur.execute(query, table_names_tuple)

    # Fetch and print the indexes
    indexes = cur.fetchall()
    if indexes:
        for index in indexes:
            print(f"Index Name: {index[0]}")
            print(f"Table Name: {index[1]}")
            print(f"SQL: {index[2]}\n")
    else:
        print("No indexes found for the specified tables.")

    # Close the connection
    conn.close()


if __name__ == "__main__":
    total_start = time.time()
    # db_filename = 'airrkb_v2_optimized.db'
    db_filename = 'airrkb_v2_tilde.db'
    
    # table_names = prejoin_tables()
    # print(f"\n Total Table Prejoin Time: {time.time() - total_start:.2f} seconds")
    
    # convert_to_sqlite(table_names, db_filename)
    
    find_index_on_a_table_sqlite(db_filename)
    print(f"\n Total Migration Time: {time.time() - total_start:.2f} seconds")
    