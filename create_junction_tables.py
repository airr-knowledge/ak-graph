import psycopg
import sys
import re
import argparse
from dotenv import dotenv_values

config = dotenv_values(".env")

DB_CONFIG = {
    "host": config['POSTGRES_HOST'],
    "dbname": config['POSTGRES_DB'],
    "user": config['POSTGRES_USER'],
    "password": config['POSTGRES_PASSWORD'],
}

# To only allow simple SQL identifiers
VALID_NAME = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")
DEBUG = True

# SQL templates
DROP_SQL_QUERY = """
    DROP TABLE IF EXISTS {table_name} CASCADE;
"""

# No BIGSERIAL here — IDs will be generated with ROW_NUMBER later
CREATE_TABLE_SQL = """
    CREATE TABLE {table_name} (
        sequence_id BIGINT PRIMARY KEY,
        junction_aa TEXT NOT NULL UNIQUE
    );
"""

# Temp table for staging
CREATE_TEMP_TABLE_SQL = """
    CREATE TEMP TABLE junctions_raw_temp (
        junction_aa TEXT NOT NULL
    );
"""

POPULATE_TABLE_TCR_SQL = """
    INSERT INTO junctions_raw_temp (junction_aa)
    SELECT c.junction_aa
    FROM "TCellReceptor" tcr
    JOIN "Chain" c ON tcr.{chain_type}_chain = c.akc_id
    WHERE c.junction_aa IS NOT NULL
    AND c.species = %(species)s;
"""

POPULATE_TABLE_BCR_SQL = """
    INSERT INTO junctions_raw_temp (junction_aa)
    SELECT c.junction_aa
    FROM "BCellReceptor" bcr
    JOIN "Chain" c ON bcr.{chain_type}_chain = c.akc_id
    WHERE c.junction_aa IS NOT NULL
    AND c.species = %(species)s;
"""

# Build final table with dense consecutive IDs
BUILD_FINAL_TABLE_SQL = """
    INSERT INTO {table_name} (sequence_id, junction_aa)
    SELECT ROW_NUMBER() OVER (ORDER BY junction_aa) - 1 AS sequence_id,
        junction_aa
    FROM (SELECT DISTINCT junction_aa FROM junctions_raw_temp) t;
"""

    
# QUERY_BCR = """
#     SELECT bcr.akc_id, c.junction_aa FROM "BCellReceptor" bcr INNER JOIN "Chain" c on bcr.{chain_type}_chain = c.akc_id WHERE c.species = %(species)s LIMIT 10
# """ 



def main():
    parser = argparse.ArgumentParser("Create and populate junction table.")
    
    parser.add_argument("LOCUS", help="Locus [tra, trb, trd, trg, igh, igk, igl]")
    parser.add_argument("--SPECIES", default="NCBITAXON:9606", help="Species CURIE (default: NCBITAXON:9606)")
    parser.add_argument("--VERSION", default="v1", help="Version of the table name that will be put on the table")
    
    args = parser.parse_args()
    LOCUS = args.LOCUS.lower()
    SPECIES = args.SPECIES
    VERSION = args.VERSION

    if LOCUS not in ['tra', 'trb', 'trd', 'trg', 'igh', 'igk', 'igl']:
        parser.print_help(sys.stderr) # Prints help message to standard error
        sys.exit(1) # Exit with an error code

    # create table name joining the LOCUS and VERSION
    TABLE_NAME = f'unique_junctions_{LOCUS}_{VERSION}'
    
    print("================================================")
    print("                   Parameters                   ")
    print("================================================")
    print(f"LOCUS:          {LOCUS}")
    print(f"SPECIES:        {SPECIES}")
    print(f"VERSION:        {VERSION}")
    print(f"TABLE-NAME:     {TABLE_NAME}")
    print("================================================")
    
    if not VALID_NAME.match(TABLE_NAME):
        raise ValueError("Invalid table name")
    
    with psycopg.connect(**DB_CONFIG) as conn:
        with conn.cursor() as cur:
            # Drop table if exists
            print(f"Dropping table {TABLE_NAME} if it exists...")
            cur.execute(DROP_SQL_QUERY.format(table_name=TABLE_NAME))
            
            # Create final table
            print(f"Creating table {TABLE_NAME}...")
            cur.execute(CREATE_TABLE_SQL.format(table_name=TABLE_NAME))
            
            # Create temporary staging table
            print("Creating temporary table for raw junctions...")
            cur.execute(CREATE_TEMP_TABLE_SQL)
            
            # Populate temp table
            if LOCUS in ['tra', 'trb', 'trd', 'trg']:
                sql = POPULATE_TABLE_TCR_SQL.format(chain_type=LOCUS)
            elif LOCUS in ['igh', 'igk', 'igl']:
                sql = POPULATE_TABLE_BCR_SQL.format(chain_type=LOCUS)
            
            print("Populating temporary table with junction_aa...")
            if DEBUG:
                print(sql)
            cur.execute(sql, {"species": SPECIES})
            
            # Build final table with dense IDs
            print("Building final table with dense consecutive IDs...")
            cur.execute(BUILD_FINAL_TABLE_SQL.format(table_name=TABLE_NAME))
            
        
        conn.commit()
    
    print("Done. Temporary table junctions_raw_temp is automatically dropped.")

if __name__ == "__main__":
    main()
