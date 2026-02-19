import psycopg
import sys
import re
import argparse

DB_CONFIG = {
    "host": "ak-db",
    "dbname": "airrkb_v1",
    "user": "postgres",
    "password": "example",
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
    parser.add_argument("table_name", help="Name of the table to create")
    parser.add_argument("locus_type", help="Locus type (TCR or BCR)")
    parser.add_argument("chain_type", help="Chain type (tra, trb, igh, igk etc.)")
    parser.add_argument("--species", default="NCBITAXON:9606", help="Species CURIE (default: NCBITAXON:9606)")
    
    args = parser.parse_args()
    table_name = args.table_name
    locus_type = args.locus_type
    chain_type = args.chain_type
    species = args.species

    if not VALID_NAME.match(table_name):
        raise ValueError("Invalid table name")
    
    with psycopg.connect(**DB_CONFIG) as conn:
        with conn.cursor() as cur:
            # Drop table if exists
            print(f"Dropping table {table_name} if it exists...")
            cur.execute(DROP_SQL_QUERY.format(table_name=table_name))
            
            # Create final table
            print(f"Creating table {table_name}...")
            cur.execute(CREATE_TABLE_SQL.format(table_name=table_name))
            
            # Create temporary staging table
            print("Creating temporary table for raw junctions...")
            cur.execute(CREATE_TEMP_TABLE_SQL)
            
            # Populate temp table
            if locus_type == "TCR":
                sql = POPULATE_TABLE_TCR_SQL.format(chain_type=chain_type)
            elif locus_type == "BCR":
                sql = POPULATE_TABLE_BCR_SQL.format(chain_type=chain_type)
            else:
                raise ValueError("Invalid Locus type, must be TCR or BCR")
            
            print("Populating temporary table with junction_aa...")
            if DEBUG:
                print(sql)
            cur.execute(sql, {"species": species})
            
            # Build final table with dense IDs
            print("Building final table with dense consecutive IDs...")
            cur.execute(BUILD_FINAL_TABLE_SQL.format(table_name=table_name))
            
        
        conn.commit()
    
    print("Done. Temporary table junctions_raw_temp is automatically dropped.")

if __name__ == "__main__":
    main()