import pandas as pd
import json
import logging
import time
from sqlalchemy import create_engine, inspect, Table, MetaData, text, JSON
from sqlalchemy.dialects.postgresql import JSONB
from dotenv import dotenv_values

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

def convert_to_sqlite(db_name='airrkb_test.db'):
    # Engines
    pg_url = f"postgresql+psycopg2://{DB_CONFIG['user']}:{DB_CONFIG['password']}@{DB_CONFIG['host']}/{DB_CONFIG['dbname']}"
    pg_engine = create_engine(pg_url)
    sqlite_engine = create_engine(f"sqlite:///{db_name}")

    # --- PERFORMANCE: Enable SQLite Turbo Mode ---
    with sqlite_engine.connect() as conn:
        conn.execute(text("PRAGMA synchronous = OFF;"))
        conn.execute(text("PRAGMA journal_mode = OFF;"))
        conn.execute(text("PRAGMA cache_size = 100000;"))
        conn.execute(text("PRAGMA temp_store = MEMORY;"))

    inspector = inspect(pg_engine)
    pg_meta = MetaData()
    tables = inspector.get_table_names()

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
            for i, chunk in enumerate(pd.read_sql(query, pg_engine, chunksize=100000)):
                
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

if __name__ == "__main__":
    total_start = time.time()
    db_filename = 'airrkb_v2_optimized.db'
    
    convert_to_sqlite(db_filename)
    
    print(f"\n Total Migration Time: {time.time() - total_start:.2f} seconds")