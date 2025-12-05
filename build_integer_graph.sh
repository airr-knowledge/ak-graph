#!/bin/bash
set -euo pipefail

WORKDIR="graph_data"
mkdir -p "$WORKDIR"

PAIR_FILE="pairs_with_indels.tsv"
EDGE_TMP="$WORKDIR/edges_tmp.tsv"
EDGE_FILE="$WORKDIR/edges_id_with_indels_updated.tsv"
MAP_FILE="$WORKDIR/cdr3_to_id_with_indels.tsv"
LOG_FILE="$WORKDIR/pipeline_2.log"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "$LOG_FILE"; }

log "Pipeline started"

# Step 1: Assign IDs and generate edges
START=$(date +%s)
log "Step 1: Assigning IDs and generating edge list (excluding self-loops)..."

awk -F'\t' '
BEGIN { id_counter=0 }
NR>1 {
    n1 = $6   # junction_aa_1
    n2 = $12  # junction_aa_2

    # Assign integer IDs to all nodes
    if(!(n1 in seq_to_id)) { seq_to_id[n1]=id_counter; print id_counter "\t" n1 >> "'"$MAP_FILE"'"; id_counter++ }
    if(!(n2 in seq_to_id)) { seq_to_id[n2]=id_counter; print id_counter "\t" n2 >> "'"$MAP_FILE"'"; id_counter++ }

    # Skip self-loops
    if(n1==n2) next

    id1=seq_to_id[n1]
    id2=seq_to_id[n2]

    # Canonicalize edges
    if(id1 < id2) print id1 "\t" id2 >> "'"$EDGE_TMP"'"
    else print id2 "\t" id1 >> "'"$EDGE_TMP"'"
}
' "$PAIR_FILE"

END=$(date +%s)
log "Step 1 complete in $((END-START)) seconds. Temp edge file: $EDGE_TMP, mapping file: $MAP_FILE"

# Step 2: Safe sort + deduplication
START=$(date +%s)
log "Step 2: Sorting and deduplicating edges safely..."


mkdir -p graph_data/tmp
export TMPDIR="graph_data/tmp"

LC_ALL=C sort -u --parallel=8 -S 500G "$EDGE_TMP" > "$EDGE_FILE"

# rm "$EDGE_TMP"
END=$(date +%s)
log "Step 2 complete in $((END-START)) seconds. Final edge file: $EDGE_FILE"

# Step 3: Count nodes and edges
NUM_NODES=$(wc -l < "$MAP_FILE")
NUM_EDGES=$(wc -l < "$EDGE_FILE")
log "Total nodes: $NUM_NODES"
log "Total edges: $NUM_EDGES"

log "Pipeline finished successfully"
