
# AIRR Knowledge graph analyses

# config
include .env
PG_CONN=postgresql://$(POSTGRES_USER):$(POSTGRES_PASSWORD)@$(POSTGRES_HOST)/postgres
PG_AK_CONN=postgresql://$(POSTGRES_USER):$(POSTGRES_PASSWORD)@$(POSTGRES_HOST)/$(POSTGRES_DB)
PG_DISPLAY_CONN=postgresql://$(POSTGRES_USER):XXXXXX@$(POSTGRES_HOST)/$(POSTGRES_DB)
# database connection info
export PG_AK_CONN
# data for ak-etvl
export IMPORT_DATA
# data for ak-graph
export GRAPH_DATA

LOCUS ?= trb
VERSION ?= v1
DATA_DIR ?=  /ak_graph_data
SPECIES ?= NCBITAXON:9606
N_THREADS ?= 8

BIN_DIR = ./bin
OUTPUT_DIR = ./output
COMPAIRR_PROGS = stream_query stream_query_no_output stream_query_line stream_query_threads stream_query_threads_no_output
PROG_NAMES = stream_query stream_query_threads

CXX      := g++
CXXFLAGS := -Wall -O2
LDLIBS   := -lpqxx

help:
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "                          -- AIRR Knowledge graph analysis"
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "                                 Using DB: $(PG_DISPLAY_CONN)"
	@echo "       Host folder for ak-etvl (/ak_data): $(IMPORT_DATA)"
	@echo "Host folder for ak-graph (/ak_graph_data): $(GRAPH_DATA)"
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo ""
	@echo "make docker                 		-- Build docker image"
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "                             -- Create Junction_aa tables"
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo " (run within docker)"
	@echo "make create-table				-- generate unique CDR3 sequence tables using $(LOCUS), $(VERSION), $(SPECIES)"
	@echo "make create-table-tra			-- generate unique CDR3 sequence tables using tra, v1, NCBITAXON:9606 "
	@echo "make create-table-trb			-- generate unique CDR3 sequence tables using trb, v1, NCBITAXON:9606 "
	@echo "make create-table-trd			-- generate unique CDR3 sequence tables using trd, v1, NCBITAXON:9606 "
	@echo "make create-table-trg			-- generate unique CDR3 sequence tables using trg, v1, NCBITAXON:9606 "
	@echo "make create-table-igh			-- generate unique CDR3 sequence tables using igh, v1, NCBITAXON:9606 "
	@echo "make create-table-igk			-- generate unique CDR3 sequence tables using igk, v1, NCBITAXON:9606 "
	@echo "make create-table-igl			-- generate unique CDR3 sequence tables using igl, v1, NCBITAXON:9606 "
	@echo "make create-all-table			-- generate unique CDR3 sequence tables using all locus, v1, NCBITAXON:9606"
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "               -- Extract CDR3 sequences from database and generate distance=1 graph."
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo " (run within docker)"
	@echo "make cdr3-edgelist			-- Generate distance 1 edgelist using $(LOCUS), $(VERSION)"
	@echo "make cdr3-edgelist-threads		-- Generate distance 1 edgelist using $(LOCUS), $(VERSION), $(N_THREADS)"
	@echo "make cdr3-edgelist-all			-- Generate distance 1 edgelist using all locus v1 with 8 threads"
	@echo "make clean				-- Remove binary in bin/ and files in output/"
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "                                -- Graph Analysis commands"
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo " (run within docker)"
	@echo "make build-graph     			-- Builds the graph using $(LOCUS), $(VERSION), $(SPECIES) and save it in binary format"
	@echo "make graph-stat      			-- Generate general statistics using $(LOCUS), $(VERSION) for the graph"
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo "                             -- Graph-Epitope Analysis commands"
	@echo "-------------------------------------------------------------------------------------------------------------"
	@echo ""
	@echo "make get-epitope-info			-- Get epitope information from the database and save in a file using $(LOCUS), $(VERSION)"
	@echo "make get-epitope-info-tra		-- Get epitope information from the database and save in a file using tra, v1 "
	@echo "make get-epitope-info-trb		-- Get epitope information from the database and save in a file using trb, v1 "
	@echo "make get-epitope-info-trd		-- Get epitope information from the database and save in a file using trd, v1 "
	@echo "make get-epitope-info-trg		-- Get epitope information from the database and save in a file using trg, v1 "
	@echo "make get-epitope-info-igh		-- Get epitope information from the database and save in a file using igh, v1 "
	@echo "make get-epitope-info-igk		-- Get epitope information from the database and save in a file using igk, v1 "
	@echo "make get-epitope-info-igl		-- Get epitope information from the database and save in a file using igl, v1 "
	@echo "make get-epitope-info			-- Get epitope information from the database and save in a file."
	@echo "make epitope-analysis 			-- Do epitope analysis using available IEDB epitopes "
	@echo ""
	@echo "-------------------------------------------------------------------------------------------------------------"

docker:
	docker build -t airrknowledge/ak-graph .

create-table:
	python create_junction_tables.py $(LOCUS) --VERSION $(VERSION)
create-table-tra:
	python create_junction_tables.py tra --VERSION v1
create-table-trb:
	python create_junction_tables.py trb --VERSION v1
create-table-trd:
	python create_junction_tables.py trd --VERSION v1
create-table-trg:
	python create_junction_tables.py trg --VERSION v1
create-table-igh:
	python create_junction_tables.py igh --VERSION v1
create-table-igk:
	python create_junction_tables.py igk --VERSION v1
create-table-igl:
	python create_junction_tables.py igl --VERSION v1

## Create all table at once. Skipping bcr for now as they are not available
create-all-table: create-table-tra create-table-trb create-table-trd create-table-trg

#set analysis type to build_graph
build-graph:
	python graph_analysis_using_networkit.py build_graph $(LOCUS) --VERSION $(VERSION) --DATA_DIR $(DATA_DIR)

#set analysis type to graph_stat
graph-stat:
	python graph_analysis_using_networkit.py graph_stat $(LOCUS) --VERSION $(VERSION) --DATA_DIR $(DATA_DIR)

get-epitope-info:
	python get_epitope_information.py $(LOCUS) --VERSION $(VERSION)

get-epitope-info-tra:
	python get_epitope_information.py tra --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-trb:
	python get_epitope_information.py trb --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-trd:
	python get_epitope_information.py trd --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-trg:
	python get_epitope_information.py trg --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-igh:
	python get_epitope_information.py igh --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-igk:
	python get_epitope_information.py igk --DATA_DIR $(DATA_DIR) --VERSION v1
get-epitope-info-igl:
	python get_epitope_information.py igl --DATA_DIR $(DATA_DIR) --VERSION v1

get-all-epitope-info: get-epitope-info-tra get-epitope-info-trb get-epitope-info-trd get-epitope-info-trg

epitope-analysis:
	python CDR3AndEpitopeAnlysis.py --LOCUS trb

# Build compairr library
compairr/libcompairr.a:
	$(MAKE) -C compairr install-lib

# Build stream query program
stream_query: stream_query.cpp compairr/libcompairr.a
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
	mv $@ $(BIN_DIR)/$@
# Run compairr_lib, stream_query and cdr3_dist1_edgelist together
cdr3-dist1-edgelist: stream_query
	$(BIN_DIR)/stream_query $(LOCUS) $(VERSION) $(DATA_DIR)


stream_query_threads: stream_query_threads.cpp compairr/libcompairr.a
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
	mv $@ $(BIN_DIR)/$@

cdr3-edgelist-threads: stream_query_threads
	$(BIN_DIR)/stream_query_threads $(LOCUS) $(VERSION) $(DATA_DIR) $(N_THREADS)

cdr3-edgelist-tra: stream_query_threads
	$(BIN_DIR)/stream_query_threads tra v1 $(DATA_DIR) 8
cdr3-edgelist-trb: stream_query_threads
	$(BIN_DIR)/stream_query_threads trb v1 $(DATA_DIR) 8
cdr3-edgelist-trd: stream_query_threads
	$(BIN_DIR)/stream_query_threads trd v1 $(DATA_DIR) 8
cdr3-edgelist-trg: stream_query_threads
	$(BIN_DIR)/stream_query_threads trg v1 $(DATA_DIR) 8

cdr3-edgelist-all: cdr3-edgelist-trb cdr3-edgelist-trd cdr3-edgelist-trg

# Add -lcompairr only for certain targets
$(COMPAIRR_PROGS): LDLIBS += -lcompairr

# %: %.cpp
# 	@mkdir -p $(BIN_DIR)
# 	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
# 	mv $@ $(BIN_DIR)/$@

# Build all programs with: make all
all: $(PROG_NAMES)

clean:
	rm $(BIN_DIR)/*
	rm $(OUTPUT_DIR)/*
