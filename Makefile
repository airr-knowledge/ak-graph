
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

BIN_DIR = ./bin
OUTPUT_DIR = ./output
COMPAIRR_PROGS = stream_query stream_query_no_output stream_query_line stream_query_threads stream_query_threads_no_output
PROG_NAMES = stream_query stream_query_threads

CXX      := g++
CXXFLAGS := -Wall -O2
LDLIBS   := -lpqxx

help:
	@echo ""
	@echo "----------------------------------------------------------------------------------------"
	@echo "                          AIRR Knowledge graph analysis"
	@echo "----------------------------------------------------------------------------------------"
	@echo "                                 Using DB: $(PG_DISPLAY_CONN)"
	@echo "       Host folder for ak-etvl (/ak_data): $(IMPORT_DATA)"
	@echo "Host folder for ak-graph (/ak_graph_data): $(GRAPH_DATA)"
	@echo "----------------------------------------------------------------------------------------"
	@echo ""
	@echo "make docker                              -- Build docker image"
	@echo ""
	@echo "----------------------------------------------------------------------------------------"
	@echo "  (run within docker)"
	@echo "make all                                 -- Compiles all C++ programs"
	@echo "make stream_query                        -- Query DB, perform overlap, write pairs file"
	@echo "make stream_query_no_output              -- Same as stream_query w/o output"
	@echo "make stream_query_threads                -- Same as stream_query but with threads"
	@echo "make stream_query_threads_no_output      -- Same as stream_query_threads w/o output"
	@echo "make clean                               -- Remove binary in bin/ and files in output/"
	@echo ""
	@echo "----------------------------------------------------------------------------------------"

docker:
	docker build -t airrknowledge/ak-graph .

# Add -lcompairr only for certain targets
$(COMPAIRR_PROGS): LDLIBS += -lcompairr

%: %.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
	mv $@ $(BIN_DIR)/$@

# Build all programs with: make all
all: $(PROG_NAMES)

clean:
	rm $(BIN_DIR)/*
	rm $(OUTPUT_DIR)/*
