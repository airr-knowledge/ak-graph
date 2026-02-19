BIN_DIR = ./bin
COMPAIRR_PROGS = run_overlap stream_query stream_query_no_output stream_query_line stream_query_threads stream_query_threads_no_output
PROG_NAMES = $(COMPAIRR_PROGS) test get_repertoire
FILE_NAMES = output_matrix.tsv output_pairs.tsv repertoire.tsv compairr.log

CXX      := g++
CXXFLAGS := -Wall -O2
LDLIBS   := -lpqxx

help:
	@echo ""
	@echo "------------------------------------------------------------"
	@echo "AIRR Knowledge graph analysis"
	@echo "------------------------------------------------------------"
	@echo ""
	@echo "make docker             -- Build docker image"
	@echo ""
	@echo "------------------------------------------------------------"
	@echo "  (run within docker)"
	@echo "make all                   -- compiles all C++ programs"
	@echo "make graph_analysis        -- Read compairr output, generate and save graph"
	@echo "make connected_components  -- Compute connected components"
	@echo "make thread_graph          -- thread example"
	@echo ""
	@echo "------------------------------------------------------------"

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
	rm -f $(addprefix $(BIN_DIR)/,$(PROG_NAMES))
	rm -f $(FILE_NAMES)
