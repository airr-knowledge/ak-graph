BIN_DIR = ./bin
COMPAIRR_PROGS = run_overlap stream_query stream_query_no_output stream_query_line stream_query_threads stream_query_threads_no_output
PROG_NAMES = $(COMPAIRR_PROGS) test get_repertoire
FILE_NAMES = output_matrix.tsv output_pairs.tsv repertoire.tsv compairr.log

CXX      := g++
CXXFLAGS := -Wall -O2
LDLIBS   := -lpqxx

# Add -lcompairr only for certain targets
$(COMPAIRR_PROGS): LDLIBS += -lcompairr

%: %.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)
	mv $@ $(BIN_DIR)/$@

# Build all programs with: make all
all: $(PROG_NAMES)

docker:
	docker build -t airrknowledge/ak-graph .

clean:
	rm -f $(addprefix $(BIN_DIR)/,$(PROG_NAMES))
	rm -f $(FILE_NAMES)
