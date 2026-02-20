BIN_DIR = ./bin
OUTPUT_DIR = ./output
COMPAIRR_PROGS = run_overlap stream_query stream_query_no_output stream_query_line stream_query_threads stream_query_threads_no_output
PROG_NAMES = $(COMPAIRR_PROGS) test get_repertoire

CXX      := g++
CXXFLAGS := -Wall -O2
LDLIBS   := -lpqxx

help:
	@echo ""
	@echo "----------------------------------------------------------------------------------------"
	@echo "                          AIRR Knowledge graph analysis"
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
