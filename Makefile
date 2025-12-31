
docker:
	docker build -t airrknowledge/ak-graph .

test: test.cpp
	g++ -o test test.cpp -lpqxx

# Compile get_repertoire.cpp
get_repertoire: get_repertoire.cpp
	g++ -o get_repertoire get_repertoire.cpp -lpqxx

# Compile run_overlap.cpp
run_overlap: run_overlap.cpp
	g++ -o run_overlap run_overlap.cpp -lpqxx -lcompairr

# Compile stream_query.cpp
stream_query: stream_query.cpp
	g++ -o stream_query stream_query.cpp -lpqxx -lcompairr

# Optional: clean compiled binaries
clean:
	rm -f test get_repertoire run_overlap stream_query output_matrix.tsv output_pairs.tsv repertoire.tsv compairr.log
