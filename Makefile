
docker:
	docker build -t airrknowledge/ak-graph .

test: test.cpp
	g++ -o test test.cpp -lpqxx

# Compile get_repertoire.cpp
get_repertoire: get_repertoire.cpp
	g++ -o get_repertoire get_repertoire.cpp -lpqxx

# Optional: clean compiled binaries
clean:
	rm -f test get_repertoire
