
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
	@echo "make test                  -- read code"
	@echo "make graph_analysis        -- Read compairr output, generate and save graph"
	@echo "make connected_components  -- Compute connected components"
	@echo "make thread_graph          -- thread example"
	@echo ""
	@echo "------------------------------------------------------------"

docker:
	docker build -t airrknowledge/ak-graph .

test: test.cpp
	g++ -o test test.cpp -lpqxx

graph_analysis: graph_analysis.cpp
	g++ -std=c++20 -o graph_analysis graph_analysis.cpp -lpqxx -lboost_serialization -lboost_system

connected_components: connected_components.cpp
	g++ -std=c++20 -o connected_components connected_components.cpp -lpqxx -lboost_serialization -lboost_system

thread_graph: thread_graph.cpp
	g++ -std=c++20 -o thread_graph thread_graph.cpp -lpqxx -lboost_serialization -lboost_thread -lboost_system

# Compile get_repertoire.cpp
get_repertoire: get_repertoire.cpp
	g++ -o get_repertoire get_repertoire.cpp -lpqxx

# Optional: clean compiled binaries
clean:
	rm -f test get_repertoire
