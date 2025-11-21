
docker:
	docker build -t airrknowledge/ak-graph .

test: test.cpp
	g++ -o test test.cpp -lpqxx
