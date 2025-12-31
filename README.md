# ak-graph
AK Graph Algorithms

# Usage
- Install `compairr` from [GitHub](https://github.com/airr-knowledge/compairr.git)
- Run `docker run --network ak-db-network -v $PWD:/work -it airrknowledge/ak-graph bash` to start `ak-graph` container.
- CD into `compairr` and switch to `akc-enhancements` branch, then run `make install-lib` to install the compairr library.
- CD into `ak-graph` and run `make stream_query.cpp` and finally `./stream_query` to gerneate `output_pairs.tsv`.