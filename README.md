# ak-graph
AK Graph Algorithms

# Compile and setup
- `git clone https://github.com/airr-knowledge/ak-graph.git`
- CD into ak-graph and init submodules: `git submodule update --init --recursive`
- cp .env.defaults .env
- nano .env # set database and other deployment config
- Run `make` to show help.
- Most make targets are performed within a docker container.
- Run `make docker` to build docker image. The docker will compile all programs.

# docker info
- Compiled programs will be in /ak-graph directory inside docker, but use /work for development.
- /ak-graph/bin and /work/bin are added to PATH

# development, start shell inside docker, graph data at /ak_graph_data is stored in work
- Run `docker run --network ak-db-network -v $PWD:/work -v ./graph_data:/ak_graph_data -it airrknowledge/ak-graph bash` to start `ak-graph` container.
- Create a bash alias to help `alias docker-ak-graph-dev="docker run --network ak-db-network -v $PWD:/work ./graph_data:/ak_graph_data -it airrknowledge/ak-graph bash"`

# production runs
- GRAPH_DATA environment variables is the host folder for the ak-graph data
- Run `docker run --network ak-db-network -v ${GRAPH_DATA}:/ak_graph_data -it airrknowledge/ak-graph bash` for shell in container.

- Run `make create-all-table` to generate unique CDR3 sequence tables in the database.
- Run `make cdr3-edgelist-threads` to extract CDR3 sequences from database and generate distance=1 graph.
- Run `make build-graph` to generate the graph from the edge list and save it in networkit format.
- Run `make graph-stat` to perform some general statitical analysis on graph.
- Run `make epitope-analysis` to perform epitope analysis using available IEDB epitopes in conjunction with the junction_aa graph.

# Timing Analysis
| Program | Threads | CHUNK | real | user | sys | query_and_stream() | run_compairr() | db_read() | compute_overlap() |
|---------|---------|-------|------|------|-----|--------------------|----------------|-----------|-------------------|
| stream_query | 1 | 1000000 | 117m53.018s | 110m32.805s | 7m14.195s | 0m38.7467s | 117m14.116s | 0m29.2801s | 116m43.7309s |
| stream_query_threads | 16 | 1000000 | 81m45.693s | 128m23.770s | 7m56.293s | 0m39.0097s | 81m5.9671s | 0m31.7337s | 80m34.2332s |
| stream_query_threads | 48 | 1000000 | 91m28.135s | 156m2.117s | 17m14.116s | 0m39.2061s | 90m48.3653s | 0m29.2814s | 90m19.0837s |
| stream_query_threads | 64 | 1000000 | 94m53.812s | 182m44.747s | 16m26.075s | 0m38.8788s | 94m14.318 | 0m29.0449s | 93m14.318s |
| stream_query_threads | 96 | 1000000 | 93m55.463s | 241m57.294s | 13m48.724s | 0m39.0387s | 93m15.6882s | 0m28.963s |  92m46.7249s |
| No Output |   |  |  |  |  |  |  |  |  |
| stream_query_no_output | 1 | 1000000 | 51m18.046s | 50m21.783s | 0m51.085s | 0m41.3267s | 50m36.1499s | 0m29.4679s | 50m6.68193s |
| stream_query_threads_no_output | 1 | 1000000 | 51m25.474s | 50m37.785s | 0m42.631s | 0m38.9758s | 50m46.0623s | 0m29.0679s | 50m16.9943s |
| strean_query_threads_no_output | 2 | 1000000 | 32m5.100s | 59m21.272s | 3m5.270s | 0m39.1839s |31m25.4246s | 0m29.4821s | 30m55.9424s |
| stream_query_threads_no_output | 4 | 1000000 | 15m5.802s | 53m24.392s | 1m33.136s | 0m42.0846s | 14m23.194s | 0m30.0385s | 13m53.1553s |
| stream_query_threads_no_output | 8 | 1000000 | 9m37.324s | 63m19.513s | 1m16.326s | 0m39.4549s | 8m57.3962s | 0m29.5826s | 8m27.8135s |
| stream_query_threads_no_output | 16 | 1000000 | 6m28.831s | 75m51.107s | 1m9.708s | 0m39.3275s | 5m48.8363 | 0m29.7976s | 5m19.0386s |
| stream_query_threads_no_output | 32 | 1000000 | 4m42.837 | 97m36.692s | 1m30.149s | 0m39.2919s | 4m3.04755s | 0m29.27s | 3m33.7774s | 
| stream_query_threads_no_output | 48 | 1000000 | 4m25.823s | 124m0.592s | 1m38.509s | 0m38.9672s | 3m46.2524s | 0m29.8667s | 3m16.3855s |
| stream_query_threads_no_output | 64 | 1000000 | 4m18.874s | 160m34.867 | 1m32.924s | 0m39.1994s | 3m39.1776s | 0m29.6685s | 3m9.50895s |
| stream_query_threads_no_output | 96 | 1000000 | 4m3.507s | 207m18.392s | 1m55.971s | 0m39.5649s | 3m23.421s | 0m29.4086s | 2m54.0122s |
