#include <boost/graph/adjacency_list.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp> // Required for graph serialization
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <chrono>
#include <ctime>

//#define CHUNK 10000
#define CHUNK 10000000

struct VertexProperty {
    int value;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty> Graph;

boost::mutex graph_mutex;
boost::mutex map_mutex;
Graph g;

void show_time() {
  // Get the current time as a time_point
  auto now = std::chrono::system_clock::now();
 
  // Convert to time_t for human-readable format
  std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
 
  // Print the stored time
  std::cout << "Current time: " << std::ctime(&currentTime);
}

void show_elapsed(auto start_time, auto last_time) {
  // Get the current time as a time_point
  auto end_time = std::chrono::high_resolution_clock::now();

  auto duration = end_time - last_time;
  double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
  std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << std::endl;

  duration = end_time - start_time;
  elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
  std::cout << "Elapsed time from start: " << elapsed_seconds << " seconds" << std::endl;
}

void add_vertex_and_edge(int val1, int val2) {
    boost::mutex::scoped_lock lock(graph_mutex); // Acquire lock before modifying graph
    Graph::vertex_descriptor u = boost::add_vertex(VertexProperty{val1}, g);
    Graph::vertex_descriptor v = boost::add_vertex(VertexProperty{val2}, g);
    boost::add_edge(u, v, g);
    //std::cout << "Added vertices " << val1 << " and " << val2 << " and an edge between them." << std::endl;
}

void generate_graph(int size, int modulus, int period) {
  show_time();
  auto start_time = std::chrono::high_resolution_clock::now();
  auto last_time = std::chrono::high_resolution_clock::now();

  // compare output
  std::ifstream infile("/data/data1/s234499/projects/ak-graph/compairr_edges.tsv");
  std::string line;

  // read nodes
  // Skip header if present
  std::getline(infile, line);

  unsigned long cnt = 0;
  unsigned long insert_cnt = 0;
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string rep1, seq_id1, dup1, v_call1, j_call1, junction1;
    std::string rep2, seq_id2, dup2, v_call2, j_call2, junction2;
    int distance;

    std::getline(ss, rep1, '\t');
    std::getline(ss, seq_id1, '\t');
    std::getline(ss, dup1, '\t');
    std::getline(ss, v_call1, '\t');
    std::getline(ss, j_call1, '\t');
    std::getline(ss, junction1, '\t');

    std::getline(ss, rep2, '\t');
    std::getline(ss, seq_id2, '\t');
    std::getline(ss, dup2, '\t');
    std::getline(ss, v_call2, '\t');
    std::getline(ss, j_call2, '\t');
    std::getline(ss, junction2, '\t');

    ss >> distance;
    ++cnt;
    if (distance == 0) continue;

    if ((cnt % period) == modulus) {
      ++insert_cnt;
      add_vertex_and_edge(10, distance);
    }

    if (cnt % CHUNK == 0) {
      std::cout << "chunk: " << CHUNK << ", cnt: " << cnt << ", insert: " << insert_cnt << std::endl;
      std::cout << "Graph contains " << boost::num_vertices(g) << " vertices and " << boost::num_edges(g) << " edges." << std::endl;
    }
    if (cnt == 4 * CHUNK) break;
    //std::cout << cnt << std::endl;
  }
}

int main() {
  int size = 100;
  int period = 16;
  int num_threads = period;
  std::vector<boost::thread> threads;

  for (int m = 0; m < num_threads; ++m) {
    std::cout << "start thread: " << m << std::endl;
    threads.push_back(boost::thread(generate_graph, size, m, period));
  }

  // Wait for all threads to complete
  for (auto& t : threads) {
    t.join();
  }

  std::cout << "All threads finished." << std::endl;

  std::cout << "Graph contains " << boost::num_vertices(g) << " vertices and " << boost::num_edges(g) << " edges." << std::endl;

  return 0;
}
