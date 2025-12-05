#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp> // Required for graph serialization
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <chrono>
#include <ctime>

//#define CHUNK 10000
#define CHUNK 10000000

struct VertexProperties {
  std::string sequence;

  // Serialization function for VertexProperties
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & sequence;
  }
};

struct EdgeProperties {
  int distance;

  // Serialization function for VertexProperties
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & distance;
  }
};

typedef boost::adjacency_list<
    boost::setS,
    boost::vecS,
    boost::undirectedS,
    VertexProperties,
    EdgeProperties
> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;


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

int main() {
  show_time();
  auto start_time = std::chrono::high_resolution_clock::now();
  auto last_time = std::chrono::high_resolution_clock::now();

  Graph g;
  std::map<std::string, Vertex> seq_to_vertex;

  // compare output
  std::ifstream infile("/data/data1/s234499/projects/ak-graph/compairr_edges.tsv");
  std::string line;

  // read nodes
  // Skip header if present
  std::getline(infile, line);

  unsigned long cnt = 0;
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
    if (!distance) continue;

    // Add vertices if not present
    if (seq_to_vertex.find(junction1) == seq_to_vertex.end()) {
      std::string* j1 = new std::string(junction1); // Allocated on heap
      Vertex v = boost::add_vertex(VertexProperties{*j1}, g);
      //g[v].sequence = junction1;
      seq_to_vertex[junction1] = v;
    }
    if (seq_to_vertex.find(junction2) == seq_to_vertex.end()) {
      std::string* j2 = new std::string(junction2); // Allocated on heap
      Vertex v = boost::add_vertex(VertexProperties{*j2}, g);
      //g[v].sequence = junction2;
      seq_to_vertex[junction2] = v;
    }

    // Add edge
    //std::cout << seq_to_vertex[junction1] << ", " << junction1 << std::endl;
    boost::add_edge(seq_to_vertex[junction1], seq_to_vertex[junction2], EdgeProperties{distance}, g);


    if ((cnt % CHUNK) == 0) {
      show_elapsed(start_time, last_time);
      last_time = std::chrono::high_resolution_clock::now();
      std::cout << "lines: " << cnt << std::endl;
      std::cout << "Graph built with " 
		<< boost::num_vertices(g) << " vertices and " 
		<< boost::num_edges(g) << " edges." << std::endl;
    }
  }

  show_elapsed(start_time, last_time);
  std::cout << "lines: " << cnt << std::endl;
  std::cout << "Graph built with " 
	    << boost::num_vertices(g) << " vertices and " 
	    << boost::num_edges(g) << " edges." << std::endl;

  show_time();

  // save graph
  std::cout << "Serializing graph." << std::endl;
  last_time = std::chrono::high_resolution_clock::now();
  std::ofstream ofs("graph.bin", std::ios::binary);
  if (!ofs.is_open()) {
    std::cerr << "Error opening file for writing." << std::endl;
    return 1;
  }
  boost::archive::binary_oarchive oa(ofs);
  oa << g;
  ofs.close();
  std::cout << "Graph serialized to graph.bin" << std::endl;
  show_elapsed(start_time, last_time);

  // load graph
  std::cout << "Loading graph." << std::endl;
  last_time = std::chrono::high_resolution_clock::now();
  Graph loaded_g;
  {
    std::ifstream ifs("graph.bin", std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ia >> loaded_g; // Deserialize the graph
  }
  std::cout << "Graph deserialized from graph.bin" << std::endl;
  show_elapsed(start_time, last_time);

  std::cout << std::endl << "Graph built with " 
	    << boost::num_vertices(loaded_g) << " vertices and " 
	    << boost::num_edges(loaded_g) << " edges." << std::endl;

  return 0;
}
