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
#include <boost/config.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/connected_components.hpp>

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

  last_time = std::chrono::high_resolution_clock::now();
  using namespace boost;
  {
    std::vector< int > component(num_vertices(loaded_g));
    auto num = connected_components(loaded_g, &component[0]);

    std::vector< int >::size_type i;
    std::cout << "Total number of components: " << num << std::endl;
    //for (i = 0; i != component.size(); ++i)
    //  cout << "Vertex " << i << " is in component " << component[i]
    //   << endl;
    //cout << endl;
  }
  show_elapsed(start_time, last_time);

  return 0;
}
