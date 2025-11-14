#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp> // For print_graph
#include <pqxx/pqxx>

/// Query employees from database.  Return result.
pqxx::result query()
{
  pqxx::connection cx{"postgresql://postgres:example@ak-db.airr-knowledge.org/airrkb_v1"};
  pqxx::work tx{cx};

  // Silly example: Add up all salaries.  Normally you'd let the database do
  // this for you.
//   long total = 0;
//   for (auto [salary] : tx.query("SELECT salary FROM Employee"))
//     total += salary;
//   std::cout << "Total salary: " << total << '\n';

  // Execute and process some data.
  pqxx::result r{tx.exec("SELECT akc_id, trb_chain FROM TCellReceptor")};
  for (auto row: r) {
    std::cout
      // Address column by name.  Use c_str() to get C-style string.
      << row["akc_id"].c_str()
      << " makes "
      // Address column by zero-based index.  Use as<int>() to parse as int.
      << row["tra_chain"].c_str()
      << "."
      << std::endl;
      
      break;
    }

  // Not really needed, since we made no changes, but good habit to be
  // explicit about when the transaction is done.
  tx.commit();

  // Connection object goes out of scope here.  It closes automatically.
  // But the result object remains valid.
  return r;
}

int main() {
    // Define the graph type:
    // adjacency_list<EdgeList, VertexList, Directedness, VertexProperties, EdgeProperties>
    // In this example:
    // - listS for EdgeList (fast edge insertion/deletion)
    // - vecS for VertexList (fast vertex access by index)
    // - directedS for a directed graph
    // - no_property for vertex and edge properties (no custom properties)
    using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS>;

    // Create a graph object
    Graph g;

    // Add vertices
    boost::add_vertex(g); // Vertex 0
    boost::add_vertex(g); // Vertex 1
    boost::add_vertex(g); // Vertex 2
    boost::add_vertex(g); // Vertex 3

    // Add edges
    // add_edge(source, target, graph_object)
    boost::add_edge(0, 1, g);
    boost::add_edge(0, 2, g);
    boost::add_edge(1, 3, g);
    boost::add_edge(2, 3, g);

    // Print the graph structure
    std::cout << "Graph structure:" << std::endl;
    boost::print_graph(g, std::cout);

    // Iterate over vertices and print their indices
    std::cout << "\nVertices:" << std::endl;
    boost::graph_traits<Graph>::vertex_iterator vi, vend;
    for (boost::tie(vi, vend) = boost::vertices(g); vi != vend; ++vi) {
        std::cout << *vi << " ";
    }
    std::cout << std::endl;

    // Iterate over edges and print source and target
    std::cout << "\nEdges:" << std::endl;
    boost::graph_traits<Graph>::edge_iterator ei, eend;
    for (boost::tie(ei, eend) = boost::edges(g); ei != eend; ++ei) {
        std::cout << "(" << boost::source(*ei, g) << ", " << boost::target(*ei, g) << ") ";
    }
    std::cout << std::endl;

    query();
    
    return 0;
}
