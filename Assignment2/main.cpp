/*
* A main function for you to build and run your
* own tests with.
* This file is not part of the marking, so you
* can do anything you want here.
*/
#include <iostream>

#include "directed_graph_algorithms.cpp"

int main() {
    directed_graph<int> g;
    g.add_vertex(1);
    g.add_vertex(2);
    g.add_vertex(3);
    g.add_edge(1, 2);
    g.add_edge(1, 3);

    for(auto it = g.nbegin(1); it != g.nend(1); ++it){
        std::cout << "neighbour " << *it << std::endl;
    }
    std::list<int> msorted = topological_sort(g);
    // for(int node : msorted){
    // 	std::cout << node << " ";
    // }
    // std::cout << std::endl;
    // directed_graph<int> graph(g);
// 	std::cout << graph.num_vertices() << std::endl;
// 	std::cout << "=====" << std::endl;
// 	std::cout << graph.in_degree(2) << std::endl;
// 	std::cout << "=====" << std::endl;
// 	std::cout << graph.in_degree(1) << std::endl;

// 	graph.remove_edge(1, 2);
// 	std::cout << graph.in_degree(2) << std::endl;
// 	std::cout << "=====" << std::endl;
// 	std::cout << g.in_degree(2) << std::endl;


}
