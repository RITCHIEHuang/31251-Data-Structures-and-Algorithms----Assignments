/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

/*
* dfs helper for DAG judge
* define the state of vertex: 1 -> visiting,  -1 -> visited
*/
template <typename vertex>
void dag_helper(std::unordered_map<vertex, int> & visited, const vertex & v, const directed_graph<vertex> & d, bool & is_dag){
	// mark vertex v as visiting
	visited[v] = 1; 
	for(typename directed_graph<vertex>::const_neighbour_iterator it = d.nbegin(v); it != d.nend(v); ++it){
		if(visited[*it] == 1){ 	// if the neighbour vertex of v is in state 1, the graph is acyclic
			is_dag = false;
			return;
		}else if(visited[*it] == -1){ // state -1, can't decide
      continue;
    }else{
      dag_helper(visited, *it, d, is_dag); // state 0, dfs the neighbour vertex
    }
	}
  visited[v] = -1; // mark vertex v as visited
}
/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to idfs t
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  bool res = true;

  std::unordered_map<vertex, int> visited;

  // dfs all vertices to decide whether the graph is DAG
  for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
    if(visited.count(*it) == 0) {
      dag_helper(visited, *it, d, res);
  	  if(!res) return false;
    }
  }	
  return res;
}
/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
  std::list<vertex> s; 
  // copy the original graph
  directed_graph<vertex> copy_graph(d); 
	std::queue<vertex> q;
	// initialize queue with vertex whose in degree equal to 0
	for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
		if(d.in_degree(*it) == 0){
			q.push(*it);
		}
	}
  
	while(!q.empty()){
		vertex v = q.front();
		q.pop();
		s.push_back(v);

		for(typename directed_graph<vertex>::const_neighbour_iterator nit = d.nbegin(v); nit != d.nend(v); ++nit){
			vertex w = *nit;
      // remove the edge from v to w
      copy_graph.remove_edge(v, w); 
      // if the in_degree of w is 0, push w into queue
			if(copy_graph.in_degree(w) == 0){ 
				q.push(w);
			}
		}
	}
	return s;
}


/**
* backtracking method based on DFS
*
*/
template <typename vertex>
bool hamiltonian_helper(const vertex & v, std::unordered_set<vertex> & paths, const directed_graph<vertex> & d){
	// all vertices visited, the graph is hamiltonian
  if(paths.size() == d.num_vertices()) return true; 

// old_paths stores the paths before recursion
	std::unordered_set<vertex> old_paths(paths); 
	for(typename directed_graph<vertex>::const_neighbour_iterator nit = d.nbegin(v); nit != d.nend(v); ++nit){
		// neighbour is not visited
    if(paths.count(*nit) == 0) { 
      // add neighbour vertex to paths
			paths.insert(*nit); 
			if(hamiltonian_helper(*nit, paths, d)) return true; 
			paths = old_paths; // backtracking
		}
    // do nothing if neighbour has been visited
	}
	return false; 
}


/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
  // deal with boundary 
  if(d.num_vertices() <= 1) return true; 
  // DFS from each vertex, so that we can deal with all cases
	for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){ 
     	std::unordered_set<vertex> paths;
      paths.insert(*it);
		if(hamiltonian_helper(*it, paths, d)) return true;
	}
	return false;
}

/**
*
* for root vertex u, DFS the component
*/
template <typename vertex>
void components_helper(const vertex & u, std::vector<vertex> & temp, std::unordered_map<vertex, bool> & visited, const directed_graph<vertex> & d){
	// mark vertex s as visited
  visited[u] = true;
  temp.push_back(u);
  // add u to component
	for(typename directed_graph<vertex>::const_neighbour_iterator nit = d.nbegin(u); nit != d.nend(u); ++nit){
		if(!visited[*nit]){
			components_helper(*nit, temp, visited, d);
		}
	}
}

template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
	directed_graph<vertex> graph(d);
	for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
		vertex v = *it;
		for(typename directed_graph<vertex>::const_neighbour_iterator nit = d.nbegin(v); nit != d.nend(v); ++nit){
			vertex w = *nit;
      if(d.adjacent(v, w) && !d.adjacent(w, v)){ 
				graph.add_edge(w, v); 
			}
		}
	}

	std::vector<std::vector<vertex>> res; 
	std::unordered_map<vertex, bool> visited; 
	
  // traverse all vertices in graph
	for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
    	// dfs the component whose root is vertex *it
      if(visited.count(*it) == 0){ 
	     	std::vector<vertex> temp; 
		  	components_helper(*it, temp, visited, graph);
		  	res.push_back(temp);
	    }
	}

  return res;
}

/* 
* implementation of Targan Algorithm
*/
template <typename vertex>
void ssc_helper(const vertex & v, int &current_index, 
	std::unordered_map<vertex, int> & index, 
	std::unordered_map<vertex, int> & low, 
	std::unordered_map<vertex, bool> & on_stack, 
	std::stack<vertex> & s,
  std::vector<std::vector<vertex>> & components,
	const directed_graph<vertex> & d){

	index[v] = low[v] = ++current_index;
	s.push(v);
	on_stack[v] = true;
  vertex w;
	
	for(typename directed_graph<vertex>::const_neighbour_iterator it = d.nbegin(v); it != d.nend(v); ++it){
    vertex u = *it;
    // if index of vertex u is not defined
		if(index.count(u) == 0){
			ssc_helper(u, current_index, index, low, on_stack, s, components, d);
      
      // take the minimum of low[v] and low[u]
			low[v] = std::min(low[v], low[u]);
		}else if(on_stack[u]){ // if vertex u is in the stack, take the minimum of low[v] and index[u]
			low[v] = std::min(low[v], index[u]);
		}
	}

  // deal with the component
	if(low[v] == index[v]){ 
		std::vector<vertex> component;
		do{
			w = s.top();
			s.pop();
			on_stack[w] = false;
			component.push_back(w);
		}while(w != v);
    components.push_back(component);
	}
} 
/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
	
  std::vector<std::vector<vertex>> res;
	std::stack<vertex> s;
  std::unordered_map<vertex, bool> on_stack;
	std::unordered_map<vertex, int> index_map;
	std::unordered_map<vertex, int> low_map;

  int index = 0;
	for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
		// if index of vertex *it is not defined, dfs the SSC
    if(index_map.count(*it) == 0){
			ssc_helper(*it, index, index_map, low_map, on_stack, s, res, d);
		}
	}
  return res;
}


/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  
  std::unordered_map<vertex, std::size_t> res; 
  std::unordered_map<vertex, bool> visited;
  for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
    // initialize the distance to u is 0 and other distances is the number of vertices in d plus 1
    if(*it == u) res[*it] = 0;
    else res[*it] = d.num_vertices() + 1; 
  }
  
  std::queue<vertex> q;
  q.push(u);
  visited[u] = true;
  
  // BFS algorithm
  while(!q.empty()){
    vertex v = q.front();
    q.pop();
    
  
    for(typename directed_graph<vertex>::const_neighbour_iterator it = d.nbegin(v); it != d.nend(v); ++it){
      if(visited.count(*it) == 0){
        visited[*it] = true;
        res[*it] = res[v] + 1;
        q.push(*it);
      }
    }
  }
  return res;
}






