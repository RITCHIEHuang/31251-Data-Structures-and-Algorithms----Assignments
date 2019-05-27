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


/**
* using dfs to help solve the first problem.
* define vector visited to notify the state of vertex: 0 -> not visited, 1 -> visiting, -1 -> visited
*/
template <typename vertex>
void dfs(std::unordered_map<vertex, int> & visited, const vertex & v, const directed_graph<vertex> & d, bool & is_dag){
	visited[v] = 1;
	for(typename directed_graph<vertex>::const_neighbour_iterator it = d.nbegin(v); it != d.nend(v); ++it){
		if(visited[*it] == 1){ 	// if the neighbour vertex of v is in state 1, the graph is acyclic
			is_dag = false;
			return;
		}else if(visited[*it] == -1){ // state -1, can't decide
      continue;
    }else{
      dfs(visited, *it, d, is_dag); // state 0, dfs the neighbour vertex
    }
	}
  visited[v] = -1; // mark vertex v as visited
}
/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  bool res = true;
  // using dfs to check whether the graph is DAG

  std::unordered_map<vertex, int> visited;

  // initialize all vertices as not visited
  for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
  	visited[*it] = 0; 
  }

  // traverse all vertices to decide whether the graph is acyclic
  for(typename directed_graph<vertex>::const_vertex_iterator it = d.begin(); it != d.end(); ++it){
    if(visited[*it] == 0) {
      dfs(visited, *it, d, res);
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
  std::list<vertex> ans; // ans stores the solution
  int num = d.num_vertices();
  // the algorithm will probably change the graph, so I copy it to another graph
  // the basic strategy
  directed_graph<vertex> graph(d); 
	std::queue<vertex> q;
	// initialize queue with vertex whose in degree equal to 0
	for(auto it = d.begin(); it != d.end(); ++it){
		if(d.in_degree(*it) == 0){
			q.push(*it);
		}
	}
  
	while(!q.empty()){
		auto v = q.front();
		ans.push_back(v);
		q.pop();

		for(auto it = d.nbegin(v); it != d.nend(v); ++it){
			graph.remove_edge(v, *it); // remove the edge from v to its' neighbour *it
			if(graph.in_degree(*it) == 0){ // if the neighbour's in degree equal to 0, push it into the queue
				q.push(*it);
			}
		}
	}
	// only DAG can be topological
	return ans.size() == num ? ans : std::list<vertex>();
}


/**
* using backtracking method to solve this problem.record the paths and decide whether it's hamiltonian.
*
*/
template <typename vertex>
bool hamiltonian_helper(const vertex & start, std::unordered_set<vertex> & paths, const directed_graph<vertex> & d){
	if(paths.size() == d.num_vertices()) return true; // find valid paths. return true

	std::unordered_set<vertex> old_paths(paths); // old_paths stores the paths before recursion
	for(auto it = d.nbegin(start); it != d.nend(start); ++it){
		if(paths.count(*it) == 0) { // neighbour is not visited
			paths.insert(*it); // add neighbour vertex to paths
			if(hamiltonian_helper(*it, paths, d)) return true; // dfs 
			paths = old_paths; // backtracking
		}
	}
	return false; // not valid hamiltonian DAG
}
/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
  if(d.num_vertices() <= 0) return true; // empty graph is hamiltonian DAG
	for(auto it = d.begin(); it != d.end(); ++it){ // pick the start vertex
  	std::unordered_set<vertex> paths;
    paths.insert(*it);
		if(hamiltonian_helper(*it, paths, d)) return true;
	}
	return false;
}

/**
* dfs helper to solve the components problem.From start vertex s, dfs the component and store to temp.	
*/
template <typename vertex>
void components_helper(const vertex & s, std::vector<vertex> & temp, std::unordered_map<vertex, bool> & visited, const directed_graph<vertex> & d){
	visited[s] = true; // mark vertex s as visited
  temp.push_back(s); // add s to component s 
	for(auto it = d.nbegin(s); it != d.nend(s); ++it){
		if(!visited[*it]){
			components_helper(*it, temp, visited, d);
		}
	}
}

template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
	// transform the directed_graph to undirected_graph, so we can solve the problem.
	directed_graph<vertex> graph(d);
	for(auto it = d.begin(); it != d.end(); ++it){
		vertex v = *it;
		for(auto nit = d.nbegin(v); nit != d.nend(v); ++nit){
			if(d.adjacent(v, *nit) && !d.adjacent(*nit, v)){
				graph.add_edge(*nit, v);
			}
		}
	}

	std::vector<std::vector<vertex>> res;
	std::unordered_map<vertex, bool> visited;
	for(auto it = d.begin(); it != d.end(); ++it){
		visited[*it] = false;
	}
	for(auto it = d.begin(); it != d.end(); ++it){
    if(!visited[*it]){
      std::vector<vertex> temp;
		  components_helper(*it, temp, visited, graph);
		  res.push_back(temp);
		  temp.clear();
    }
	}

  return res;
}


template <typename vertex>
void ssc_helper(const vertex & v, int &index, 
	std::unordered_map<vertex, int> & indices, 
	std::unordered_map<vertex, int> & lows, 
	std::unordered_map<vertex, bool> & on_stack, 
	std::stack<vertex> & s,
  std::vector<std::vector<vertex>> & components,
	const directed_graph<vertex> & d){

	indices[v] = index;
	lows[v] = index;
	++index;
	s.push(v);
	on_stack[v] = true;
  vertex w;
	
	for(auto it = d.nbegin(v); it != d.nend(v); ++it){
		if(indices.count(*it) == 0){
			ssc_helper(*it, index, indices, lows, on_stack, s, components, d);
			lows[v] = std::min(lows[v], lows[*it]);
		}else if(on_stack[*it]){
			lows[v] = std::min(lows[v], indices[*it]);
		}
	}

	if(lows[v] == indices[v]){
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
	int index = 0;
	std::stack<vertex> s;
	std::unordered_map<vertex, int> indices;
	std::unordered_map<vertex, int> lows;
	std::unordered_map<vertex, bool> on_stack;

	std::vector<std::vector<vertex>> res;


	for(auto it = d.begin(); it != d.end(); ++it){
		if(indices.count(*it) == 0){
			ssc_helper(*it, index, indices, lows, on_stack, s, res, d);
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
  size_t vertices_num = d.num_vertices();
  for(auto it = d.begin(); it != d.end(); ++it){
    res[*it] = vertices_num + 1;
  }
  res[u] = 0;
  
  std::queue<vertex> q;
  q.push(u);
  
  while(!q.empty()){
    auto v = q.front();
    q.pop();
    
    for(auto it = d.nbegin(v); it != d.nend(v); ++it){
      if(res[*it] == vertices_num + 1){
        res[*it] = res[v] + 1;
        q.push(*it);
      } 
    }
  }
  return res;
}





