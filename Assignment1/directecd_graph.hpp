#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//A large sehection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
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

//Forward declarations for classes below so they can be used below without worrying too much about the ordering.
template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:

    int edges_count; // num of edges in the directed graph
    std::vector<std::vector<int> > adj_matrix; // adjacent matrix of the directed graph


    int get_index(const vertex&) const; // get the index of the corresponding vertex
    bool is_index_valid(const int&, const int&) const;  // check if the edge constructed by the two indexes is valid


public:
    std::vector<vertex> vertices; // vector of vertex, put as public member to make it accessible to class vertex_iterator

    directed_graph(); //A constructor for directed_graph. The graph should start empty.
    ~directed_graph(); //A destructor. Depending on how you do things, this may
    //not be necessary.

    bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.

    bool adjacent(const vertex&, const vertex&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

    void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
    void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

    void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
    void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.

    std::size_t in_degree(const vertex&) const; //Returns number of edges coming in to a vertex.
    std::size_t out_degree(const vertex&) const; //Returns the number of edges leaving a vertex.
    std::size_t degree(const vertex&) const; //Returns the degree of the vertex (both in and out edges).

    std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
    std::size_t num_edges() const; //Returns the total number of edges in the graph.

    std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
    std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

    vertex_iterator<vertex> begin(); //Returns a graph_iterator pointing to the start of the vertex set.
    vertex_iterator<vertex> end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

    neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
    neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

    std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
    std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

    directed_graph<vertex> out_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
    directed_graph<vertex> in_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

    bool reachable(const vertex&, const vertex&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.

};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator {

private:

    directed_graph<vertex> graph; // the corresponding directed graph
    int position; // the iterator position

public:
    vertex_iterator(const vertex_iterator<vertex>&);
    vertex_iterator(const directed_graph<vertex>&, std::size_t);
    ~vertex_iterator();
    vertex_iterator<vertex> operator=(const vertex_iterator<vertex>&);
    bool operator==(const vertex_iterator<vertex>&) const;
    bool operator!=(const vertex_iterator<vertex>&) const;
    vertex_iterator<vertex> operator++();
    vertex_iterator<vertex> operator++(int);
    vertex operator*();
    vertex* operator->();
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator {

private:

    directed_graph<vertex> graph; // the corresponding directed graph
    std::vector<vertex> neighbours; // the vector of all neighbour vertices
    int position; // the position of neighbour iterator
public:
    neighbour_iterator(const neighbour_iterator<vertex>&);
    neighbour_iterator(const directed_graph<vertex>&, const vertex&, std::size_t);
    ~neighbour_iterator();
    neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>&);
    bool operator==(const neighbour_iterator<vertex>&) const;
    bool operator!=(const neighbour_iterator<vertex>&) const;
    neighbour_iterator<vertex> operator++();
    neighbour_iterator<vertex> operator++(int);
    vertex operator*();
    vertex* operator->();
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename vertex> int directed_graph<vertex>::get_index(const vertex& u) const {
    // for each vertex
    for (unsigned i = 0; i < vertices.size(); ++i) {
        // if the vertex equal to the vertex we are looking for
        if (vertices[i] == u)
            // vertex has been found, and return the index of the vertex position
            return i;
    }
    // else return "vertex does not exist" flag
    return -1;
}

template <typename vertex> bool directed_graph<vertex>::is_index_valid(const int& u, const int& v) const {
    // the two indexes must be greater than or equal to zero, and not equal each other.
    return (u >= 0) && (v >= 0) && (u != v);
}

template <typename vertex> directed_graph<vertex>::directed_graph() {
    //initialize the directed with empty
    edges_count = 0;
}

template <typename vertex> directed_graph<vertex>::~directed_graph() {}

template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const {
    // traverse all the vertices to check whether containes vertex u
    for(unsigned i = 0; i < num_vertices(); ++i){
        if(u == vertices[i]) return true; // contains u
    }
    return false; // u is not in the vertices of directed graph
}
template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) const {
    // if vertex u and vertex v are in the directed graph, just check the adjacent matrix element adj_matrix[u_pos][v_pos]
    // if adj_matrix[u_pos][v_pos] > 0, return true
    // 	  else return false
    // if vertex u or vertex v is not in the directed graph, just return false
    return (contains(u) && contains(v)) ? adj_matrix[get_index(u)][get_index(v)] > 0 : false;
}

template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
    //if vertex u is not in the graph, then add it
    //else do nothing
    if(!contains(u)){
        vertices.push_back(u); // add vertex u to vertices
        for(auto &row : adj_matrix){ // add column to each row of the adjacent matrix
            row.push_back(0);
        }
        adj_matrix.push_back(std::vector<int>(vertices.size(), 0)); // add row to the adjacent matrix
    }
}

template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
    int u_index = get_index(u), v_index = get_index(v); // first find the indexes of u and v
    // check if edge (u, v) already exists in the directed graph, if NOT, add it
    // else do nothing
    if(is_index_valid(u_index, v_index) && adj_matrix[u_index][v_index] == 0){
        adj_matrix[u_index][v_index] = 1; // put the corresponding element of adjacent matrix to 1
        ++edges_count; // number of edges increment by 1
    }
}

template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
    // get index of vertex
    int u_index = get_index(u);
    // if index is valid
    if (u_index >= 0) {
        // remove edges and edge weights from edge and weight count variables
        for (unsigned i = 0; i < adj_matrix[u_index].size(); ++i) {
            if (adj_matrix[u_index][i] > 0) --edges_count;
        }
        // remove vertex from vertex list
        vertices.erase(vertices.begin() + u_index);
        // remove vertex row from adj_matrix
        adj_matrix.erase(adj_matrix.begin() + u_index);
        // remove corresponding vertex values from each row
        for (unsigned i = 0; i < adj_matrix.size(); i++) {
            adj_matrix[i].erase(adj_matrix[i].begin() + u_index);
        }
    }
}

// is just on the contrary to add_edge, see add_edge
template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
    int u_index = get_index(u), v_index = get_index(v); // find the indexes of u and v
    // check if edge(u, v) is in the directed graph, if YES then remove it,
    // else do nothing
    if(is_index_valid(u_index, v_index) && adj_matrix[u_index][v_index] > 0){
        adj_matrix[u_index][v_index] = 0; // modify the element of adjacent matrix
        --edges_count; // number of edges decrease by 1
    }
}

template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) const {
    int ans = 0;
    int u_index = get_index(u); // get the index of vertex u
    //u_index is valid;
    if(u_index >= 0){
        // check the u_index row of adjacent matrix count how many 1 s
        for(unsigned i = 0; i < adj_matrix[u_index].size(); ++i){
            ans += (adj_matrix[i][u_index] > 0) ? 1 : 0; // add 1 if adjacent matrix element is 1, else add 0
        }
    }
    return ans;
}

template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u) const {
    int ans = 0;
    int u_index = get_index(u); // get the index of vertex u
    //u_index is valid;
    if(u_index >= 0){
        // observe the corresponding column of the matrix , see how many 1 s
        for(unsigned i = 0; i < adj_matrix.size(); ++i){
            ans += (adj_matrix[u_index][i] > 0) ? 1 : 0; // add 1 if the element equals 1, else add 0
        }
    }
    return ans;
}

template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u) const {
    // degree = in_degree + out_degree, just call the methods implemented above
    return in_degree(u) + out_degree(u);
}

template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const {
    // the size of vertices
    return vertices.size();
}
template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const {
    // the number of edges
    return edges_count;
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() {
    // return the vertices
    return vertices;
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) {
    std::vector<vertex> neighbours;// container to record neighbour vertices
    if(contains(u)){ // if u is in the directed graph
        int u_index = get_index(u); // get the index of u
        // observe the u_index row of adjacent matrix, see how many 1 s
        for(unsigned i = 0; i < adj_matrix.size(); ++i){
            if(adj_matrix[u_index][i] > 0){
                neighbours.push_back(vertices[i]); // put to neighbours vector
            }
        }
    }
    return neighbours;
}

template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() {
    return vertex_iterator<vertex>(*this, 0);
}
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() {
    return vertex_iterator<vertex>(*this, vertices.size());
}

template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) {
    return neighbour_iterator<vertex>(*this, u, 0);
}

template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) {
    return neighbour_iterator<vertex>(*this, u, get_neighbours(u).size());
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) {
    std::vector<bool> visited(vertices.size(), false);
    std::vector<vertex> ans;
    std::stack<vertex> s;

    if(get_index(u) >= 0){
        s.push(u);

        while(!s.empty()){
            int index = get_index(s.top());
            s.pop();

            if(!visited[index]){
                visited[index] = true;
                ans.push_back(vertices[index]);

                // traverse all the neighbours of vertices[index]
                for(unsigned i = vertices.size(); i > 0; --i){
                    if(adj_matrix[index][i - 1] > 0){
                        s.push(vertices[i - 1]);
                    }
                }
            }
        }
    }

    return ans;
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) {
    std::vector<bool> visited(vertices.size(), false);
    std::vector<vertex> ans;
    std::queue<vertex> q;

    if(get_index(u) >= 0){
        q.push(u);

        while(!q.empty()){
            int index = get_index(q.front());
            q.pop();

            if(!visited[index]){
                visited[index] = true;
                ans.push_back(vertices[index]);

                // traverse all the neighbours of vertices[index]
                for(unsigned i = 0; i < vertices.size(); ++i){
                    if(adj_matrix[index][i] > 0){
                        q.push(vertices[i]);
                    }
                }
            }
        }
    }

    return ans;
}


template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) {

    directed_graph<vertex> g;
    std::vector<bool> visited(vertices.size(), false);
    std::queue<std::pair<int, int> > q;

    int u_index = get_index(u);
    if(u_index >= 0){
        q.push({u_index, u_index});
        while (!q.empty()){
            std::pair<int, int> n = q.front();
            q.pop();
            if (!visited[n.first]){
                visited[n.first] = true;
                if(n.first != n.second){
                    g.add_vertex(vertices[n.first]);
                    g.add_vertex(vertices[n.second]);
                    g.add_edge(n.second, n.first);
                }
                for (unsigned i = 0; i < adj_matrix.size(); ++i){
                    if (adj_matrix[n.first][i]){
                        q.push({i, n.first});
                    }
                }
            }
        }
    }


    return g;

}

template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) {

    directed_graph<vertex> g;
    std::vector<bool> visited(vertices.size(), false);
    std::queue<std::pair<int, int> > q;

    int u_index = get_index(u);
    if(u_index >= 0){
        q.push({u_index, u_index});
        while (!q.empty()){
            std::pair<int, int> n = q.front();
            q.pop();
            if (!visited[n.first]){
                visited[n.first] = true;
                if(n.first != n.second){
                    g.add_vertex(vertices[n.first]);
                    g.add_vertex(vertices[n.second]);
                    g.add_edge(n.first, n.second);
                }
                for (unsigned i = 0; i < adj_matrix.size(); ++i){
                    if (adj_matrix[i][n.first]){
                        q.push({i, n.first});
                    }
                }
            }
        }
    }


    return g;

}

template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) const {
    std::vector<bool> visited(vertices.size(), false);
    std::vector<vertex> ans;
    std::queue<vertex> q;

    if(get_index(u) >= 0){
        q.push(u);

        while(!q.empty()){
            vertex t = q.front();

            if(t == v) return true;
            int index = get_index(t);
            q.pop();

            if(!visited[index]){
                visited[index] = true;
                ans.push_back(vertices[index]);

                // traverse all the neighbours of vertices[index]
                for(unsigned i = 0; i < vertices.size(); ++i){
                    if(adj_matrix[index][i] > 0){
                        q.push(vertices[i]);
                    }
                }
            }
        }
    }

    return false;
}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const vertex_iterator<vertex>& other) {
    this->graph = other.graph;
    this->position = other.position;
}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const directed_graph<vertex>& graph, std::size_t position) {
    this->graph = graph;
    this->position = position;
}

template <typename vertex> vertex_iterator<vertex>::~vertex_iterator() {
    position = 0;
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator=(const vertex_iterator<vertex>& other) {

    this->graph = other.graph;
    this->position = other.position;

    return *this;
}

template <typename vertex> bool vertex_iterator<vertex>::operator==(const vertex_iterator<vertex>& other) const {

    return this->poistion == other.position;
}

template <typename vertex> bool vertex_iterator<vertex>::operator!=(const vertex_iterator<vertex>& other) const {

    return !(this->position == other.position);
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++() {

    position++;
    return *this;

}


template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++(int) {

    ++position;
    return *this;

}


template <typename vertex> vertex vertex_iterator<vertex>::operator*() {
    auto v = graph.vertices[position];
    return v;

}
template <typename vertex> vertex* vertex_iterator<vertex>::operator->() {
    auto v = graph.vertices[position];
    return &v;
}



template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const neighbour_iterator<vertex>& other) {
    this->graph = other.graph;
    this->position = other.position;
    this->neighbours = other.neighbours;
}

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position) {
    this->graph = graph;
    this->position = position;
    for(auto v : this->graph.get_neighbours(u)){
        neighbours.push_back(v);
    }
}

template <typename vertex> neighbour_iterator<vertex>::~neighbour_iterator() {
    position = 0;
    neighbours.clear();
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator=(const neighbour_iterator<vertex>& other) {
    this->neighbours = other.neighbours;
    this->position = other.position;
    this->graph = other.graph;
    return *this;
}

template <typename vertex> bool neighbour_iterator<vertex>::operator==(const neighbour_iterator<vertex>& other) const {
    return this->position == other.position;
}
template <typename vertex> bool neighbour_iterator<vertex>::operator!=(const neighbour_iterator<vertex>& other) const {

    return !(this->position == other.position);

}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++() {
    ++position;
    return *this;
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++(int) {
    position++;
    return *this;
}

template <typename vertex> vertex neighbour_iterator<vertex>::operator*() {
    auto v = neighbours[position];
    return v;
}

template <typename vertex> vertex* neighbour_iterator<vertex>::operator->() {
    auto v = neighbours[position];
    return &v;
}


#endif
