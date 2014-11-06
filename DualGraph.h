#pragma once
#include <vector>
#include <map>

class GraphEdge;

class GraphNode {
public:
    GraphNode(int id, double weight);
    ~GraphNode();
    int  id() const { return m_id; }
    void addEdge(GraphEdge* e) { m_edges.push_back(e); }
	double & weight() { return m_weight; }
	double weight() const { return m_weight; }

    //print
	void print();

    //iterator
	struct NodeIter {
        NodeIter(GraphNode *node) : m_node(node), m_it(node->m_edges.begin()), m_end(node->m_edges.end()) {}
        void operator++() { ++m_it; }
        bool end() { return m_it == m_end; }
		GraphNode* operator*();// { return (*m_it)->other(m_node); }
        std::vector<GraphEdge*>::iterator m_it, m_end;
		GraphNode *m_node;
	};

    struct EdgeIter {
        EdgeIter(GraphNode *node) : m_it(node->m_edges.begin()), m_end(node->m_edges.end()) {}
        void operator++() { ++m_it; }
        bool end() { return m_it == m_end; }
        GraphEdge* operator*() { return *m_it; }
        std::vector<GraphEdge*>::iterator m_it, m_end;
	};
    
private:
    int m_id;
    double m_weight;
    std::vector<GraphEdge*> m_edges;
};

inline std::ostream& operator<< (std::ostream& os, const GraphNode* node) {
	if (node) {
		os << "v " << node->id() << " ";
	}
	else {
		os << "v (null) ";
	}
	return os;
}

class GraphEdge {
public:
    GraphEdge(GraphNode *n1, GraphNode *n2, double weight);
    ~GraphEdge();
	GraphNode *& from() { return m_from;}
	GraphNode *& to() { return m_to;}
	double & weight() { return m_weight; }
	double weight() const { return m_weight; }
	GraphNode * other(GraphNode *n) { 
		if (m_from == n) {
			return m_to;
		}
		else if (m_to == n) {
			return m_from;
		}
		return nullptr;
		//return p = (m_from == n ? m_to : m_from)  ? p : nullptr; 
	}
	void print();
private:
    GraphNode *m_from, *m_to;
	double m_weight;
};

class DualGraph
{
public:
	DualGraph(void);
	~DualGraph(void);

    //node id start from 1
    GraphNode *getNode(int id) { return m_nodes[id - 1]; }
    GraphEdge *getEdge(int nid0, int nid1); 

    GraphNode *addNode(double weight = 1.0);
    GraphEdge *addEdge(GraphNode *from, GraphNode *to, double weight = 1.0);

    int nodeSize() const { return m_nodes.size(); }
    //print
	void print();

    //in and out

    void saveGraph(const char *filename);
	void saveGV(const char *filename);

    void loadGraph(const char *filename);
	void loadMesh(const char *filename);
	void loadTetMesh(const char *filename);

	void saveMetis(const char *filename);

    //iterators
    struct NodeIter {
        NodeIter(DualGraph *graph) : m_it(graph->m_nodes.begin()), m_end(graph->m_nodes.end()) {}
        void operator++() { ++m_it; }
        bool end() { return m_it == m_end; }
        GraphNode* operator*() { return *m_it; }
        std::vector<GraphNode*>::iterator m_it, m_end;
	};

    struct EdgeIter {
        EdgeIter(DualGraph *graph) : m_it(graph->m_edges.begin()), m_end(graph->m_edges.end()) {}
        void operator++() { ++m_it; }
        bool end() { return m_it == m_end; }
        GraphEdge* operator*() { return (*m_it).second; }
        std::map<std::pair<int, int>, GraphEdge*>::iterator m_it, m_end;
	};




private:
	void reset();
    std::vector<GraphNode*> m_nodes;
    std::map<std::pair<int, int>, GraphEdge*> m_edges;
};


void dijkstra_shortest_path(DualGraph *g, GraphNode *s,
	std::map<GraphNode*, double>& dist, std::map<GraphNode*, GraphNode*>& previous);

void modified_dijkstra_shortest_path(DualGraph *g, GraphNode *s,
	std::map<GraphNode*, double>& dist, std::map<GraphNode*, GraphNode*>& previous);

DualGraph* shortest_path_tree(DualGraph *g, std::map<GraphNode*, GraphNode*>& previous);




	
	
	


