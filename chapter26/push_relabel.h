#include <vector>

#ifndef PUSH_RELABEL_H_
#define PUSH_RELABEL_H_

struct Node {
  int _h = 0;
  int _e = 0;
};

class PushRelabel {
public:
  PushRelabel(std::vector<std::vector<int>>& edges, int vertice_count, 
			  int source, int target):
	_edges(edges), _source(source), _target(target) {
	_nodes.resize(vertice_count);
	_residual_network.resize(vertice_count);
	for (int i = 0; i < vertice_count; ++i) {
	  _nodes[i] = new Node;
	  _residual_network[i].resize(vertice_count);
	}
  }
  int Exec();
private:  
  void InitializePreflow();
  void Push(int u, int v);
  void Relabel(int u);
  int Capacity(int u, int v);
  
private:
  std::vector<Node *> _nodes;
  std::vector<std::vector<int>> _edges;
  std::vector<std::vector<int>> _residual_network;
  int _source;
  int _target;
};

#endif    // PUSH_RELABEL_H_
