// compile:
// clang++ -g -std=c++11 push_relabel.cc -o push_relabel -DDebug

#include "push_relabel.h"

#include <limits>
#include <algorithm>

using std::vector;

int PushRelabel::exec() {
  InitializePreflow();
  do {
	bool relabel_flag = false;
	for (int i = 0; i < _nodes.size(); ++i) {
	  if (i == _source || i == _target || _nodes[i]->_e == 0) continue;
	  bool need_relabel_flag = false;
	  for (int j = 0; j < _nodes.size(); ++j) {
		if (capacity(i, j) == 0) continue;
		if (_nodes[i]->_h > _nodes[j]->_h) {
		  need_relabel_flag = false;
		  break;
		} else {
		  need_relabel_flag = true;
		}
	  }

	  if (need_relabel_flag) {
		relabel_flag = true;
		Relabel(i);
	  }
	}

	bool push_flag = false;
	for (int i = 0; i < _nodes.size(); ++i) {
	  if (_nodes[i]->_e == 0) continue;
	  for (int j = 0; j < _nodes.size(); ++j) {
		if (capacity(i, j) == 0 || _nodes[i]->_h != _nodes[j]->_h + 1)
			continue;
		Push(i, j);
		push_flag = true;
	  }
	}

	if (!relabel_flag && !push_flag) break;
  } while (1); 
	
  return _nodes[_target]->_e;
}

void PushRelabel::InitializePreflow() {
  _nodes[_source]->_h = _nodes.size();
  for (int i = 0; i < _nodes.size(); ++i) {
	if (_edges[_source][i] == 0) continue;
	_nodes[i]->_e = _edges[_source][i];
	_nodes[_source]->_e -= _edges[_source][i];
	_residual_network[_source][i] = _edges[_source][i];
  }
}

void PushRelabel::Push(int u, int v) {
  double delta = std::min(_nodes[u]->_e, capacity(u, v));
  if (_edges[u][v] > 0) _residual_network[u][v] += delta;
  else _residual_network[v][u] -= delta;
  _nodes[u]->_e -= delta;
  _nodes[v]->_e += delta;
}

void PushRelabel::Relabel(int u) {
  int label = std::numeric_limits<int>::max();
  for (int i = 0; i < _nodes.size(); ++i) {
	if (capacity(u, i) > 0) label = std::min(label, _nodes[i]->_h);
  }

  _nodes[u]->_h = label + 1;
}

int PushRelabel::capacity(int u, int v) {
  if (_edges[u][v] > 0) return _edges[u][v] - _residual_network[u][v];
  else if (_edges[v][u] > 0) return _residual_network[v][u];
  
  return 0;
}

#if Debug
int main(int argc, char *argv[]) {
  vector<vector<int>> edges(6, vector<int>(6, 0));
  edges[0][1] = 15;
  edges[0][3] = 4;
  edges[1][2] = 12;
  edges[4][1] = 5;
  edges[2][3] = 3;
  edges[3][4] = 10;
  edges[2][5] = 7;
  edges[4][5] = 10;
  PushRelabel flow(edges, 6, 0, 5);
  printf("%d\n", flow.exec());

  return 0;
}
#endif
