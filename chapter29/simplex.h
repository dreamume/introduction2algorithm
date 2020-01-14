#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>

struct LinearProgrammingData {
LinearProgrammingData(std::vector<std::vector<double>>& A,
                      std::vector<double>& b,
                      std::vector<double>& c)
  : _A(std::vector<std::vector<double>>(A.size(), std::vector<double>(A.size() + A[0].size(), 0.0))), 
    _b(b), _c(c) {
    for (int i = 0; i < A.size(); ++i)
      for (int j = 0; j < A[0].size(); ++j)
        _A[i][j] = A[i][j];

    for (int i = 0; i < b.size(); ++i) _b[i] = b[i];
    _c.resize(_A[0].size());
    _non_basic_variable_indexes.resize(A[0].size());
    for (int i = 0; i < A[0].size(); ++i) _non_basic_variable_indexes[i] = i;
    _basic_variable_indexes.resize(A.size());
    for (int i = 0; i < A.size(); ++i) 
      _basic_variable_indexes[i] = A[0].size() + i;
  }
LinearProgrammingData(int num_rows, int num_cols)
  : _A(std::vector<std::vector<double>>(num_rows, std::vector<double>(num_rows + num_cols, 0.0))),
    _b(std::vector<double>(num_rows, 0.0)),
    _c(std::vector<double>(num_rows + num_cols, 0.0)) {
    _non_basic_variable_indexes.resize(num_cols);
    for (int i = 0; i < num_cols; ++i) _non_basic_variable_indexes[i] = i;
    _basic_variable_indexes.resize(num_rows);
    for (int i = 0; i < num_rows; ++i) 
      _basic_variable_indexes[i] = num_cols + i;
  }
    
  std::vector<std::vector<double>> _A;
  std::vector<double> _b;
  std::vector<double> _c;
  std::vector<int> _non_basic_variable_indexes;
  std::vector<int> _basic_variable_indexes;
  double _v = 0.0;
  bool _is_feasible = true;
  bool _is_bound = true;
};

class Simplex {
 public:
  /* A, b, c is a standard form */
  Simplex(std::vector<std::vector<double>>& A, 
          std::vector<double>& b,
          std::vector<double>& c)
    : _data(A, b, c) {}
  void Run(std::vector<double>& solution);
  bool IsFeasible() { return _data._is_feasible; }
  bool IsBound() { return _data._is_bound; }
  double ObjectiveFunctionValue() { return _data._v; }

 private:
  int Initialize(std::vector<std::vector<double>>& A,
                 std::vector<double>& b,
                 std::vector<double>& c);
  void Pivot(LinearProgrammingData& data, int leave_index, int enter_index);
  void LoopInternal(LinearProgrammingData& data);

 private:
  LinearProgrammingData _data;
};

#endif  /* SIMPLEX_H */
