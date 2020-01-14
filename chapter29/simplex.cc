// compile: clang++ -DDebug -g -std=c++11 simplex.cc -o simplex

#include "simplex.h"

#include <algorithm>
#include <utility>
#include <limits>

void Simplex::Run(std::vector<double>& solution) {
  int result = Initialize(_data._A, _data._b, _data._c);
  if (result != 0) {
    if (result == -1) _data._is_bound = false;
    else if (result == -2) _data._is_feasible = false;
    return;
  }

  LoopInternal(_data);
  for (int i = 0; i < _data._A[0].size(); ++i) {
    auto it = std::find(_data._basic_variable_indexes.begin(),
                        _data._basic_variable_indexes.end(),
                        i);
    if (it != _data._basic_variable_indexes.end()) {
      int index = std::distance(_data._basic_variable_indexes.begin(), it);
      solution.push_back(_data._b[index]);
    } else {
      solution.push_back(0);
    }
  }
}

int Simplex::Initialize(std::vector<std::vector<double>>& A,
                         std::vector<double>& b,
                         std::vector<double>& c) {
  auto min_b = std::min_element(b.begin(), b.end());
  if (*min_b >= 0) return 0;
  LinearProgrammingData data(b.size(), c.size() - b.size() + 1);
  data._b = b;
  data._c[0] = -1;
  for (int i = 0; i < A.size(); ++i) {
    data._A[i][0] = -1;
    for (int j = 0; j < A[0].size(); ++j) data._A[i][j + 1] = A[i][j];
  }
  
  int index = std::distance(b.begin(), min_b);
  Pivot(data, index, 0);

  LoopInternal(data);
  if (!data._is_bound) return -1;
  auto it_0 = std::find(data._basic_variable_indexes.begin(), 
                        data._basic_variable_indexes.end(), 
                        0);
  int x0_index = 0;
  if (it_0 != data._basic_variable_indexes.end())
    x0_index = std::distance(data._basic_variable_indexes.begin(), it_0);
  if (it_0 == data._basic_variable_indexes.end() || data._b[x0_index] == 0) {
    if (it_0 != data._basic_variable_indexes.end()) {
      int enter = -1;
      for (int i = 0; i < data._non_basic_variable_indexes.size(); ++i) {
        if (data._c[data._non_basic_variable_indexes[i]] > 0) {
          enter = i;
          break;
        }
      }
      if (enter == -1) {
        for (int i = x0_index + 1; 
             i < data._non_basic_variable_indexes.size(); 
             ++i) {
          data._b[i - 1] = data._b[i];
          data._A[i - 1] = data._A[i];
          data._b.resize(data._b.size() - 1);
          data._A.resize(data._A.size() - 1);
        }
        data._basic_variable_indexes.erase(it_0);
      } else {
        Pivot(data, enter, x0_index);
        auto it_x0 = std::find(data._non_basic_variable_indexes.begin(), 
                               data._non_basic_variable_indexes.end(), 
                               0);
        data._non_basic_variable_indexes.erase(it_x0);
      }
    } else {
      auto it_x0 = std::find(data._non_basic_variable_indexes.begin(), 
                             data._non_basic_variable_indexes.end(), 
                             0);
      data._non_basic_variable_indexes.erase(it_x0);
    }
    
    for (int i = 0; i < data._A.size(); ++i)
      for (int j = 1; j < data._A[0].size(); ++j) 
        data._A[i][j] = data._A[i][j - 1];

    //for (int i = 1; i < data._c.size(); ++i) data._c[i - 1] = data._c[i];
    for (auto& index: data._basic_variable_indexes) --index;
    for (auto& index: data._non_basic_variable_indexes) --index;

    std::swap(data._A, _data._A);
    std::swap(data._b, _data._b);
    //std::swap(data._c, _data._c);
    std::swap(data._non_basic_variable_indexes, 
              _data._non_basic_variable_indexes);
    std::swap(data._basic_variable_indexes, _data._basic_variable_indexes);
    //std::swap(data._v, _data._v);
    std::swap(data._is_bound, _data._is_bound);
    std::swap(data._is_feasible, _data._is_feasible);

    return 0;
  }

  return -2;
}

void Simplex::LoopInternal(LinearProgrammingData& data) {
  // std::vector<int> count(data._non_basic_variable_indexes.size() + 
  //                        data._basic_variable_indexes.size(), 
  //                        0);
  while (true) {
    int enter = -1;
    for (int i = 0; i < data._non_basic_variable_indexes.size(); ++i) {
      if (data._c[data._non_basic_variable_indexes[i]] > 0) {
        // if (count[data._non_basic_variable_indexes[i]] == 0) {
        //   ++count[data._non_basic_variable_indexes[i]];
        // } else {
        //   data._is_bound = false;
        //   return;
        // }
        enter = i;
        break;
      }
    }
    if (enter == -1) break;

    int enter_index = data._non_basic_variable_indexes[enter];
    double delta = std::numeric_limits<double>::max();
    int leave = -1;
    for (int i = 0; i < data._basic_variable_indexes.size(); ++i) {
      if (data._A[i][enter_index] > 0) {
        if (data._b[i] / data._A[i][enter_index] < delta) {
          delta = data._b[i] / data._A[i][enter_index];
          leave = i;
        }
      }
    }

    if (delta == std::numeric_limits<double>::max()) {
      data._is_bound = false;
      return;
    }
    Pivot(data, leave, enter);
  }
}

void Simplex::Pivot(LinearProgrammingData& data, int leave, int enter) {
  /* Compute the coefficients of the equation for new basic variable x_enter */
  std::vector<std::vector<double>> new_A(data._A);
  std::vector<double> new_b(data._b);
  int leave_index = data._basic_variable_indexes[leave];
  int enter_index = data._non_basic_variable_indexes[enter];
  new_b[leave] = data._b[leave] / data._A[leave][enter_index];
  new_A[leave][enter_index] = 0;
  for (auto& j : data._non_basic_variable_indexes) {
    if (j == enter_index) continue;
    new_A[leave][j] = 
      data._A[leave][j] / data._A[leave][enter_index];
  }
  new_A[leave][leave_index] = 1 / data._A[leave][enter_index];

  // Compute the coefficients of the remaining constraints
  for (int i = 0; i < data._basic_variable_indexes.size(); ++i) {
    if (i == leave) continue;
    new_b[i] = data._b[i] - data._A[i][enter_index] * new_b[leave];
    new_A[i][enter_index] = 0;
    for (int& j : data._non_basic_variable_indexes) {
      if (j == enter_index) continue;
      new_A[i][j] = 
        data._A[i][j] - data._A[i][enter_index] * new_A[leave][j];
    }
    new_A[i][leave_index] = 
      0 - data._A[i][enter_index] * new_A[leave][leave_index];
  }

  // Compute the objective function
  data._v += data._c[enter_index] * new_b[leave];
  std::vector<double> new_c(data._c);
  new_c[enter_index] = 0;
  for (int& j : data._non_basic_variable_indexes) {
    if (j == enter_index) continue;
    new_c[j] = data._c[j] - data._c[enter_index] * new_A[leave][j];
  }
  new_c[leave_index] = 
    0 - data._c[enter_index] * new_A[leave][leave_index];

  // Compute new sets of basic and nonbasic variables.
  data._basic_variable_indexes[leave] = enter_index;
  data._non_basic_variable_indexes[enter] = leave_index;
  
  std::swap(data._b, new_b);
  std::swap(data._A, new_A);
  std::swap(data._c, new_c);
}

#if Debug
int main(int argc, char *argv[]) {
  // normal case:
  // optimal objective value should be 28
  std::vector<std::vector<double>> A{{1, 1, 3},
                                     {2, 2, 5},
                                     {4, 1, 2}};
  std::vector<double> b{30, 24, 36};
  std::vector<double> c{3, 1, 2};

  // degeneracy cases:
  // optimal objective value should be 16
  // std::vector<std::vector<double>> A{{1, 1, 0},
  //                                    {0, -1, 1}};
  // std::vector<double> b{8, 0};
  // std::vector<double> c{1, 1, 1};

  // optimal objective value should be 18
  // std::vector<std::vector<double>> A{{1, 4},
  //                                    {1, 2}};
  // std::vector<double> b{8, 4};
  // std::vector<double> c{3, 9};

  // optimal objective value should be 30750
  // std::vector<std::vector<double>> A{{1, 2, 0},
  //                                    {2, 3, -1},
  //                                    {-1, 0, 0},
  //                                    {0, -1, 0}};
  // std::vector<double> b{70, 100, -20, -25};
  // std::vector<double> c{750, 900, -450};

  // optimal objective value should be 2
  // std::vector<std::vector<double>> A{{1, 1, 0},
  //                                    {0, -1, 1}};
  // std::vector<double> b{1, 0};
  // std::vector<double> c{1, 1, 1};

  Simplex simplex(A, b, c);
  std::vector<double> solution;
  simplex.Run(solution);

  if (!simplex.IsFeasible()) {
    printf("the solution is not feasible!\n");
  } else if (!simplex.IsBound()) {
    printf("the solution is not boudned!\n");
  } else {
    double result = 0.0;
    for (int i = 0; i < c.size(); ++i) result += c[i] * solution[i];
    printf("the objective function value is %lf\n", result);
    for (int i = 0; i < solution.size(); ++i)
      printf("the %d variable is %lf\n", i, solution[i]);
  }

  return 0;
}
#endif
