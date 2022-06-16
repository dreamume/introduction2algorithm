#include "strassen.h"

#include <cassert>
#include <cstdio>
#include <algorithm>

using std::vector;
using std::max;

template <typename T> vector<vector<T>> Strassen<T>::compute(const vector<vector<T>>& a, 
                                                             const vector<vector<T>>& b) {
    assert(a.size() == b.size() && a[0].size() == b[0].size());

    int row = a.size();
    int col = a[0].size();
    int new_row = 1;
    int new_col = 1;
    
    while (new_row < row) new_row <<= 1;
    while (new_col < col) new_col <<= 1;
    if (new_row == row && new_col == col && row == col) {
        return computeInternal(a, b);
    } else {
        int new_len = max(new_row, new_col);
        vector<vector<T>> adapt_a(new_len, vector<T>(new_len, 0));
        vector<vector<T>> adapt_b(new_len, vector<T>(new_len, 0));
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                adapt_a[i][j] = a[i][j];
                adapt_b[i][j] = b[i][j];
            }
        }
        vector<vector<T>> m = computeInternal(adapt_a, adapt_b);
        vector<vector<T>> res(row, vector<T>(col));
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                res[i][j] = m[i][j];
            }
        }

        return res;
    }
}

template <typename T> void add(const vector<vector<T>>& a, const vector<vector<T>>& b, vector<vector<T>>& res) {
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a[0].size(); ++j) {
            res[i][j] = a[i][j] + b[i][j];
        }
    }
}

template <typename T> void sub(const vector<vector<T>>& a, const vector<vector<T>>& b, vector<vector<T>>& res) {
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a[0].size(); ++j) {
            res[i][j] = a[i][j] - b[i][j];
        }
    }
}

template <typename T> vector<vector<T>> Strassen<T>::computeInternal(const vector<vector<T>>& a, 
                                                                     const vector<vector<T>>& b) {
    vector<vector<T>> res(a.size(), vector<T>(a.size(), 0));
    if (a.size() == 1) {
        res[0][0] = a[0][0] * b[0][0];
    } else {
        int half = a.size() >> 1;
        vector<T> half_vec(half, 0);
        vector<vector<T>> c11(half, half_vec);
        vector<vector<T>> c12(half, half_vec);
        vector<vector<T>> c21(half, half_vec);
        vector<vector<T>> c22(half, half_vec);

        vector<vector<T>> s1(half, half_vec);
        vector<vector<T>> s2(half, half_vec);
        vector<vector<T>> s3(half, half_vec);
        vector<vector<T>> s4(half, half_vec);
        vector<vector<T>> s5(half, half_vec);
        vector<vector<T>> s6(half, half_vec);
        vector<vector<T>> s7(half, half_vec);
        vector<vector<T>> s8(half, half_vec);
        vector<vector<T>> s9(half, half_vec);
        vector<vector<T>> s10(half, half_vec);
        
        vector<vector<T>> a11(half, half_vec);
        vector<vector<T>> a12(half, half_vec);
        vector<vector<T>> a21(half, half_vec);
        vector<vector<T>> a22(half, half_vec);
        vector<vector<T>> b11(half, half_vec);
        vector<vector<T>> b12(half, half_vec);
        vector<vector<T>> b21(half, half_vec);
        vector<vector<T>> b22(half, half_vec);

        for (int i = 0; i < half; ++i) {
            for (int j = 0; j < half; ++j) {
                a11[i][j] = a[i][j];
                a12[i][j] = a[i][j + half];
                a21[i][j] = a[i + half][j];
                a22[i][j] = a[i + half][j + half];
                b11[i][j] = b[i][j];
                b12[i][j] = b[i][j + half];
                b21[i][j] = b[i + half][j];
                b22[i][j] = b[i + half][j + half];
            }
        }

        sub<T>(b12, b22, s1);
        add<T>(a11, a12, s2);
        add<T>(a21, a22, s3);
        sub<T>(b21, b11, s4);
        add<T>(a11, a22, s5);
        add<T>(b11, b22, s6);
        sub<T>(a12, a22, s7);
        add<T>(b21, b22, s8);
        sub<T>(a11, a21, s9);
        add<T>(b11, b12, s10);

        vector<vector<T>> p1 = computeInternal(a11, s1);
        vector<vector<T>> p2 = computeInternal(s2, b22);
        vector<vector<T>> p3 = computeInternal(s3, b11);
        vector<vector<T>> p4 = computeInternal(a22, s4);
        vector<vector<T>> p5 = computeInternal(s5, s6);
        vector<vector<T>> p6 = computeInternal(s7, s8);
        vector<vector<T>> p7 = computeInternal(s9, s10);

        add<T>(p1, p2, c12);
        add<T>(p3, p4, c21);
        add<T>(p5, p4, c11);
        sub<T>(c11, p2, c11);
        add<T>(c11, p6, c11);
        add<T>(p5, p1, c22);
        sub<T>(c22, p3, c22);
        sub<T>(c22, p7, c22);
        
        for (int i = 0; i < half; ++i) {
            for (int j = 0; j < half; ++j) {
                res[i][j] = c11[i][j];
                res[i][j + half] = c12[i][j];
                res[i + half][j] = c21[i][j];
                res[i + half][j + half] = c22[i][j];
            }
        }
    }

    return res;
}

#ifdef Debug
int main() {
    Strassen<int> s;
    vector<vector<int>> a = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    vector<vector<int>> b = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
    // vector<vector<int>> a{{1,3}, {7,5}};
    // vector<vector<int>> b{{6,8}, {4,2}};
    vector<vector<int>> res = s.compute(a, b);
    for (int i = 0; i < res.size(); ++i) {
        for (int j = 0; j < res[0].size(); ++j) {
            printf("%d%c", res[i][j], j == res[0].size() - 1 ? '\n' : ' ');
        }
    }

    return 0;
}
#endif
