#include "strassen.h"

#include <cassert>
#include <cstdio>
#include <algorithm>

using std::vector;
using std::max;

vector<vector<int>> Strassen::compute(const vector<vector<int>>& a, 
                                      const vector<vector<int>>& b) {
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
        vector<vector<int>> adapt_a(new_len, vector<int>(new_len, 0));
        vector<vector<int>> adapt_b(new_len, vector<int>(new_len, 0));
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                adapt_a[i][j] = a[i][j];
                adapt_b[i][j] = b[i][j];
            }
        }
        vector<vector<int>> m = computeInternal(adapt_a, adapt_b);
        vector<vector<int>> res(row, vector<int>(col));
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                res[i][j] = m[i][j];
            }
        }

        return res;
    }
}

void add(const vector<vector<int>>& a, const vector<vector<int>>& b, vector<vector<int>>& res) {
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a[0].size(); ++j) {
            res[i][j] = a[i][j] + b[i][j];
        }
    }
}

void sub(const vector<vector<int>>& a, const vector<vector<int>>& b, vector<vector<int>>& res) {
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a[0].size(); ++j) {
            res[i][j] = a[i][j] - b[i][j];
        }
    }
}

vector<vector<int>> Strassen::computeInternal(const vector<vector<int>>& a, 
                                              const vector<vector<int>>& b) {
    vector<vector<int>> res(a.size(), vector<int>(a.size(), 0));
    if (a.size() == 1) {
        res[0][0] = a[0][0] * b[0][0];
    } else {
        int half = a.size() >> 1;
        vector<int> half_vec(half, 0);
        vector<vector<int>> c11(half, half_vec);
        vector<vector<int>> c12(half, half_vec);
        vector<vector<int>> c21(half, half_vec);
        vector<vector<int>> c22(half, half_vec);

        vector<vector<int>> s1(half, half_vec);
        vector<vector<int>> s2(half, half_vec);
        vector<vector<int>> s3(half, half_vec);
        vector<vector<int>> s4(half, half_vec);
        vector<vector<int>> s5(half, half_vec);
        vector<vector<int>> s6(half, half_vec);
        vector<vector<int>> s7(half, half_vec);
        vector<vector<int>> s8(half, half_vec);
        vector<vector<int>> s9(half, half_vec);
        vector<vector<int>> s10(half, half_vec);
        
        vector<vector<int>> a11(half, half_vec);
        vector<vector<int>> a12(half, half_vec);
        vector<vector<int>> a21(half, half_vec);
        vector<vector<int>> a22(half, half_vec);
        vector<vector<int>> b11(half, half_vec);
        vector<vector<int>> b12(half, half_vec);
        vector<vector<int>> b21(half, half_vec);
        vector<vector<int>> b22(half, half_vec);

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

        sub(b12, b22, s1);
        add(a11, a12, s2);
        add(a21, a22, s3);
        sub(b21, b11, s4);
        add(a11, a22, s5);
        add(b11, b22, s6);
        sub(a12, a22, s7);
        add(b21, b22, s8);
        sub(a11, a21, s9);
        add(b11, b12, s10);

        vector<vector<int>> p1 = computeInternal(a11, s1);
        vector<vector<int>> p2 = computeInternal(s2, b22);
        vector<vector<int>> p3 = computeInternal(s3, b11);
        vector<vector<int>> p4 = computeInternal(a22, s4);
        vector<vector<int>> p5 = computeInternal(s5, s6);
        vector<vector<int>> p6 = computeInternal(s7, s8);
        vector<vector<int>> p7 = computeInternal(s9, s10);

        add(p1, p2, c12);
        add(p3, p4, c21);
        add(p5, p4, c11);
        sub(c11, p2, c11);
        add(c11, p6, c11);
        add(p5, p1, c22);
        sub(c22, p3, c22);
        sub(c22, p7, c22);
        
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
    Strassen s;
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
