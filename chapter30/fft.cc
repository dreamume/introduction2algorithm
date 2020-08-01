// compile: clang++ -DDebug -g -std=c++11 fft.cc -o fft

#include "fft.h"

#include <cmath>

std::vector<std::complex<double>> recursiveFFT(std::vector<double>& a) {
    if (a.size() == 1) return { std::complex<double>(a.front(), 0) };

    double PI = std::acos(-1);
    std::complex<double> wn(0, 2);
    wn *= PI / a.size();
    wn = std::exp(wn);
    std::complex<double> w(1.0f, 0.0f);
    std::vector<double> a0;
    std::vector<double> a1;
    for (int i = 0; i < a.size(); i += 2) {
        a0.push_back(a[i]);
        a1.push_back(a[i + 1]);
    }
    std::vector<std::complex<double>> y0 = recursiveFFT(a0);
    std::vector<std::complex<double>> y1 = recursiveFFT(a1);
    std::vector<std::complex<double>> res(a.size(), std::complex<double>(0, 0));
    for (int i = 0; i < a.size() / 2; ++i) {
        res[i] = y0[i] + w * y1[i];
        res[i + a.size() / 2] = y0[i] - w * y1[i];
        w *= wn;
    }

    return res;
}

#if Debug

#include <iostream>

int main(int argc, char *argv[]) {
    // normal case:
    // 求得的值似乎跟实际值偏差较大
    std::vector<double> a{-10, 1, -1, 7, 0, 0, 0, 0};
    std::vector<double> b{3, -6, 0, 8, 0, 0, 0, 0};

    std::vector<std::complex<double>> y0 = recursiveFFT(a);
    std::vector<std::complex<double>> y1 = recursiveFFT(b);
    std::vector<std::complex<double>> y;
    for (int i = 0; i < y0.size(); ++i) y.push_back(y0[i] * y1[i]);
    std::complex<double> w(1, 0);
    std::complex<double> wn(0, 2);
    double PI = std::acos(-1);
    wn *= PI / y.size();
    std::vector<std::complex<double>> wlist(y.size(), std::complex<double>(0, 0));
    wlist[0] = std::complex<double>(1, 0);
    for (int i = 1; i < y.size(); ++i) wlist[i] = wlist[i - 1] * wn;
    std::vector<std::complex<double>> worigin(wlist);
    std::vector<std::complex<double>> alist(y.size(), std::complex<double>(0, 0));
    for (int i = 0; i < y.size(); ++i) alist[0] += y[i];
    alist[0] /= y.size();
    for (int i = 1; i < y.size(); ++i) {
        alist[i] += y[0];
        for (int j = 1; j < y.size(); ++j) {
            alist[i] += y[j] * wlist[j];
            wlist[j] *= worigin[j];
        }
        alist[i] /= y.size();
    }

    for (int i = 0; i < alist.size(); ++i) std::cout << alist[i] << '\n';

    return 0;
}

#endif
