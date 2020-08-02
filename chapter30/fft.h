#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

std::vector<std::complex<double>> recursiveFFT(std::vector<double>& a);
std::vector<std::complex<double>> iterativeFFT(std::vector<double>& a);

#endif  /* FFT_H */
