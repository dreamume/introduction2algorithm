// clang++ -g -DDebug -std=c++11 euclid.cc -o euclid

#include "euclid.h"

#include <cstdio>
#include <vector>

void CEuclid::extended_euclid(int a, int b) {
  if (b == 0) {
    d_ = a;
    x_ = 1;
    y_ = 0;
    return;
  }

  extended_euclid(b, a % b);
  int x = y_;
  int y = x_ - (a / b) * y_;
  x_ = x;
  y_ = y;
}

void CEuclid::modular_linear_equatiion_solver(int a, int b, int n) {
  extended_euclid(a, n);
  if (b % d_ == 0) {
    int res = x_ * (b / d_) % n;
    while (res < 0) res += n;
    for (int i = 0; i < d_; ++i)
      printf("%d ", ((res + i * (n / d_)) % n));
    printf("\n");
  } else {
    printf("no solutions\n");
  }
}

int CEuclid::modular_exponentiation(int a, int b, int n) {
  int c = 0;
  int d = 1;
  std::vector<bool> bits;
  for (int i = 0; i <= 32; ++i) {
    if (((long long)1 << i) <= b) bits.push_back((long long)b & ((long long)1 << i));
    else break;
  }

  for (int i = bits.size() - 1; i >= 0; --i) {
    c = 2 * c;
    d = (d * d) % n;
    if (bits[i]) {
      c = c + 1;
      d = (d * a) % n;
    }
  }

  return d;
}

int CEuclid::modular_exponentiation2(int a, int b, int n) {
  int e = a % n;
  int d = 1;

  int index = 0;
  while ((1ULL << index) <= b) {
    if ((0x1 << index) & b) d = (d * e) % n;
    e = (e * e) % n;
    ++index;
  }

  return d;
}

#if Debug
int main(int argc, char *argv[]) {
  CEuclid euclid;
  //euclid.extended_euclid(899, 493);
  // printf("d is %d, x is %d, y is %d, ax + by = %d\n", 
  //        euclid.d(), euclid.x(), euclid.y(), 899 * euclid.x() + 493 * euclid.y());
  euclid.modular_linear_equatiion_solver(14, 30, 100);
  printf("module exponentiation a is 7, b is 560, n = 561, result is %d\n",
         euclid.modular_exponentiation2(7, 560, 561));
  
  return 0;
}
#endif
