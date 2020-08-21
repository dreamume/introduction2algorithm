// clang++ -g -DDebug -std=c++11 euclid.cc -o euclid

#include "euclid.h"

#include <cstdio>

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

#if Debug
int main(int argc, char *argv[]) {
  CEuclid euclid;
  //euclid.extended_euclid(899, 493);
  // printf("d is %d, x is %d, y is %d, ax + by = %d\n", 
  //        euclid.d(), euclid.x(), euclid.y(), 899 * euclid.x() + 493 * euclid.y());
  euclid.modular_linear_equatiion_solver(14, 30, 100);
  
  return 0;
}
#endif
