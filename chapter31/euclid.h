class CEuclid {
 public:
  CEuclid() {}
  void extended_euclid(int a, int b);
  void modular_linear_equatiion_solver(int a, int b, int n);
  int modular_exponentiation(int a, int b, int n);
  int modular_exponentiation2(int a, int b, int n);
  int gcd(int a, int b);
  int x() { return x_; }
  int y() { return y_; }
  int d() { return d_; }
 private:
  int x_;
  int y_;
  int d_;
};
