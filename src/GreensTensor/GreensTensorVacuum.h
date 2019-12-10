#ifndef GREENSTENSORVACUUM_H
#define GREENSTENSORVACUUM_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"

class GreensTensorVacuum : public GreensTensor
{

public:

  GreensTensorVacuum(double v, double beta);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, double omega, double kv, Options opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, double omega, Options opts);
  static double integrand_k_1d(double kv, void* options);

};


#endif // GREENSTENSORVACUUM_H
