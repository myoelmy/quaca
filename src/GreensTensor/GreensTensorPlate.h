#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include "GreensTensor.h"
#include "Permittivity/PermittivityFactory.h"
#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>
// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

//! The class of the Green's tensor above a flat macroscopic surface
class GreensTensorPlate : public GreensTensor {
private:
  // permittivity is needed to describe the surface's response
  Permittivity *permittivity;

  double delta_cut; // numerical cut-off of the kappa integration
  vec::fixed<2> rel_err = {NAN, NAN};
  double za;

public:
  // constructors
  GreensTensorPlate(double v, double za, double beta,
                    Permittivity *permittivity, double delta_cut,
                    vec::fixed<2> rel_err);
  GreensTensorPlate(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrands
  static double integrand_1d_k(double kx, void *opts);
  static double integrand_2d_k(double ky, void *opts);

  // getter functions
  std::complex<double> get_epsilon(double omega) {
    return this->permittivity->epsilon(omega);
  };
  double get_za() { return this->za; }
  double get_delta_cut() { return this->delta_cut; }
  double get_rel_err_0() { return this->rel_err(0); }
  double get_rel_err_1() { return this->rel_err(1); }
};

#endif // GREENSTENSORPLATE_H
