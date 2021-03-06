#ifndef PERMITTIVITYDRUDE_H
#define PERMITTIVITYDRUDE_H

#include "Permittivity.h"
#include <complex>

//! A Drude model permittivity
class PermittivityDrude : public Permittivity {
private:
  double omega_p; // plasma frequency
  double gamma;   // damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  explicit PermittivityDrude(const std::string &input_file);

  // calculate the permittivity
  std::complex<double> calculate(double omega) const override;

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> calculate_times_omega(double omega) const override;

  // getter methods
  double get_gamma() const { return this->gamma; };
  double get_omega_p() const { return this->omega_p; };

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // PERMITTIVITYDRUDE_H
