// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"

PermittivityDrude::PermittivityDrude(double gamma, double omega_p) {
  // set parameters
  this->gamma = gamma;
  this->omega_p = omega_p;
};

// constructor for drude model from .ini file
PermittivityDrude::PermittivityDrude(std::string input_file)
    : Permittivity(input_file) {

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "drude");

  // read parameters
  this->gamma = root.get<double>("Permittivity.gamma");
  this->omega_p = root.get<double>("Permittivity.omega_p");
};

// calculate the permittivity
std::complex<double> PermittivityDrude::epsilon(double omega) {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = 1.0 - omega_p * omega_p /
                     std::complex<double>(omega * omega, gamma * omega);

  return result;
};

// calculate the permittivity scaled by omega
std::complex<double> PermittivityDrude::epsilon_omega(double omega) {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = omega - omega_p * omega_p / (omega + I * gamma);

  return result;
};
// getter method for damping coefficient
double PermittivityDrude::get_gamma() { return this->gamma; };

// getter method for plasma frequency
double PermittivityDrude::get_omega_p() { return this->omega_p; }
