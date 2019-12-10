// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "OhmicMemoryKernel.h"

// constructor for ohmic memory kernel
OhmicMemoryKernel::OhmicMemoryKernel(double gamma): gamma(gamma)
{
};

// constructor for ohmic memory kernel from .ini file
OhmicMemoryKernel::OhmicMemoryKernel(std::string input_file)
{
    // Create a root
    pt::ptree root;

    // Load the ini file in this ptree
    pt::read_ini(input_file, root);

    // read damping coefficient
    this->gamma = root.get<double>("MemoryKernel.gamma");
};

// return mu(omega) for defined memory kernel
std::complex<double> OhmicMemoryKernel::mu(double omega)
{
    const std::complex<double> gammac(this->gamma,0E0);
    return gammac;
};

// getter method for damping kernel
double OhmicMemoryKernel::get_gamma()
{
    return this->gamma;
};
