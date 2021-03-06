#ifndef LOOPERZA_H
#define LOOPERZA_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperZa : public Looper {
public:
  // constructors
  LooperZa(double start, double end, int number_of_steps,
           const std::string &scale);
  LooperZa(const std::string &input_file);

  // calculate the the value of quantum friction
  double
  calculate_value(int step,
                  std::shared_ptr<Friction> quantum_friction) const override;

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // LOOPERZA_H
