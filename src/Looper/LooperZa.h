#ifndef LOOPERZA_H
#define LOOPERZA_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperZa : public Looper {
public:
  // constructors
  LooperZa(double start, double end, int number_of_steps, std::string scale);
  LooperZa(std::string input_file);

  // calculate the the value of quantum friction
  double calculate_value(int step, Friction* quantum_friction);
};

#endif // LOOPERZA_H