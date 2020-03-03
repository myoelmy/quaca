#ifndef LOOPERFACTORY_H
#define LOOPERFACTORY_H

#include "Looper.h"

class LooperFactory {
public:
  static Looper *create(std::string type);
};

#endif // LOOPERFACTORY_H