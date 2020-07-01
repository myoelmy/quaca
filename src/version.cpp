#include "version.h"


void print_quaca_build_info(const char *prefix, std::ostream &stream) {
#ifdef GIT_HEAD
  stream << prefix << " git head: " << GIT_HEAD << "\n";
#endif
#ifdef GIT_BRANCH
  stream << prefix << " git branch: " << GIT_BRANCH << "\n";
#endif
}
