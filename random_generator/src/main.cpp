#include "Engine.hpp"
#include "versioning.h"
#include <iostream>

int main(int argc, char *argv[]) {

  std::cout << PROJECT_NAME << " version " << PROJECT_VER << std::endl;
  Engine engine(argc, argv);
  return engine.exec();
}
