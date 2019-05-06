#include "RooUnfoldHelpers.h"

#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>

namespace RooUnfolding {
  void printVector(const char* name, const TVectorD& vec){
    std::cout << name << std::endl;
    for(int i=0; i<vec.GetNrows(); ++i){
      std::cout << vec[i] << std::endl;
    }
  }

}
