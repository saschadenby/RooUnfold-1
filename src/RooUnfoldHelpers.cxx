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

  void add(TMatrixD& target, const TMatrixD& addition){
    for(int i=0; i<addition.GetNrows(); ++i){
      for(int j=0; j<addition.GetNcols(); ++j){
        target(i,j) = target(i,j) + addition(i,j);
      }
    }
  }
}
