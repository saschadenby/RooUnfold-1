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

  void printMatrix(const TMatrixD& m, const char* name, const char* format, Int_t cols_per_sheet){
  // Print the matrix as a table of elements.
  // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
  // By default the format "%11.4g" is used to print one element.
  // One can specify an alternative format with eg
  //  format ="%6.2f  "

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i = 1; i <= nrows; i++) {
      printf("%4d |",i+rowlwb-1);
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
        printf(format,m(i+rowlwb-1,j+collwb-1));
      printf("\n");
    }
  }
  printf("\n");
}


TMatrixD& ABAT (const TMatrixD& a, const TMatrixD& b, TMatrixD& c){
  // Fills C such that C = A * B * A^T. Note that C cannot be the same object as A.
  TMatrixD d (b, TMatrixD::kMultTranspose, a);
  c.Mult (a, d);
  return c;
}

TMatrixD& ABAT (const TMatrixD& a, const TVectorD& b, TMatrixD& c){
  // Fills C such that C = A * B * A^T, where B is a diagonal matrix specified by the vector.
  // Note that C cannot be the same object as A.
  TMatrixD d (TMatrixD::kTransposed, a);
  d.NormByColumn (b, "M");
  c.Mult (a, d);
  return c;
}




}
