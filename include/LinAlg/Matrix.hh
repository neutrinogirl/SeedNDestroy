//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_INCLUDE_SND_MATRIX_HH_
#define SND_INCLUDE_SND_MATRIX_HH_

#include <vector>
#include <iostream>

typedef std::vector<std::vector<double> > Matrix_t;

struct Matrix {

  Matrix_t M;
  std::size_t nrows, ncols;

  Matrix(std::size_t n, std::size_t m)
	  : nrows(n), ncols(m) {
	M = Matrix_t (n, std::vector<double>(m, 0));
  }

  std::vector<double>& operator[](std::size_t idx) {return M[idx];}

  void Print(){
	for(auto i=0; i<nrows; i++){
	  for(auto j=0; j<ncols; j++){
		std::cout << M[i][j] << " ";
	  }
	  std::cout << std::endl;
	}
  }

};
typedef struct Matrix Matrix;


struct DiagMatrix {

  Matrix_t dM;
  std::size_t dim;

  explicit DiagMatrix(std::size_t n)
	  : dim(n) {
	dM = Matrix_t (n, std::vector<double>(n, 0));
  }

  double& operator[](std::size_t idx) {return dM[idx][idx];}

  void Print(){
	for(auto i=0; i<dim; i++){
	  for(auto j=0; j<dim; j++){
		std::cout << dM[i][j] << " ";
	  }
	  std::cout << std::endl;
	}
  }

};
typedef struct DiagMatrix DiagMatrix;


#endif //SND_INCLUDE_SND_MATRIX_HH_
