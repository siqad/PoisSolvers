extern "C" {
  #include <cblas.h>
}
#include <iostream>

int main( void )
{
  int i=0;
  double A[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};
  double B[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};
  double C[9] = {.5,.5,.5,.5,.5,.5,.5,.5,.5};
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,3,3,2,1,A, 3, B, 3,2,C,3);

  for(i=0; i<9; i++){
    std::cout << " " << C[i] << " ";
  }
  std::cout << std::endl;
}

//compiles with: g++ -o test test_cblas_dgemm.cpp -I/usr/include/ -L/usr/lib -lopenblas -lpthread -lgfortran -llapack
