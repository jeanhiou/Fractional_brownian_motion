#include<iostream>
#include "fbm.h"

int main(){

  random_device rd;
  mt19937_64 gen(rd());
  double H = 0.5 ;
  double L = 1 ;
  int cum = 1 ;
  long n = 10 ;

  cout << covariance(1,0.5) << endl;

  cout << "methode hosking" << endl;
  path<double> fbm1 =  hosking(gen,n) ;
  save_fichier_path("fbm_hosking_brown.txt",fbm1);
  path<double> fbm2 =  hosking(gen,n,0.3) ;
  save_fichier_path("fbm_hosking_brownh3.txt",fbm2);
  path<double> fbm3 =  hosking(gen,n,0.1) ;
  save_fichier_path("fbm_hosking_brownh1.txt",fbm3);
  cout << endl;
  cout << " methode Cholesky " << endl;
  path<double> cho = Cholesky(gen,n,0.2);
  save_fichier_path("fbm_cholesky.txt",cho);
  path<double> cho2 = Cholesky(gen,n,0.2);
  save_fichier_path("fbm_cholesky2.txt",cho2);
  path<double> cho3 = Cholesky(gen,n,0.2);
  save_fichier_path("fbm_cholesky3.txt",cho3);


  return 0;
};
