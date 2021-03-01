#include<iostream>
#include "fbm.h"

int main(){

  random_device rd;
  mt19937_64 gen(rd());
  double H = 0.5 ;
  double L = 1 ;
  int cum = 1 ;
  long n = 10 ;

  cout << "methode hosking" << endl;
  path<double> fbm1 =  hosking(gen,n) ;
  save_fichier_path("fbm_hosking_brown.txt",fbm1);
  path<double> fbm2 =  hosking(gen,n,0.3) ;
  save_fichier_path("fbm_hosking_brownh3.txt",fbm2);
  path<double> fbm3 =  hosking(gen,n,0.1) ;
  save_fichier_path("fbm_hosking_brownh1.txt",fbm3);

  return 0;
};
