#include<iostream>
#include "fbm.h"

int main(){

  random_device rd;
  mt19937_64 gen(rd());
  double H = 0.5 ;
  double L = 1 ;
  int cum = 1 ;
  long n = 10 ;
  path<double> p  = David_and_Harte(gen,n);
  save_fichier_path("David_and_Harte_essai.txt",p);


  return 0;
};
