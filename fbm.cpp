#include<iostream>
#include "fbm.h"

int main(){

  random_device rd;
  mt19937_64 gen(rd());

  int N = 10;
  double H = 0.5;

  std::cout <<  Hosking_method(gen,N,H) << std::endl;

  return 0;
};
