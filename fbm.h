#include<iostream>
#include<random>
#include <cmath>
#include <eigen3/Eigen/dense>

using namespace std;
using namespace Eigen;

double covariance(int i, double H){
  if (i ==0) return 1;
  else return 0.5*(pow(i-1,2*H)-2*pow(i,2*H)+pow(i+1,2*H));
};

Vectorxd c_n(int n, double H = 0.5){
  Vectorxd c(n+1);
  for (int i = 0;i<n+1;i++){
    c(i) = covariance(i+1,H);
  };
  return c;
};

double indicatrice(int i, int j, int n){
  if ( i = n - j ) return 1.;
  else return 0;
};

MatrixXd F_matrix(int n ){
  MatrixXd F_ma(n+1,n+1);
  for (int i=0;i<n+1;i++){
    for (int j=0;j<n+1;j++){
      F_ma(i,j) = indicatrice(i,j,n);
    };
  };
  return F_ma;
};

double tau(Vectorxd cn,MatrixXd F_n, Vectorxd dn){
  return cn.dot( F_n * dn);
};

Vectorxd d_n( Vectorxd d_pre, MatrixXd F_n,double gamma,int n){
  d_plus(0) = gamma;
  if (n==0){
    Vectorxd d_n(n+1);
    d_n(0)=1;
    return d_n;
  }
  else{
    Vectorxd new_d_n(n+1);
    Vectorxd d_plus(1);
    d_plus(0) = gamma;
    d_pre = d_pre - gamma * F_n * d_pre;
    new_d_n << d_pre, d_plus;


  }
}

template<tynename TGen>
double Hosking_method(Tgen &gen, int N, double H = 0.5){

  double cova0 = covariance(1,H);
  double sigma_n = 1- pow(cova0,2);
  double mu_n = cova0;
  double tau_n = pow(cova0,2);

  MatrixXd F_n(1);
  F_n(0) = 1;
  Vectorxd d_pre(1);
  d_pre(0) = 1;


  for (int j = 1;j<N;j++){
    double gamma = (covariance(j+1,H)-tau_n)/sigma_n;
    sigma_n = sigma_n  - gamma/sigma_n;
    MatrixXd F_n = F_matrix(j);
    VectorXd d_pre = d_n(d_pre,F_n,gamma,j);
    VectorXd cn = c_n(j,H);
    double tau_n = tau(cn,F_n,d_pre);
  };
  return 1.;
};








    }
  }


}
