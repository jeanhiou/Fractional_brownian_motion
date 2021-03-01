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

VectorXd c_n(int n, double H = 0.5){
  VectorXd c(n+1);
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

struct recu{

  recu(VectorXd d_pre): d_pre(d_pre) {};

  recu operator()(MatrixXd,double,VectorXd,int );

  VectorXd operator()();

private:
  VectorXd d_pre;
};

recu recu::operator()(MatrixXd F_n, double gamma,VectorXd d_prec,int n)
{
    VectorXd new_d_n(n+1);
    VectorXd d_plus(1);
    d_plus(0) = gamma;
    d_pre = d_pre - gamma * F_n * d_pre;
    new_d_n << d_pre, d_plus;
    recu new_d(new_d_n);
    return new_d;
};

VectorXd recu::operator()(){
  return d_pre;
};


double tau(VectorXd cn,MatrixXd F_n, VectorXd dn){
  return cn.dot( F_n * dn);
};


template <typename TGen>
double Hosking_method(TGen &gen, int N, double H = 0.5){

  double cova0 = covariance(1,H);
  double sigma_n = 1- pow(cova0,2);
  double mu_n = cova0;
  double tau0 = pow(cova0,2);
  double tau_n = tau0;
  MatrixXd F_n(1,1);
  F_n << 1;
  VectorXd d_pre;
  d_pre <<  1;
  VectorXd cn(1);
  cn << 1.;
  recu d(d_pre);
  double gamma = (cova0-tau0)/sigma_n;
  for (int j = 1;j<N;j++){
    sigma_n = sigma_n  - gamma/sigma_n;
    tau_n = tau(c_n(j,H),F_matrix(j),d(F_matrix(j),gamma,d(),j)());
    gamma = (covariance(j+1,H)-tau_n)/sigma_n;

  };
  return 1.;
};
