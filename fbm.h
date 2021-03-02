#include<iostream>
#include <fstream>

#include "fft.h"
#include <math.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include<random>
#include <complex>

using namespace std;
using namespace std::complex_literals;


template <typename T>
struct state {
    double time;
    T value;
};

template <typename T>
std::ostream & operator<<(std::ostream & o, state<T> const & s) {
    return o << s.time << "\t" << s.value;
}

template <typename T>
struct path : protected std::vector<state<T>> {
    using vec = std::vector<state<T>>;  // alias de nom
    using vec::vec;             // constructeur de la classe vector utilisable
    using vec::operator[];      // opérateur [] utilisable (public)
    using vec::begin;           // itérateurs utilisables (for-range loop)
    using vec::end;
    using vec::size;            // utile !
};

template <typename T>
std::ostream & operator<<(std::ostream & o, path<T> const & p) {
    for (auto const & st : p)
        o << st << std::endl;
    return o << std::endl;
}

template<typename T>
void save_fichier_path(string filename, path<T> & p)
{
    ofstream of(filename); // ouverture du fichier
    for (int i=0; i<p.size() ; i++)
        of << p[i].time << " " << p[i].value << endl;
    of.close();           // fermature du fichier
};


double covariance(int i, double H) {
  if (i == 0) return 1;
  else return (pow(i-1,2*H)-2*pow(i,2*H)+pow(i+1,2*H))/2;
};


template<typename TGen>
inline unsigned p2(unsigned n) { unsigned u = 1; return u <<= n; }

template<typename TGen>
path<double> hosking(TGen & gen, long n, double H = 0.5, double L = 1. , int cum = 1) {

  long i, j ;
  int  m = pow(2,n);

  vector<double> output(m);

  double *phi = (double *) calloc(m, sizeof(double));
  double *psi = (double *) calloc(m, sizeof(double));
  double *cov = (double *) calloc(m, sizeof(double));
  double v, scaling;

  std::normal_distribution<double> G(0,1);

  output[0] = G(gen);
  v = 1;
  phi[0] = 0;
  for (i=0; i<m; i++)
    cov[i] = covariance(i,H);

  /* simulation */
  for(i=1; i<m; i++) {
    phi[i-1] = cov[i];
    for (j=0; j<i-1; j++) {
      psi[j] = phi[j];
      phi[i-1] -= psi[j]*cov[i-j-1];
    }
    phi[i-1] /= v;
    for (j=0; j<i-1; j++) {
      phi[j] = psi[j] - phi[i-1]*psi[i-j-2];
    }
    v *= (1-phi[i-1]*phi[i-1]);

    output[i] = 0;
    for (j=0; j<i; j++) {
      output[i] += phi[j]*output[i-j-1];
    }
    output[i] += sqrt(v)*G(gen);
  }

  /* rescale to obtain a sample of size 2^(*n) on [0,L] */
  scaling = pow(L/m,H);
  for(i=0;i<m;i++) {
    output[i] = scaling*(output[i]);
    if (cum && i>0) {
      output[i] += output[i-1]; } ;
    }
  path<double> sortie(m);
  for (int i = 0;i<m;i++){
    sortie[i] = { i / float(m-1) , output[i] };
  };

  free(phi);
  free(psi);
  free(cov);
  return sortie;
};

template<typename TGen>
path<double> Cholesky(TGen & gen, long n, double H = 0.5, double L = 1. , int cum = 1){
  int m = pow(2,n);
  std::normal_distribution<double> G(0,1);
  Eigen::MatrixXd M_final=Eigen::MatrixXd::Zero(m,m);
  M_final(0,0) = sqrt(covariance(0,H));
  Eigen::VectorXd Vector_simul_V = Eigen::VectorXd::Zero(m);
  Eigen::VectorXd output = Eigen::VectorXd::Zero(m);
  Eigen::VectorXd new_line = Eigen::VectorXd::Zero(2);
  for (int i = 0; i< m ; i++)
  {
    Vector_simul_V(i) = G(gen);
  };

  M_final(1,0) = covariance(1,H)/M_final(0,0);
  M_final(1,1) = sqrt( covariance(0,H) - M_final(1,0));

  output(0) = Vector_simul_V(0);
  output(1) = Vector_simul_V(0)*M_final(1,0) + Vector_simul_V(1)*M_final(1,1);

  for (int i = 2 ; i < m ; i++ ){
    new_line.resize(i+1);
    // creation de la nouvelle ligne //
    new_line(0) = covariance(i+1,H)/M_final(0,0);
    for (int j = 1 ; j<i;j++){
      double somme = 0;
      for (int k = 0; k <j;k++){
        somme += M_final(i,k)*M_final(j,k);
      };
      new_line(j) = 1/M_final(j,j) * ( covariance(i+1-j,H) - somme);
    };
    double somme = 0;
    for (int k = 0 ; k<i+1; k++){
      somme += pow(new_line(k),2);
    };
    new_line(i) = sqrt( covariance(0,H) - somme ) ;

    // final simulation  //
    double somme_final = 0;
    for (int j = 0; j<i+1;j++){
      somme_final += new_line(j) * Vector_simul_V(j);
      M_final(i,j) = new_line(j);
    };
    output(i) = somme_final ;
  };
  double scaling = pow(L/m,H);
  for(int i=0;i<m;i++) {
    output[i] = scaling*(output[i]);
    if (cum && i>0) {
      output[i] += output[i-1]; } ;
  };
  path<double> sortie(m);
  for (int i = 0;i<m;i++){
    sortie[i] = { i / float(m-1) , output[i] };
  };
  return sortie;
};

Eigen::VectorXd eigenvalues(long n, double H) {
  int size = pow(2,n+1);
  Complex ligne_C[size];
  Complex z1 = 1i;
  for (int i = 0 ; i<size;i++){
    if (i<=pow(2,n))
    {
      ligne_C[i] = 0. *z1 +  covariance(i,H);
    }
    else
    {
      ligne_C[i] = ligne_C[size-i];
    };
  };
  CArray data(ligne_C,size);
  fft(data);
  Eigen::VectorXd eigen(size);
  for (int i = 0; i< size; i++){
    eigen(i) = real(data[i]);
  }
  return eigen;
};

template<typename TGen>
path<double> David_and_Harte(TGen & gen, long n, double H = 0.5, double L = 1. , int cum = 1)
{
  long m = pow(2,n);

  Eigen::VectorXd eigenvalue = eigenvalues(n,H);

  Complex W[m];

  std::normal_distribution<double> G(0,1);
  std::complex<double> z1 = 1i;

  W[0] = G(gen);
  W[m-1] = G(gen);
  for (int i = 1 ; i< m/2 ; i++ ){
    double v1 = G(gen);
    double v2 = G(gen);
    W[i]   = 1/sqrt(2) * ( v1 + z1 *v2);
    W[m-i] = 1/sqrt(2) * ( v1 - z1 *v2);
  };

  for (int i = 0;i<m;i++){
  W[i] = eigenvalue[i] * W[i] / sqrt(m);
  }

  CArray data(W,m);

  ifft(data);
  Eigen::VectorXd output(m/2);
  double scaling = pow(L*2/m,H);
  for(int i=0; i<m/2; i++) {
    output[i] = scaling * (real(data[i]));
    if (cum && i>0)
    {
      output[i] += output[i-1];
    };
  };
  path<double> path_david(m/2);
  for (int i = 0; i<m/2;i++){
    path_david[i] = {i*2./m, output[i]};
  };
return path_david;
};
