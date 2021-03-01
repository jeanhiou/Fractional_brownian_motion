#include<iostream>
#include <fstream>

#include <math.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include<random>
using namespace std;

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
