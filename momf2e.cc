#include <iostream>
#include <complex>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include "geometry.hh"
#include "gamma.hh"
#include "momf2e.hh"

using namespace boost::numeric::ublas;
using std::complex;



void moment_f2e(double cor[],
		vector< complex<double> > &u,
		vector< complex<double> > &v,
		matrix< complex<double> > &config,
		const int T, const int L) {


  const int svol = L*L*L;
  const int volume = svol*T;
  vector< complex<double> > eta (volume*12);
  vector< complex<double> > chi (volume*12);

  long double C[T];
  for(size_t t = 0; t < T; t++) {
    C[t] = 0.;
  }
  vector< complex<double> > vtmp(12);
  for(size_t mu = 0; mu < 4; mu++) {
    // now we compute chi_mu(x+mu) = gamma_mu U^+_x,mu psi(x)
    // and eta_mu(x) = gamma_mu U_x,mu psi(x+mu)
    // and contract with bar(psi) for each t
    complex<double> phase(1,0);
    if(mu == 0) {
      // antiperiodic boundary conditions in time for the fermion fields
      // need to be compensated for...
      phase = std::exp(complex<double>(0,1)*3.141593/T);
    }
    int dir[4]={0,0,0,0};
    dir[mu]=1;
    double res_eta[T], res_chi[T];
    for(size_t t = 0; t < T; t++) {
      res_eta[t] = 0.;
      res_chi[t] = 0.;
      for(size_t x = 0; x < L; x++) {
        for(size_t y = 0; y < L; y++) {
          for(size_t z = 0; z < L; z++) {
            size_t ix = get_index(t, x, y, z, T, L);
	    size_t ixu = get_index(t+dir[0], x+dir[1], y+dir[2], z+dir[3], T, L);
	    matrix_range<matrix< complex<double> > > U (config, range(ggi(ix, mu), ggi(ix,mu)+3), range(0,3));

	    size_t k = gsi(ix);
	    size_t ku = gsi(ixu);

	    for(size_t i = 0; i < 12; i+=3) {
              size_t j = gamma_permutation[mu][i/3];
	      noalias(project(eta, range(k+i, k+i+3))) = phase*gamma_factor[5][i/3]*gamma_factor[mu][j]*prod(U, project(v, range(ku+j*3, ku+j*3+3)));
	    }

	    for(size_t i = 0; i < 12; i+=3) {
              size_t j = gamma_permutation[mu][i/3];
              noalias(project(chi, range(ku+i, ku+i+3))) = conj(phase)*gamma_factor[5][i/3]*gamma_factor[mu][j]*prod(herm(U), project(v, range(k+j*3, k+j*3+3)));
	    }
          }
        }
      }
    }

    for(size_t t = 0; t < T; t++) {
      vector_range< vector< complex<double> > > vr(u, range(svol*12*t, svol*12*(t+1)));
      vector_range< vector< complex<double> > > etavr(eta, range(svol*12*t, svol*12*(t+1)));
      vector_range< vector< complex<double> > > chivr(chi, range(svol*12*t, svol*12*(t+1)));
      
      res_eta[t] = real(inner_prod(conj(vr), etavr));
      res_chi[t] = real(inner_prod(conj(vr), chivr));
      if(mu != 0) C[t] -= 1./3.*(res_chi[t] - res_eta[t]);
    }
    if(mu == 0) {
      for(size_t t = 0; t < T; t++) {
        size_t tu = (t+1+T)%T;
        size_t td = (t-1+T)%T;
        C[t] += 0.5*(res_chi[tu] - res_eta[td] - res_eta[t] + res_chi[t]);
      }
    }
  }
  for(size_t t = 0; t < T; t++) {
    cor[t] = double(C[t])/double(svol);
  }
  return;
}
