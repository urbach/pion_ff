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
#include "DmuDnu.hh"

using namespace boost::numeric::ublas;
using std::complex;


void DmuDnu(double cor[],
	    vector< complex<double> > &u,
	    vector< complex<double> > &v,
	    matrix< complex<double> > &config,
	    const int T, const int L, size_t kappa = 0, 
	    size_t mu = 0, size_t nu = 0) {


  const int svol = L*L*L;
  const int volume = svol*T;
  vector< complex<double> > eta (volume*12);
  vector< complex<double> > chi (volume*12);

  long double C[T];
  for(size_t t = 0; t < T; t++) {
    C[t] = 0.;
  }

  vector< complex<double> > vtmp(12);

  // now we compute chi(x+nu) = U^+_x,nu psi(x)
  // and eta(x) = U_x,nu psi(x+nu)

  complex<double> phase(1,0);
  if(nu == 0) {
    // antiperiodic boundary conditions in time for the fermion fields
    // need to be compensated for...
    phase = std::exp(complex<double>(0,1)*3.141593/T);
  }
  int dir[4]={0,0,0,0};
  dir[nu] = 1;
  double res_eta1[T], res_chi1[T], res_eta2[T], res_chi2[T];
  for(size_t t = 0; t < T; t++) {
    for(size_t x = 0; x < L; x++) {
      for(size_t y = 0; y < L; y++) {
	for(size_t z = 0; z < L; z++) {
	  size_t ix = get_index(t, x, y, z, T, L);
	  size_t ixu = get_index(t+dir[0], x+dir[1], y+dir[2], z+dir[3], T, L);
	  matrix_range<matrix< complex<double> > > U (config, range(ggi(ix, nu), ggi(ix,nu)+3), range(0,3));
	  
	  size_t k = gsi(ix);
	  size_t ku = gsi(ixu);
	  
	  for(size_t i = 0; i < 12; i+=3) {
	    noalias(project(eta, range(k+i, k+i+3))) = phase*prod(U, project(v, range(ku+i, ku+i+3)));
	  }
	  
	  for(size_t i = 0; i < 12; i+=3) {
	    noalias(project(chi, range(ku+i, ku+i+3))) = conj(phase)*prod(herm(U), project(v, range(k+i, k+i+3)));
	  }
	}
      }
    }
  }
  dir[nu] = 0;
  dir[mu] = 1;

  if(mu == 0) {
    // antiperiodic boundary conditions in time for the fermion fields
    // need to be compensated for...
    phase = std::exp(complex<double>(0,1)*3.141593/T);
  }

  vector< complex<double> > eta1 (volume*12);
  vector< complex<double> > eta2 (volume*12);
  vector< complex<double> > chi1 (volume*12);
  vector< complex<double> > chi2 (volume*12);

  // and eta1(x) = U_x,mu eta(x+mu)
  // and chi1(x) = U_x,mu chi(x+mu)
  // and eta2(x+mu) = U^+_x,mu eta(x)
  // and chi2(x+mu) = U^+_x,mu chi(x)

  for(size_t t = 0; t < T; t++) {
    for(size_t x = 0; x < L; x++) {
      for(size_t y = 0; y < L; y++) {
	for(size_t z = 0; z < L; z++) {
	  size_t ix = get_index(t, x, y, z, T, L);
	  size_t ixu = get_index(t+dir[0], x+dir[1], y+dir[2], z+dir[3], T, L);
	  size_t ixd = get_index(t-dir[0], x-dir[1], y-dir[2], z-dir[3], T, L);
	  matrix_range<matrix< complex<double> > > U (config, range(ggi(ix, mu), ggi(ix, mu)+3), range(0,3));
	  
	  size_t k = gsi(ix);
	  size_t kd = gsi(ixd);
	  size_t ku = gsi(ixu);
	  
	  for(size_t i = 0; i < 12; i+=3) {
	    noalias(project(eta1, range(k+i, k+i+3))) = phase*prod(U, project(eta, range(ku+i, ku+i+3)));
	    noalias(project(chi1, range(k+i, k+i+3))) = phase*prod(U, project(chi, range(ku+i, ku+i+3)));
	  }
	  
	  for(size_t i = 0; i < 12; i+=3) {
	    noalias(project(eta2, range(ku+i, ku+i+3))) = conj(phase)*prod(herm(U), project(eta, range(k+i, k+i+3)));
	    noalias(project(chi2, range(ku+i, ku+i+3))) = conj(phase)*prod(herm(U), project(chi, range(k+i, k+i+3)));
	  }

	  // this is tested to be the same
	  // 	  matrix_range<matrix< complex<double> > > Ud (config, range(ggi(ixd, mu), ggi(ixd, mu)+3), range(0,3));
	  // 	  for(size_t i = 0; i < 12; i+=3) {
	  // 	    noalias(project(eta2, range(k+i, k+i+3))) = conj(phase)*prod(herm(Ud), project(eta, range(kd+i, kd+i+3)));
	  // 	    noalias(project(chi2, range(k+i, k+i+3))) = conj(phase)*prod(herm(Ud), project(chi, range(kd+i, kd+i+3)));
	  // 	  }
	}
      }
    }
  }

  for(size_t t = 0; t < T; t++) {
    vector_range< vector< complex<double> > > vr(u, range(svol*12*t, svol*12*(t+1)));
    vector_range< vector< complex<double> > > eta1vr(eta1, range(svol*12*t, svol*12*(t+1)));
    vector_range< vector< complex<double> > > chi1vr(chi1, range(svol*12*t, svol*12*(t+1)));
    vector_range< vector< complex<double> > > eta2vr(eta2, range(svol*12*t, svol*12*(t+1)));
    vector_range< vector< complex<double> > > chi2vr(chi2, range(svol*12*t, svol*12*(t+1)));
    
    res_eta1[t] = real(inner_prod(conj(vr), eta1vr));
    res_chi1[t] = real(inner_prod(conj(vr), chi1vr));
    res_eta2[t] = real(inner_prod(conj(vr), eta2vr));
    res_chi2[t] = real(inner_prod(conj(vr), chi2vr));
    C[t] = 1./4.*(res_eta1[t] - res_chi1[t] - res_eta2[t] + res_chi2[t]);
    std::cout << C[t] << std::endl;
  }

  for(size_t t = 0; t < T; t++) {
    cor[t] = double(C[t])/double(svol);
  }

  return;
}
