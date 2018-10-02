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
#include "momf2edisc.hh"

using namespace boost::numeric::ublas;
using std::complex;


extern compressed_matrix<complex<double> > gamma5;
extern compressed_matrix<complex<double> > gamma0;
extern compressed_matrix<complex<double> > gamma1;
extern compressed_matrix<complex<double> > gamma2;
extern compressed_matrix<complex<double> > gamma3;


void moment_f2e_disc(double cor[],
		     vector< complex<double> > &u,
		     vector< complex<double> > &v,
		     matrix< complex<double> > &config,
		     const int T, const int L) {

  compressed_matrix< complex<double> > gamma(4*12, 12, 4*12);
  project(gamma, range(0,12), range(0,12)) = gamma0;
  project(gamma, range(12,24), range(0,12)) = gamma1;
  project(gamma, range(24,36), range(0,12)) = gamma2;
  project(gamma, range(36,48), range(0,12)) = gamma3;

//   for(int mu=0; mu < 4; mu++) {
//     cout << project(gamma, range(mu*12,mu*12+12), range(0,12)) << endl << endl;
//   }

  vector< complex<double> > vtmp(12);
  for(int t = 0; t < T; t++) {
    double result[4] = {0.,0.,0.,0.};
    for(int x = 0; x < L; x++) {
      for(int y = 0; y < L; y++) {
	for(int z = 0; z < L; z++) {
	  int ix = get_index(t, x, y, z, T, L);
	  for(int mu = 0; mu < 4; mu++) {
	    int dir[4]={0,0,0,0};
	    dir[mu]=1;
	    int ixd = get_index(t-dir[0], x-dir[1], y-dir[2], z-dir[3], T, L);
	    int ixu = get_index(t+dir[0], x+dir[1], y+dir[2], z+dir[3], T, L);
	    //	    cout <<t<<" "<<x<<" "<<y<<" "<<z<<" "<<mu<<" "<<ix<<" "<<ixu<<" "<<ixd<<endl;
	    matrix_range<matrix< complex<double> > > U (config, range(ggi(ix, mu), ggi(ix,mu)+3), range(0,3));
	    matrix_range<matrix< complex<double> > > Ud (config, range(ggi(ixd,mu), ggi(ixd,mu)+3), range(0,3));
	    matrix_range<compressed_matrix< complex<double> > > g(gamma, range(mu*12, (mu+1)*12), range(0,12));
	    int k = gsi(ix);
	    int kd = gsi(ixd);
	    int ku = gsi(ixu);
	    // multiply G(x+mu) with U(x,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //vector_range< complex<double> > vr(v, range(ku+i, ku+i+3));
	      //vector_range< complex<double> > vtmpr(vtmp, range(i, i+3));
	      //axpy_prod(U, vr, vtmpr, true);
	      noalias(project(vtmp, range(i, i+3))) = prod(U, project(v, range(ku+i, ku+i+3)));
	    }
	    // S^\dagger(x) \gamma_mu times vtmp
	    result[mu] += real(inner_prod(conj(project(u, range(k, k+12))), prod(g, vtmp)));

	    // multiply G(x-mu) with U^\dagger(x-mu,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(herm(Ud), project(v, range(kd+i, kd+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(herm(Ud), project(v, range(kd+i, kd+i+3)));
	    }
	    // -S^\dagger(x) \gamma_mu times vtmp
	    result[mu] -= real(inner_prod(conj(project(u, range(k, k+12))), prod(g, vtmp)));

	    // multiply G(x) with U^\dagger(x,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(herm(U), project(v, range(k+i, k+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(herm(U), project(v, range(k+i, k+i+3)));
	    }
	    // -S^\dagger(x+mu) \gamma_mu times vtmp
	    result[mu] -= real(inner_prod(conj(project(u, range(ku, ku+12))), prod(g, vtmp)));

	    // multiply G(x) with U(x-mu,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(Ud, project(v, range(k+i, k+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(Ud, project(v, range(k+i, k+i+3)));
	    }
	    // S^\dagger(x-mu) \gamma_mu times vtmp
	    result[mu] += real(inner_prod(conj(project(u, range(kd, kd+12))), prod(g, vtmp)));
	  }
	}
      }
    }
    cor[t] = -0.5*(result[0] - (result[1]+result[2]+result[3])/3.)/L/L/L;
  }
  return;
}
