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
#include "plaquette.hh"

using namespace boost::numeric::ublas;
using std::complex;

double plaquette(matrix< complex<double> > &config,
		 const int T, const int L) {
  double res=0.;
  matrix< complex<double> > tmp1(3,3);
  matrix< complex<double> > tmp2(3,3);
  matrix< complex<double> > tmp3(3,3);
  for(int t = 0; t < T; t++) {
    for(int x = 0; x < L; x++) {
      for(int y = 0; y < L; y++) {
	for(int z = 0; z < L; z++) {
	  int ix = get_index(t, x, y, z, T, L);
	  for(int mu1 = 0; mu1 < 3; mu1++) {
	    int dir1[4]={0,0,0,0};
	    dir1[mu1]=1;
	    int ix1 = get_index(t+dir1[0], x+dir1[1], y+dir1[2], z+dir1[3], T, L); 
	    matrix_range<matrix<complex<double> > > c1(config, range(ggi(ix,mu1), ggi(ix, mu1)+3), range(0,3));
	    for(int mu2 = mu1+1; mu2 < 4; mu2++) {
	      int dir2[4]={0,0,0,0};
	      dir2[mu2]=1;
	      int ix2 = get_index(t+dir2[0], x+dir2[1], y+dir2[2], z+dir2[3], T, L); 
	      matrix_range<matrix<complex<double> > > c2(config, range(ggi(ix1,mu2), ggi(ix1, mu2)+3), range(0,3));
	      matrix_range<matrix<complex<double> > > c3(config, range(ggi(ix,mu2), ggi(ix,mu2)+3), range(0,3));
	      matrix_range<matrix<complex<double> > > c4(config, range(ggi(ix2,mu1), ggi(ix2,mu1)+3), range(0,3));
	      noalias(tmp2) = herm(prod(c3,c4));
	      noalias(tmp1) = prod(c1,c2);
	      res += real(prod(tmp1, tmp2)(0,0)) + real(prod(tmp1, tmp2)(1,1)) + real(prod(tmp1, tmp2)(2,2));
	    }
	  }
	}
      }
    }
  }
  return(res/18./L/L/L/T);
}
