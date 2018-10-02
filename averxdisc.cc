#include <iostream>
#include <complex>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
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
#include "io.hh"

using namespace boost::numeric::ublas;
namespace po = boost::program_options;
using std::endl;
using std::ends;
using std::cout;
using std::cerr;
using std::complex;
using std::ofstream;
using std::string;
using std::ostringstream;
using std::setw;
using std::setprecision;

compressed_matrix<complex<double> > gamma5 (12, 12, 12);
compressed_matrix<complex<double> > gamma0 (12, 12, 12);
compressed_matrix<complex<double> > gamma1 (12, 12, 12);
compressed_matrix<complex<double> > gamma2 (12, 12, 12);
compressed_matrix<complex<double> > gamma3 (12, 12, 12);

void moment_f2e_disc(complex<double> cor[],
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
    complex<double> result[4] = {0.,0.,0.,0.};
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
	    result[mu] += inner_prod(conj(project(u, range(k, k+12))), prod(g, vtmp));

	    // multiply G(x-mu) with U^\dagger(x-mu,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(herm(Ud), project(v, range(kd+i, kd+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(herm(Ud), project(v, range(kd+i, kd+i+3)));
	    }
	    // -S^\dagger(x) \gamma_mu times vtmp
	    result[mu] -= inner_prod(conj(project(u, range(k, k+12))), prod(g, vtmp));

	    // multiply G(x) with U^\dagger(x,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(herm(U), project(v, range(k+i, k+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(herm(U), project(v, range(k+i, k+i+3)));
	    }
	    // -S^\dagger(x+mu) \gamma_mu times vtmp
	    result[mu] -= inner_prod(conj(project(u, range(ku, ku+12))), prod(g, vtmp));

	    // multiply G(x) with U(x-mu,mu) -> vtmp
	    for(int i = 0; i < 12; i+=3) {
	      //axpy_prod(Ud, project(v, range(k+i, k+i+3)), project(vtmp, range(i, i+3)), true);
	      noalias(project(vtmp, range(i, i+3))) = prod(Ud, project(v, range(k+i, k+i+3)));
	    }
	    // S^\dagger(x-mu) \gamma_mu times vtmp
	    result[mu] += inner_prod(conj(project(u, range(kd, kd+12))), prod(g, vtmp));
	  }
	}
      }
    }
    cor[t] = -0.5*(result[0] - (result[1]+result[2]+result[3])/3.)/double(L)/double(L)/double(L);
  }
  return;
}



int main (int ac, char* av[]) {

  int T, L, nstore, s;
  string propfilename, gpropfilename, configfilename, postfix;
  double kappa, mu;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce this help message")
      ("spatialsize,L", po::value<int>(&L), "spatial lattice size")
      ("temporalsize,T", po::value<int>(&T), "temporal lattice size")
      ("config-number,n", po::value<int>(&nstore)->default_value(0), "configuration number")
      ("mu,m", po::value<double>(&mu)->default_value(0.5), "twisted mass")
      ("kappa,k", po::value<double>(&kappa)->default_value(0.5), "hopping parameter value")
      ("no-samples,s", po::value<int>(&s)->default_value(1), "number of stochastic samples")
      ("propagator-prefix", po::value< string >(&propfilename), "prefix of the stochastic propagator filename")
      ("propagator-postfix", po::value< string >(&postfix)->default_value(".inverted"), "postfix of the stochastic propagator filename")
      ("config-filename,c", po::value< string >(&configfilename)->default_value("conf"), "filename of the configuration")      ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      cout << desc << endl;
      return 1;
    }
    if (!vm.count("spatialsize") && !vm.count("help")) {
      cerr << "spatial lattice size must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("temporalsize") && !vm.count("help")) {
      cerr << "temporal lattice size must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("propagator-prefix") && !vm.count("help")) {
      cerr << "prefix of the propagator filename must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("config-filename") && !vm.count("help")) {
      cerr << "configuration filename must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
  }
  catch(std::exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }
  
  const int svol = L*L*L;
  const int volume = svol*T;
  
  vector< complex<double> > v (volume*12);
  vector< complex<double> > u (volume*12);

  compressed_matrix< complex<double> > gamma(16*12, 12, 16*12);
  //compressed_matrix< complex<double> > gamma2(16*12, 12, 16*12);
  complex<double> cI(0.,1.);
  // define independent gamma matrix combinations
  define_ukqcd_gammas(gamma0, gamma1, gamma2, gamma3, gamma5);
  project(gamma, range(0,12), range(0,12)) = gamma5;
  project(gamma, range(1*12,2*12), range(0,12)) = gamma1;
  project(gamma, range(2*12,3*12), range(0,12)) = gamma2;
  project(gamma, range(3*12,4*12), range(0,12)) = gamma3;
  project(gamma, range(4*12,5*12), range(0,12)) = -cI*prod(gamma0, gamma5);
  project(gamma, range(5*12,6*12), range(0,12)) = -cI*prod(gamma0, gamma1);
  project(gamma, range(6*12,7*12), range(0,12)) = -cI*prod(gamma0, gamma2);
  project(gamma, range(7*12,8*12), range(0,12)) = -cI*prod(gamma0, gamma3);
  project(gamma, range(8*12,9*12), range(0,12)) = prod(gamma5, gamma5);
  project(gamma, range(9*12,10*12), range(0,12)) = -cI*prod(gamma5, gamma1);
  project(gamma, range(10*12,11*12), range(0,12)) = -cI*prod(gamma5, gamma2);
  project(gamma, range(11*12,12*12), range(0,12)) = -cI*prod(gamma5, gamma3);
  project(gamma, range(12*12,13*12), range(0,12)) = gamma0;
  project(gamma, range(13*12,14*12), range(0,12)) = prod(gamma5, project(gamma, range(5*12,6*12), range(0,12)));
  project(gamma, range(14*12,15*12), range(0,12)) = prod(gamma5, project(gamma, range(6*12,7*12), range(0,12)));
  project(gamma, range(15*12,16*12), range(0,12)) = prod(gamma5, project(gamma, range(7*12,8*12), range(0,12)));

  matrix< complex<double> > config (4*3*volume, 3);
  cout << "Reading configuration from file " << configfilename << endl;
  read_lime_gauge_field_doubleprec(config, configfilename.c_str(), T, L, L, L);

  ostringstream ossvv;
  ossvv << "momf2edisc.v4." << nstore << ends;
  ofstream ofs;
  ofs.open(ossvv.str().c_str());
  ofs << "# k=" << kappa << " mu=" << mu 
      << " L=" << L << " T=" << T << endl; 

  //nxyz*nc*ndim*2.0
  complex<double> addimag(0.,2*kappa*mu/sqrt(1 + 4*kappa*kappa*mu*mu)*L*L*L*3*4*2.);
  complex<double> addreal(1./sqrt(1 + 4*kappa*kappa*mu*mu)*L*L*L*3*4*2.);

  for(int sample = 0; sample < s; sample++) {
    ostringstream oss, oss3;
    char prev;
    int prev2;
    oss << propfilename.c_str() << ".";
    oss3 << propfilename.c_str() << ".";
    prev2 = oss.width(5);
    prev = oss.fill('0');
    oss3.width(5);
    oss3.fill('0');
    oss << sample+1;
    oss3 << sample+1;
    oss.fill(prev);
    oss.width(prev2);
    oss3.fill(prev);
    oss3.width(prev2);
    oss << postfix.c_str() << ends;
    oss3 << ".applied" << ends;
    define_etmc_gammas(gamma0, gamma1, gamma2, gamma3, gamma5);
    cout << "Reading stochastic propagator from file " << oss.str().c_str() << endl;
    read_cmi(v, L, T, oss.str().c_str());
    cout << "Reading D^4 applied source from file " << oss3.str().c_str() << endl;
    read_cmi(u, L, T, oss3.str().c_str());
    cout << "Rotating vectors to UKQCD gamma matrix convention" << endl;
    rotate_etmc_ukqcd(v, v, gamma0, gamma5, volume);
    rotate_etmc_ukqcd(u, u, gamma0, gamma5, volume);
    define_ukqcd_gammas(gamma0, gamma1, gamma2, gamma3, gamma5);
    cout << "Computing disconnected contribution of momf2e for sample " << sample+1 << endl;

    complex<double> cor[100];

    moment_f2e_disc(cor, u, v, config, T, L);

    for(int t = 0; t < T; t++) {
      ofs << setw(5) << nstore << setw(3) << 1 << setw(6) << sample+1 << setw(3) << t+1;
      ofs.setf(std::ios::fixed);
      ofs.setf(std::ios::showpoint);
      ofs << setw(15) << setprecision(5) << real(cor[t]/2.);
      ofs << setw(15) << setprecision(5) << imag(cor[t]/2.) ;
      ofs << endl;
      ofs.unsetf(std::ios::fixed);
      ofs.unsetf(std::ios::showpoint);
      ofs.width(prev2);
    }

  }
  ofs.close();
  return(0);
}
