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
      
      ;

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

  ostringstream ossvv, ossv4;
  ossvv << "disc.vv." << nstore << ends;
  ossv4 << "disc.v4." << nstore << ends;
  ofstream ofs;
  ofs.open(ossvv.str().c_str());
  ofstream ofsv4;
  ofsv4.open(ossv4.str().c_str());
  ofs << "# k=" << kappa << " mu=" << mu 
      << " L=" << L << " T=" << T << endl; 
  ofsv4 << "# k=" << kappa << " mu=" << mu 
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
    cout << "Computing disconnected contribution for sample " << sample+1 << endl;

    for(int t=0; t < T; t++) {
      for(int i = 0; i < 16; i++) {
	compressed_matrix< complex<double> > g(12, 12, 12);
	g = prod(gamma5, project(gamma, range(i*12,(i+1)*12), range(0,12)));
	matrix_range< compressed_matrix<complex<double> > > 
	  g2(gamma, range(i*12,(i+1)*12), range(0,12));
	complex<double> tmp = 0.;
	complex<double> tmp2 = 0.;
	for(int j = 0; j < svol; j++) {
	  vector_range< vector< complex<double> > > vr(v, range(svol*12*t+j*12, svol*12*(t+1)+(j+1)*12));
	  tmp += inner_prod(conj(vr), prod(g, vr));
	  vector_range< vector< complex<double> > > ur(u, range(svol*12*t+j*12, svol*12*(t+1)+(j+1)*12));
	  tmp2 += inner_prod(conj(ur), prod(g2, vr));
	}
	tmp *= mu*kappa;
	prev2 = ofs.width();
	ofs << setw(5) << nstore << setw(3) << i+1 << setw(6) << sample+1 << setw(3) << t+1;
	ofs.setf(std::ios::fixed);
	ofs.setf(std::ios::showpoint);
	ofs << setw(15) << setprecision(5) << imag(tmp);
	ofs << setw(15) << setprecision(5) << real(tmp) ;
	ofs << endl;
	ofs.unsetf(std::ios::fixed);
	ofs.unsetf(std::ios::showpoint);
	ofs.width(prev2);
	
	if(i == 0) {
	  tmp2 += addimag;
	}
	if(i == 8) {
	  tmp2 += addreal;
	}

	tmp2 /= 2.;
	prev2 = ofsv4.width();
	ofsv4 << setw(5) << nstore << setw(3) << i+1 << setw(6) << sample+1 << setw(3) << t+1;
	ofsv4.setf(std::ios::fixed);
	ofsv4.setf(std::ios::showpoint);
	ofsv4 << setw(15) << setprecision(5) << real(tmp2);
	ofsv4 << setw(15) << setprecision(5) << imag(tmp2) ;
	ofsv4 << endl;
	ofsv4.unsetf(std::ios::fixed);
	ofsv4.unsetf(std::ios::showpoint);
	ofsv4.width(prev2);
      }
    }
    ofs.close();
  }
  return(0);
}
