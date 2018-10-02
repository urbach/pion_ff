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
#include "momf2e.hh"
#include "gamma.hh"
#include "plaquette.hh"
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

compressed_matrix<complex<double> > gamma5 (12, 12, 12);
compressed_matrix<complex<double> > gamma0 (12, 12, 12);
compressed_matrix<complex<double> > gamma1 (12, 12, 12);
compressed_matrix<complex<double> > gamma2 (12, 12, 12);
compressed_matrix<complex<double> > gamma3 (12, 12, 12);



int main (int ac, char* av[]) {

  int T, L, t0, nstore;
  string propfilename, gpropfilename, configfilename;
  double kappa;
  int proppos, gproppos;
  bool binary = false;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce this help message")
      ("spatialsize,L", po::value<int>(&L), "spatial lattice size")
      ("temporalsize,T", po::value<int>(&T), "temporal lattice size")
      ("sourcetime,t", po::value<int>(&t0)->default_value(0), "source time slice")
      ("config-number,n", po::value<int>(&nstore)->default_value(0), "configuration number")
      ("kappa,k", po::value<double>(&kappa)->default_value(0.5), "hopping parameter value")
      ("propagator-filename", po::value< string >(&propfilename), "file name of the propagator")
      ("gen-propagator-filename", po::value< string >(&gpropfilename), "file name of the generalised propagator")
      ("config-filename", po::value< string >(&configfilename), "base filename of the configuration")
      ("propagator-pos", po::value< int >(&proppos)->default_value(0), "position of propagator binary data in lime file")
      ("gen-propagator-pos", po::value< int >(&gproppos)->default_value(0), "position of gen-propagator binary data in lime file")
      ("binary,b", "write results in binary format instead of text format")      
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      cout << desc << endl;
      return 1;
    }
    if (vm.count("binary")) {
      binary=true;
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
    if (!vm.count("propagator-filename") && !vm.count("help")) {
      cerr << "file name of the propagator must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("gen-propagator-filename") && !vm.count("help")) {
      gpropfilename = propfilename;
      if (!vm.count("gen-propagator-pos")) {
	gproppos = 1;
      }
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
  vector< complex<double> > gv (volume*12);
  define_etmc_gammas(gamma0, gamma1, gamma2, gamma3, gamma5);

  char fn[100];
  strcpy(fn, propfilename.c_str());
  cout << "Reading propagator from file " << propfilename << endl;
  //read_cmi(v, L, T, propfilename.c_str());
  if(read_lime_spinor(v, fn, proppos, -1, T, L, L, L) != 0) {
    cout << "Could not read spinor from file, aborting...!" << endl;
    exit(-1);
  }
  cout << "Reading gen. propagator from file " << gpropfilename << endl;
  //read_cmi(gv, L, T, gpropfilename.c_str());
  strcpy(fn, gpropfilename.c_str());
  if(read_lime_spinor(gv, fn, gproppos, -1, T, L, L, L) != 0) {
    cout << "Could not read spinor from file, aborting...!" << endl;
    exit(-1);
  }

  char prev;
  int prev2;

  cout << "Computing 2-pt function " << endl;
  ostringstream oss;
  oss << "ppcor.";
  prev2 = oss.width(2);
  prev = oss.fill('0'); 
  oss << t0;
  oss.fill(prev);
  oss.width(prev2);
  oss << ".";
  oss.width(4);
  oss.fill('0'); 
  oss << nstore;
  oss.fill(prev);
  oss.width(prev2);
  oss << ends;
  ofstream ofs(oss.str().c_str());
  for(int tt = 0; tt < T; tt++) {
    int t = (tt + t0)%T;
    //int t = tt;
    vector_range< vector< complex<double> > > vr(v, range(svol*12*t, svol*12*(t+1)));
    ofs << tt << " " << real(inner_prod(conj(vr), vr))/svol*2*2*kappa*kappa << endl;
  }
  ofs.close();

  matrix< complex<double> > config (4*3*volume, 3);
  cout << "Reading configuration from file " << configfilename << endl;
  read_lime_gauge_field_doubleprec(config, configfilename.c_str(), T, L, L, L);
  cout << "Computing plaquette..." << endl;
  cout << "plaquette = " << plaquette(config, T, L) << endl;

  double cor[T];
  cout << "Computing 3-pt function for O_44" << endl;
  moment_f2e(cor, v, gv, config, T, L);
  ostringstream oss2;
  oss2 << "momf2e.";
  prev2 = oss2.width(2);
  prev = oss2.fill('0'); 
  oss2 << t0;
  oss2.fill(prev);
  oss2.width(prev2);
  oss2 << ".";
  oss2.width(4);
  oss2.fill('0'); 
  oss2 << nstore;
  oss2.fill(prev);
  oss2.width(prev2);
  oss2 << ends;
  ofs.open(oss2.str().c_str());
  for(int tt=0; tt < T; tt++) {
    int t = (tt+t0)%T;
    ofs << tt << " " << cor[t]*2*2*2*kappa*kappa*kappa << endl;
  }
  ofs.close();

  return(0);
}
