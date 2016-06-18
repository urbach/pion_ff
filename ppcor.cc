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

int main (int ac, char* av[]) {

  int T, L, t0 = -1, nstore;
  string propfilename, configfilename, bpropfilename;
  double kappa;
  int proppos, bproppos, momentum, samples;
  bool binary = false;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce this help message")
      ("spatialsize,L", po::value<int>(&L), "spatial lattice size")
      ("temporalsize,T", po::value<int>(&T), "temporal lattice size")
      ("config-number,n", po::value<int>(&nstore)->default_value(0), "configuration number")
      ("kappa,k", po::value<double>(&kappa)->default_value(0.5), "hopping parameter value")
      ("propagator-filename", po::value< string >(&propfilename), "basefile name of the propagator")
      ("b-propagator-filename", po::value< string >(&bpropfilename), "basefile name of the b-propagator")
      ("config-filename", po::value< string >(&configfilename), "configuration filename")
      ("propagator-pos", po::value< int >(&proppos)->default_value(0), "position of propagator binary data in lime file")
      ("b-propagator-pos", po::value< int >(&bproppos)->default_value(0), "position of propagator binary data in lime file")
      ("momentum,p",  po::value< int >(&momentum)->default_value(0), "momentum counter for naming output files")
      ("binary,b", "write results in binary format instead of text format")
      ("samples,s", po::value< int >(&samples)->default_value(1), "number of samples per gauge")
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
    if (!vm.count("b-propagator-filename") && !vm.count("help")) {
      cerr << "file name of the b-propagator must be given!" << endl;
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
  vector< complex<double> > bv (volume*12);

  double *ppcor = new double[T];
  // set correlators to zero
  for(int t = 0; t < T; t++) {
    ppcor[t] = 0.;
  }

  for(size_t s = 0; s < samples; s++) {

    FILE * ifs = NULL;
    char source_filename[400];
    
    // we first need to determine the time slice for this sample
    for(size_t t = 0; t < T; t++) {
      sprintf(source_filename, "%s.%.4d.%.5lu.%.2lu.inverted", propfilename.c_str(), nstore, s, t);
      if( (ifs = fopen(source_filename, "r")) != NULL) {
        fclose(ifs);
        t0 = t;
        break;
      }
    }
    if(t0 == -1) {
      cerr << "could not determine time slice, Aborting...!" << endl;
      exit(0);
    }

    char fn[400];
    sprintf(fn, "%s.%.4d.%.5lu.%.2d.inverted", propfilename.c_str(), nstore, s, t0);
    cout << "Reading propagator from file " << fn << endl;
    if(read_lime_spinor(v, fn, proppos, -1, T, L, L, L) != 0) {
      cout << "Could not read spinor from file, aborting...!" << endl;
      exit(-1);
    }
    
    sprintf(fn, "%s.%.4d.%.5lu.%.2d.inverted", bpropfilename.c_str(), nstore, s, t0);
    cout << "Reading b-propagator from file " << fn << endl;
    if(read_lime_spinor(bv, fn, bproppos, -1, T, L, L, L) != 0) {
      cout << "Could not read spinor from file, aborting...!" << endl;
      exit(-1);
    }
    
    cout << "... reading for sample " << s << " done...!" << endl;
    
    cout << "Computing 2-pt function " << endl;
    for(int tt = 0; tt < T; tt++) {
      int t = (tt + t0)%T;
      //int t = tt;
      vector_range< vector< complex<double> > > vr(v, range(svol*12*t, svol*12*(t+1)));
      vector_range< vector< complex<double> > > bvr(bv, range(svol*12*t, svol*12*(t+1)));
      ppcor[tt] += real(inner_prod(conj(bvr), vr));
    }
  }

  // now the output

  char prev;
  int prev2;
  ostringstream oss;
  oss << "ppcor.";
  prev2 = oss.width(2);
  prev = oss.fill('0'); 
  oss << momentum;
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
  for(int t = 0; t < T; t++) {
    ofs << t << " " << ppcor[t]/svol*2*2*kappa*kappa/double(samples) << endl;
  }
  ofs.close();
  
  return(0);
}
