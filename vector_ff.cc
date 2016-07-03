#include <iostream>
#include <complex>
#include <fstream>
#include <iterator>
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include "geometry.hh"
#include "local_vec.hh"
#include "scalar.hh"
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
using std::ifstream;
using std::string;
using std::ios;

compressed_matrix<complex<double> > gamma5 (12, 12, 12);
compressed_matrix<complex<double> > gamma0 (12, 12, 12);
compressed_matrix<complex<double> > gamma1 (12, 12, 12);
compressed_matrix<complex<double> > gamma2 (12, 12, 12);
compressed_matrix<complex<double> > gamma3 (12, 12, 12);



int main (int ac, char* av[]) {

  int T, L, t0 = -1, nstore;
  string propfilename, gpropfilename, bpropfilename;
  int proppos, gproppos, bproppos, momentum, samples;
  bool binary = false, sampleout = false;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce this help message")
      ("spatialsize,L", po::value<int>(&L), "spatial lattice size")
      ("temporalsize,T", po::value<int>(&T), "temporal lattice size")
      ("config-number,n", po::value<int>(&nstore)->default_value(0), "configuration number")
      ("propagator-filename", po::value< string >(&propfilename), "basefile name of the forward propagator")
      ("b-propagator-filename", po::value< string >(&bpropfilename), "basefile name of the second forward propagator")
      ("gen-propagator-filename", po::value< string >(&gpropfilename), "basefile name of the sequential propagator")
      ("propagator-pos", po::value< int >(&proppos)->default_value(0), "position of forward propagator binary data in lime file")
      ("gen-propagator-pos", po::value< int >(&gproppos)->default_value(0), "position of sequential propagator binary data in lime file")
      ("b-propagator-pos", po::value< int >(&bproppos)->default_value(0), "position of second forward propagator binary data in lime file")
      ("momentum,p",  po::value< int >(&momentum)->default_value(0), "momentum counter for naming output files")
      ("binary,b", "write results in binary format instead of text format")
      ("sample-output,S", "write results additionally sample-wise in binary format") 
      ("samples,s", po::value< int >(&samples)->default_value(1), "number of samples per gauge")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      cout << desc << endl;
      return 1;
    }
    if (vm.count("sample-output") || (vm.count("S"))) {
      sampleout = true;
    }
    if (vm.count("binary") || (vm.count("b"))) {
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
      cerr << "file name of the forward propagator must be given!" << endl;
      cout << endl << desc << endl;
      return 1;
    }
    if (!vm.count("b-propagator-filename") && !vm.count("help")) {
      cerr << "file name of the second forward propagator must be given!" << endl;
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
  vector< complex<double> > bv (volume*12);
  define_etmc_gammas(gamma0, gamma1, gamma2, gamma3, gamma5);

  double *ppcor = new double[T];
  double *cor = new double[T];
  complex<double> *scor = new complex<double>[T];
  double *tmp = new double[T];
  double *ptmp = new double[T];
  complex<double> *stmp = new complex<double>[T];
  // set correlators to zero
  for(int t = 0; t < T; t++) {
    ppcor[t] = 0.;
    cor[t] = 0.;
    scor[t] = 0.;
  }

  ofstream opp, ovec, osca;
  if(sampleout) {
    char fn[400];
    sprintf(fn, "ppcor.samples.%.2d.%.4d", momentum, nstore);
    opp.open(fn, ios::out|ios::binary);
    sprintf(fn, "vector_ff.samples.%.2d.%.4d", momentum, nstore);
    ovec.open(fn, ios::out|ios::binary);
    sprintf(fn, "scalar_ff.samples.%.2d.%.4d", momentum, nstore);
    osca.open(fn, ios::out|ios::binary);
  }

  for(size_t s = 0; s < samples; s++) {

    ifstream ifs;
    char source_filename[400];
    
    // we first need to determine the time slice for this sample
    for(size_t t = 0; t < T; t++) {
      sprintf(source_filename, "%s.%.4d.%.5lu.%.2lu.inverted", propfilename.c_str(), nstore, s, t);
      ifs.open(source_filename);
      if(ifs.is_open()) {
        ifs.close();
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
    
    sprintf(fn, "%s.%.4d.%.5lu.%.2d.inverted", gpropfilename.c_str(), nstore, s, t0);
    cout << "Reading gen. propagator from file " << fn << endl;
    if(read_lime_spinor(gv, fn, gproppos, -1, T, L, L, L) != 0) {
      cout << "Could not read spinor from file, aborting...!" << endl;
      exit(-1);
    }
    
    cout << "... reading for sample " << s << " done...!" << endl;
    
    cout << "Computing 2-pt function " << endl;
    for(int tt = 0; tt < T; tt++) {
      int t = (tt + t0)%T;
      vector_range< vector< complex<double> > > vr(v, range(svol*12*t, svol*12*(t+1)));
      vector_range< vector< complex<double> > > bvr(bv, range(svol*12*t, svol*12*(t+1)));
      ptmp[tt] = real(inner_prod(conj(bvr), vr));
      ppcor[tt] += ptmp[tt];
    }
    // sample wise output
    if(sampleout && opp.is_open()) {
      opp.write((char*) ptmp, (size_t) T*sizeof(double));
    }
    
    cout << "Computing 3-pt function for local vector V_0" << endl;
    local_vec(tmp, v, gv, (long unsigned int) T, (long unsigned int) L);
    for(int tt = 0; tt < T; tt++) {
      int t = (tt - t0+T)%T;
      cor[t] += tmp[tt];
    }
    // sample wise output
    if(sampleout && ovec.is_open()) {
      for(int tt = 0; tt < T; tt++) {
        int t = (tt + t0)%T;
        ovec.write((char*) &tmp[t], (size_t) sizeof(double));
      }
    }

    cout << "Computing 3-pt function for scalar" << endl;
    scalar(stmp, v, gv, (long unsigned int) T, (long unsigned int) L);
    for(int tt = 0; tt < T; tt++) {
      int t = (tt - t0+T)%T;
      scor[t] += stmp[tt];
    }
    // sample wise output
    if(sampleout && ovec.is_open()) {
      for(int tt = 0; tt < T; tt++) {
        int t = (tt + t0)%T;
        osca.write((char*) &stmp[t], (size_t) sizeof(complex<double>));
      }
    }
  }

  if(sampleout) {
    opp.close();
    ovec.close();
    osca.close();
  }
  for(int t = 0; t < T; t++) {
    ppcor[t] = ppcor[t]/svol/double(samples);
    cor[t]   = cor[t]/double(samples);
    scor[t]  = scor[t]/double(samples);
  }


  // now the output
  char fn[400];
  sprintf(fn, "ppcor.%.2d.%.4d", momentum, nstore);
  ofstream ofs(fn);
  if(binary) {
    ofs.write((char*) ppcor, T*sizeof(double));
  }
  else {
    for(int t = 0; t < T; t++) {
      ofs << t << " " << ppcor[t] << endl;
    }
  }
  ofs.close();
  
  sprintf(fn, "vector_ff.%.2d.%.4d", momentum, nstore);
  ofs.open(fn);
  if(binary) {
    ofs.write((char*) cor, T*sizeof(double));
  }
  else {
    for(int t = 0; t < T; t++) {
      ofs << t << " " << cor[t] << endl;
    }
  }
  ofs.close();

  sprintf(fn, "scalar_ff.%.2d.%.4d", momentum, nstore);
  ofs.open(fn);
  if(binary) {
    ofs.write((char*) scor, T*sizeof(complex<double>));
  }
  else {
    for(int t = 0; t < T; t++) {
      ofs << t << " " << real(scor[t]) << " " << imag(scor[t]) << endl;
    }
  }
  ofs.close();
  
  return(0);
}
