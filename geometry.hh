#ifndef _GEOMETRY_HH
#define _GEOMETRY_HH

inline int get_index(const int t, const int x, const int y, const int z,
	      const int T, const int L) {
  int tt = (t+T)%T;
  int xx = (x+L)%L;
  int yy = (y+L)%L;
  int zz = (z+L)%L;
  return(((tt*L+xx)*L+yy)*L+zz);
}

inline int ggi(const int ix, const int mu) {
  
  return((4*ix+mu)*3);
}

inline int gsi(const int ix) {
  return(ix*12);
}

#endif
