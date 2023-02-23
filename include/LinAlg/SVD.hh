//
// Created by zsoldos on 7/30/20.
//

#ifndef SND_INCLUDE_SND_SVD_HH_
#define SND_INCLUDE_SND_SVD_HH_

#include <limits>
#include <cmath>

#include "LinAlg/Matrix.hh"

// macro-like inline functions

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
{return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
{return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
{return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
{return b < a ? float(b) : (a);}

struct SVD {
  int m, n;
  Matrix u, v;
  DiagMatrix w;
  double eps, tsh;

  SVD(const Matrix &A)
	  : m(A.nrows), n(A.ncols), u(A), v(n,n), w(n) {
	eps = std::numeric_limits<double>::epsilon();
	decompose();
	reorder();
	tsh = 0.5*sqrt(m+n+1)*w[0]*eps;
  };

  void solve(DiagMatrix &b, DiagMatrix &x, double thresh=-1.);
  void solve(Matrix &b, Matrix &x, double thresh=-1.);

  int rank(double thresh=-1.);
  int nullity(double thresh=-1.);
  Matrix range(double thresh=-1.);
  Matrix nullspace(double thresh=-1.);

  double inv_condition() {
	return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
  }

  void decompose();
  void reorder();
  double pythag(const double a, const double b);
};

typedef struct SVD SVD;


#endif // SND_INCLUDE_SND_SVD_HH_
