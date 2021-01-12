//
// Created by zsoldos on 6/17/20.
//

#ifndef _MATHUTILS_HH_
#define _MATHUTILS_HH_

#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>

#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>

#define PI 3.14159265359
#define SQRT5 2.2360679775
#define SQRT2 1.41421356237
#define RINDEX_WATER 1.33
#define SOL_VACUUM 299.792
#define SOL SOL_VACUUM/RINDEX_WATER

static double EvalLL(double nObs, double nPred){
  return nObs*TMath::Log(nPred);
}

static double EvalExtendedLL(double nObs, double nPred){
  return nObs*TMath::Log(nPred) - nPred;
}

static double EvalNLL(double nObs, double nPred){
  double L;
  if(nObs>0 && nPred>0)
	L=nObs*TMath::Log(nObs/nPred) + nPred-nObs;
  else
	L=nPred;
  return -L;
}

template <typename T>
double CalculateLL(T const *hPDF, T const *hExp, bool isNormalized = true){

  auto nBinsX = hPDF->GetNbinsX();
  auto nBinsY = hPDF->GetNbinsY();

  auto normPDF = 1.;
  auto normExp = 1.;
  if(!isNormalized){
	normPDF = hPDF->Integral("width");
	normExp = hExp->Integral("width");
  }

  auto Chi2 = 0.;

  auto NonNullBin = 0;

  for(auto iBinX=1; iBinX<nBinsX+1; iBinX++){
	for(auto iBinY=1; iBinY<nBinsY+1; iBinY++) {

	  auto nObs  = hExp->GetBinContent(hExp->GetBin(iBinX, iBinY))/normExp;
	  auto nPred = hPDF->GetBinContent(hPDF->GetBin(iBinX, iBinY))/normPDF;

	  // Chi2 += EvalLL(nObs, nPred);
	  // Chi2 += EvalExtendedLL(nObs, nPred);
	  Chi2 += EvalNLL(nObs, nPred);
	}
  }

  return Chi2 /*/ static_cast<double>(nBinsX + nBinsY)*/;

}

TVector3 ScaleVector(const TVector3& v, const double& scale){
  return TVector3(v.x()*scale, v.y()*scale, v.z()*scale);
}

std::vector<double> TranslateStdVector(const std::vector<double>& vSeed,
									   const std::vector<double>& vOrig){

  if(vSeed.size() != vOrig.size()){
	std::cerr << "std::vector<double> TranslateStdVector() has different dimensions" << std::endl;
	return vSeed;
  } else {
	auto nDim = vSeed.size();
	std::vector<double> vTranslated(nDim, 0.);
	double seed, orig;
	unsigned iDim = 0;
	BOOST_FOREACH(boost::tie(seed, orig), boost::combine(vSeed, vOrig)){
			vTranslated[iDim++] = seed+orig;
		  }

	return vTranslated;

  }

}

double GetXYRho(const TVector3& Pos){
  return sqrt(Pos.x()*Pos.x() + Pos.y()*Pos.y());
}

template <typename T>
double GetXYRho(const std::vector<T>& Pos){
  if(Pos.size()>2){
	return sqrt(Pos[0]*Pos[0] + Pos[1]*Pos[1]);
  } else {
	return -1;
  }
}

struct CylVec {

  double rho, theta, z;

  CylVec() {
	rho = 0;
	theta = 0;
	z = 0;
  }

  CylVec(double rho, double theta, double z)
	  : rho(rho), theta(theta), z(z) {}

  explicit CylVec(const TVector3 &v) {

	rho = sqrt(v.x() * v.x() + v.y() * v.y());
	theta = atan2(v.y(), v.x());
	z = v.z();

  }

  explicit CylVec(const std::vector<double> &v) {

	rho = sqrt(v[0] * v[0] + v[1] * v[1]);
	theta = atan2(v[1], v[0]);
	z = v[2];

  }

  void Print() const {
	std::cout << "CylVec() with rho=" << rho << " theta=" << theta << " z=" << z << std::endl;
  }

  TVector3 GetVec3() const {
	return TVector3(rho*cos(theta), rho*sin(theta), z);
  }

  double GetMag2() const {
	return rho*rho+z*z;
  }

  double GetMag() const {
	return sqrt(GetMag2());
  }

};

static bool operator<(const CylVec& v1, const CylVec& v2){
  return v1.theta < v2.theta;
}

static bool operator>(const CylVec& v1, const CylVec& v2){
  return v1.theta > v2.theta;
}

static CylVec operator+(const CylVec& v1, const CylVec& v2){
  return {v1.rho+v2.rho, v1.theta+v2.theta, v1.z+v2.z};
}

static CylVec operator-(const CylVec& v1, const CylVec& v2){
  return {v1.rho-v2.rho, v1.theta-v2.theta, v1.z-v2.z};
}

static double GetDCylVec(const CylVec& v1, const CylVec& v2){
  return sqrt(std::pow(v2.rho-v1.rho,2) + std::pow(v2.z-v1.z,2));
}

static CylVec GetScaledCylVec(const CylVec& cv, const double& scale){
  return {cv.rho/scale, cv.theta, cv.z/scale};
}

static CylVec GetRepCylVec(const CylVec& cv, const double& scale){
  auto scv = GetScaledCylVec(cv, scale);
  return {scv.rho, scv.theta, scv.z/2 + 0.5};
}

typedef std::vector<std::vector<double> > Matrix_t;

struct Matrix {

  Matrix_t M;
  std::size_t nrows, ncols;

  Matrix(std::size_t n, std::size_t m)
	  : nrows(n), ncols(m) {
	M = Matrix_t (n, std::vector<double>(m, 0));
  }

  std::vector<double>& operator[](std::size_t idx) {return M[idx];}

  void Print(){
	for(auto i=0; i<nrows; i++){
	  for(auto j=0; j<ncols; j++){
		std::cout << M[i][j] << " ";
	  }
	  std::cout << std::endl;
	}
  }

};
typedef struct Matrix Matrix;


struct DiagMatrix {

  Matrix_t dM;
  std::size_t dim;

  explicit DiagMatrix(std::size_t n)
	  : dim(n) {
	dM = Matrix_t (n, std::vector<double>(n, 0));
  }

  double& operator[](std::size_t idx) {return dM[idx][idx];}

  void Print(){
	for(auto i=0; i<dim; i++){
	  for(auto j=0; j<dim; j++){
		std::cout << dM[i][j] << " ";
	  }
	  std::cout << std::endl;
	}
  }

};
typedef struct DiagMatrix DiagMatrix;

static double GetSignificance(const double& S, const double& B){
  return S+B != 0 ? S / std::sqrt(S + B) : 0;
}
static double ErrRate(const double& S, const double& T){
  return std::sqrt(S) / T + S / std::pow(T,2);
}

static double GetSigStatUn(const double& S, const double& B, const double& dS, const double& dB){
  if (S+B > 0) {
	double dSigdS, dSigdB;
	dSigdS = (2*B + S) / (2*std::pow(S+B, 3/2));
	dSigdB = - S / (2*std::pow(S+B, 3/2));
	return std::sqrt(std::pow(dSigdS*dS, 2) + std::pow(dSigdB*dB, 2));
  } else {
	return 0.;
  };
}

static double GetSigStatUn(const double& S, const double& B, const double& T){
  if (S+B > 0) {
	double dSigdS, dS, dSigdB, dB;
	dSigdS = (2*B + S) / (2*std::pow(S+B, 3/2));
	dSigdB = - S / (2*std::pow(S+B, 3/2));
	dS = ErrRate(S, T);
	dB = ErrRate(B, T);
	return std::sqrt(std::pow(dSigdS*dS, 2) + std::pow(dSigdB*dB, 2));
  } else {
	return 0.;
  };
}

static double PoissonScaled(const double& x, const double& lambda,
							const double& fScale){

  return std::pow(lambda / fScale, x / fScale) * std::exp(-lambda/fScale) / std::tgamma(x / fScale);

}

Double_t FitPoisson(const Double_t *x, const Double_t *par){

  Double_t xx = x[0]; Double_t lambda = par[0];

  return PoissonScaled(xx, lambda, 1000.);

}

Double_t FitGaus(const Double_t *x, const Double_t *par){

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

Double_t FitExp(const Double_t *x, const Double_t *par){

  Double_t xx = x[0];
  Double_t A = par[0]; Double_t t0 = par[1]; Double_t tau = par[2];

  Double_t arg = 0;
  if (tau!=0) arg = (xx - t0)/tau;
  Double_t fitval = par[0]*TMath::Exp(-arg);
  return fitval;

}

typedef struct zAxis {
  int nBins;
  double min, max;

  zAxis() = default;

  zAxis(int n, double mmin, double mmax)
	  : nBins(n), min(mmin), max(mmax){}

  explicit zAxis(TAxis* ax)
	  : nBins(ax->GetNbins()), min(ax->GetXmin()), max(ax->GetXmax()){

  }

  void Set(TAxis* ax){
	nBins = ax->GetNbins();
	min   = ax->GetXmin();
	max   = ax->GetXmax();
  }

} zAxis ;

std::vector<double> GetIntSpace(const unsigned int& nSteps, const double& min=-1, const double& max=1.){
  const double step = (max-min) / (double)(nSteps);
  std::vector<double> v(nSteps);
  for(auto i=0; i<nSteps;i++)
    v[i]=min+i*step;
  v.emplace_back(max);
  return v;
}

#endif //_MATHUTILS_HH_
