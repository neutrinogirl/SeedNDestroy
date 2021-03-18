//
// Created by zsoldos on 6/17/20.
//

#ifndef _MATHUTILS_HH_
#define _MATHUTILS_HH_

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>

#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>
#include <utility>

#define PI 3.14159265359
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
	normPDF = hPDF->Integral();
	normExp = hExp->Integral();
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

  std::vector<double> GetStdVec() const{
    std::vector<double> v(nBins+1);
    const double width = (max-min) / static_cast<double>(nBins);
    std::iota(v.begin(), v.end(), 0);
	std::transform(v.begin(), v.end(), v.begin(), [&](const double& val){
	  return min + (val+0.5)*width;
	});
	return v;
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

template<typename T>
void ScaleHist(T *hist, const double& norm = 0){
  if(norm > 0){
	hist->Scale(1./norm);
  } else {
	hist->Scale(1./hist->Integral());
  }
}

template <typename T>
struct Bnd{

  T min, max;

  bool IsIn(const T& val) const {
    return val > min && val < max;
  }

};

template <typename T, typename U>
struct Bnds {

  std::vector<Bnd<T>> vPos;
  Bnd<U> vT;

  virtual bool IsInPos(const std::vector<T>& v) const {

    if(v.size() != vPos.size()){
	  std::cerr << "Bnds::InInPos() DIM BOUNDS != DIM V" << std::endl;
	  return false;
    }

    for(auto i=0; i<v.size(); i++){
      if(!vPos[i].IsIn(v[i]))
        return false;
    }
	return true;

  }

  virtual bool IsInT(const U& val) const {
    return vT.IsIn(val);
  }

  virtual bool IsInPos(const TVector3& v) const { return false; }
  virtual T GetDWall(const TVector3& v) const { return -std::numeric_limits<T>::max(); }
  virtual TVector3 GetTVector3() const { return TVector3(); }
  virtual void Print() const {
    for(const auto& bnd: vPos){
      std::cout << "[" << bnd.min << ", " << bnd.max << "]mm ";
    }
    std::cout << "[" << vT.min << ", " << vT.max << "]ns" << std::endl;
  }
  virtual T GetMaxDWall() const {
	T maxdwall = std::numeric_limits<T>::max();
	for (const auto &bnd: vPos) {
	  maxdwall = std::min(maxdwall, std::min(std::abs(bnd.min), std::abs(bnd.max)));
	}
	return maxdwall;
  }
  virtual std::vector<T> GetVDetBnds() const {
    std::vector<T> v;
	for(const auto& bnd: vPos) {
	  v.push_back(std::max(std::abs(bnd.min), std::abs(bnd.max)));
	}
	return v;
  }

  virtual std::vector<T> GetVLB() const {
	std::vector<T> v;
	for(const auto& bnd: vPos) {
	  v.push_back(bnd.min);
	}
	v.push_back(-vT.max);
	return v;
  }

  virtual std::vector<T> GetVUB() const {
	std::vector<T> v;
	for(const auto& bnd: vPos) {
	  v.push_back(bnd.max);
	}
	v.push_back(vT.min);
	return v;
  }

};

typedef struct Bnds<double, double> bnds;

struct CylBnds : public bnds{

  CylBnds(){
	vPos = { {0, 0}, {0, 0} };
	vT = {0, 0};
  };
  CylBnds(const double& radius, const double& hheight){
	vPos = { {0, radius}, {0, hheight} };
	vT = {0, sqrt(2*std::pow(radius, 2) + std::pow(hheight, 2)) / SOL};
  }
  CylBnds(const double& radius, const double& hheight,
		  const std::vector<double>& TBnds){
	vPos = { {0, radius}, {0, hheight} };
	vT = {TBnds[0], TBnds[1]};
  }
  CylBnds(const double& radius, const double& hheight,
		  const double& TMin, const double& TMax){
	vPos = { {0, radius}, {0, hheight} };
	vT = {TMin, TMax};
  }

  double GetRadius() const{
    return vPos[0].max;
  }

  double GetHHeight() const{
	return vPos[1].max;
  }

  bool IsInPos(const TVector3&v) const override {
    return vPos[0].IsIn(v.Perp()) && vPos[1].IsIn(std::abs(v.z()));
  }

  double GetDWall(const TVector3& v) const override {
	return std::min(GetRadius() - v.Perp(), GetHHeight() - std::abs(v.z()));
  }

  TVector3 GetTVector3() const override {
    return TVector3(GetRadius(), GetRadius(), GetHHeight());
  }

  double GetMaxDWall() const override {
    return std::min(GetRadius(), GetHHeight());
  }

  std::vector<double> GetVDetBnds() const override {
    return {GetRadius(), GetRadius(), GetHHeight()};
  }

  std::vector<double> GetVLB() const override{
    return {-GetRadius(), -GetRadius(), -GetHHeight(), - vT.max};
  }

  std::vector<double> GetVUB() const override{
	return {GetRadius(), GetRadius(), GetHHeight(), vT.min};
  }

};

struct BoxBnds : public bnds{

  BoxBnds(){
	vPos = { {0, 0}, {0, 0}, {0, 0} };
	vT = {0, 0};
  };
  BoxBnds(const std::vector<double>& vvPos) : BoxBnds(){
	if(vPos.size() != vvPos.size())
	  std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
	double max2 = 0;
	for(auto i=0; i<vvPos.size();i++){
	  vPos[i] = {0, vvPos[i]};
	  max2 += std::pow(vvPos[i], 2);
	}
	vT = {0, sqrt(max2)};
  }
  BoxBnds(const std::vector<std::vector<double>>& vvPos) : BoxBnds(){
	if(vPos.size() != vvPos.size())
	  std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
	double max2 = 0;
	for(auto i=0; i<vvPos.size();i++){
	  vPos[i] = {vvPos[i][0], vvPos[i][1]};
	  max2 += std::pow(vvPos[i][1], 2);
	}
	vT = {0, sqrt(max2)};
  }
  BoxBnds(const std::vector<std::vector<double>>& vvPos,
		  const std::vector<double>& TBnds) : BoxBnds(){
	if(vPos.size() != vvPos.size())
	  std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
	for(auto i=0; i<vvPos.size();i++){
	  vPos[i] = {vvPos[i][0], vvPos[i][1]};
	}
	vT = {TBnds[0], TBnds[1]};
  }
  BoxBnds(const std::vector<std::vector<double>>& vvPos,
		  const double& TMin, const double& TMax) : BoxBnds(){
	if(vPos.size() != vvPos.size())
	  std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
	for(auto i=0; i<vvPos.size();i++){
	  vPos[i] = {vvPos[i][0], vvPos[i][1]};
	}
	vT = {TMin, TMax};
  }

  bool IsInPos(const TVector3&v) const override {
    for(auto i=0; i<3; i++){
	  if(!vPos[i].IsIn(v[i]))
	    return false;
    }
    return true;
  }

  double GetDWall(const TVector3& v) const override {
    double dWall = std::numeric_limits<double>::max();
	for(auto i=0; i<3; i++) {
	  double absv = std::abs(v[i]);
	  double min = std::min(std::abs(vPos[0].min) - absv, std::abs(vPos[0].max) - absv);
	  dWall = min < dWall ? min : dWall;
	}
	return dWall;
  }

  TVector3 GetTVector3() const override {
	return TVector3(vPos[0].max, vPos[1].max, vPos[2].max);
  }

};


template <typename T>
T GetDWall(const TVector3& v,
		   const T& radius, const T& hheight)  {
  return std::min(radius - v.Perp(), hheight - std::abs(v.z()));
}

template <typename T>
T GetDWall(const TVector3& v,
		   const T& radius, const T& hheight,
		   std::size_t& idx)  {
  std::vector<T> vv = {radius - v.Perp(), hheight - std::abs(v.z())};
  auto min_element = std::min_element(vv.begin(), vv.end());
  idx = std::distance(vv.begin(), min_element);
  return *min_element;
}

template <typename T>
T GetDWall(const TVector3& v,
		   const std::vector<T>& vDims,
		   std::size_t& idx)  {
  std::vector<T> vv = {vDims[0] - v.Perp(), vDims[1] - std::abs(v.z())};
  auto min_element = std::min_element(vv.begin(), vv.end());
  idx = std::distance(vv.begin(), min_element);
  return *min_element;
}

#endif //_MATHUTILS_HH_
