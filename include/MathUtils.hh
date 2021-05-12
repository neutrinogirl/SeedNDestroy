//
// Created by zsoldos on 6/17/20.
//

#ifndef _MATHUTILS_HH_
#define _MATHUTILS_HH_

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <utility>
#include <bitset>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>

#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>

#define PI 3.14159265359
#define SQRT2 1.41421356237
#define RINDEX_WATER 1.33
#define SOL_VACUUM 299.792
#define SOL SOL_VACUUM/RINDEX_WATER

const double GetPI(){
  return PI;
}

const double GetSQRT2(){
  return SQRT2;
}

const double GetSOL(){
  return SOL;
}

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

Double_t FitGaus(const Double_t *x, const Double_t *par){

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
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

  Bnd() = default;
  Bnd(T min, T max) : min(min), max(max) {}
  Bnd(const Bnd<T>& rhs) {
    min = rhs.min;
    max = rhs.max;
  }
  // Bnd(Bnd<T>&& rhs) noexcept {
  // min = std::move(rhs.min);
  // max = std::move(rhs.max);
  // }

  bool IsIn(const T& val) const {
    return val > min && val < max;
  }

};

template <typename T, typename U>
struct Bnds {

  std::vector<Bnd<T>> vPos;
  Bnd<U> vT;

  Bnds() = default;
  Bnds(const std::vector<Bnd<T>> &v_pos, const Bnd<U> &v_t) : vPos(v_pos), vT(v_t) {}
  Bnds(const Bnds& rhs){
    vPos = rhs.vPos;
    vT = rhs.vT;
  }
  // Bnds(Bnds&& rhs){
  // vPos = std::move(rhs.vPos);
  // vT = std::move(rhs.vT);
  // }

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
    v.push_back(-vT.min);
    return v;
  }

};

typedef struct Bnds<double, double> bnds;

struct CylBnds : public bnds{

  CylBnds(){
    vPos = { {0., 0.}, {0., 0.} };
    vT = {0., 0.};
  };
  explicit CylBnds(const bnds& b){
    vPos = b.vPos;
    vT = b.vT;
  }
  CylBnds(const double& radius, const double& hheight, const double& SoL = GetSOL()){
    vPos = { {0, radius}, {0, hheight} };
    vT = {0, std::min(radius, hheight) / SoL};
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
    return {GetRadius(), GetRadius(), GetHHeight(), - vT.min};
  }

};

struct BoxBnds : public bnds{

  BoxBnds(){
    vPos = { {0, 0}, {0, 0}, {0, 0} };
    vT = {0, 0};
  };
  explicit BoxBnds(const bnds& b){
    vPos = b.vPos;
    vT = b.vT;
  }
  explicit BoxBnds(const std::vector<double>& vvPos, const double& SoL = GetSOL()) : BoxBnds(){
    if(vPos.size() != vvPos.size())
      std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
    double min = *std::min_element(vvPos.begin(), vvPos.end());
    for(auto i=0; i<vvPos.size();i++){
      vPos[i] = {0, vvPos[i]};
    }
    vT = {0, min/SoL};
  }
  explicit BoxBnds(const std::vector<std::vector<double>>& vvPos, const double& SoL = GetSOL()) : BoxBnds(){
    if(vPos.size() != vvPos.size())
      std::cerr << "BoxBnds::BoxBnds() NOT ALL DIMS ARE SPECIFIED" << std::endl;
    auto min =
      *std::min_element(vvPos.begin(), vvPos.end(),
			[](const std::vector<double>&v1, const std::vector<double>&v2){
			  return *std::min_element(v1.begin(), v1.end()) < *std::min_element(v2.begin(), v2.end());
			}
			);
    for(auto i=0; i<vvPos.size();i++){
      vPos[i] = {vvPos[i][0], vvPos[i][1]};
    }
    vT = {0, *std::min_element(min.begin(), min.end())/SoL};
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


static std::vector< std::vector<double> >
Get4DSmplGuess(const std::vector<double>& xGuess,
	       const std::vector<double>& vScale){

  const unsigned int nDim = 4;
  const std::size_t nPts = std::pow(2, nDim);
  std::vector< std::vector<double> > vSmplGuess; vSmplGuess.reserve(nPts+1);
  vSmplGuess.emplace_back(xGuess);

  auto GetBnd = [](const unsigned int& b, const double& x){
    return b ? x : -x;
  };

  for(auto i=0; i<nPts;i++){
    std::bitset<nDim> word(i);
    std::vector<double> v(nDim);
    for(auto iDim=0; iDim<nDim; iDim++){
      v[iDim] = xGuess[iDim] + GetBnd(word[iDim], vScale[iDim]);
    }
    vSmplGuess.emplace_back(v);
  }

  return vSmplGuess;

}
static std::vector< std::vector<double> >
Get4DSmplGuess(const std::vector<double>& xGuess,
	       const double& Scale){

  const unsigned int nDim = 4;
  const std::size_t nPts = std::pow(2, nDim);
  std::vector< std::vector<double> > vSmplGuess; vSmplGuess.reserve(nPts+1);
  vSmplGuess.emplace_back(xGuess);

  auto GetBnd = [](const unsigned int& b, const double& x){
    return b ? x : -x;
  };

  for(auto i=0; i<nPts;i++){
    std::bitset<nDim> word(i);
    std::vector<double> v(nDim);
    for(auto iDim=0; iDim<nDim; iDim++){
      v[iDim] = xGuess[iDim] + GetBnd(word[iDim], Scale);
    }
    vSmplGuess.emplace_back(v);
  }

  return vSmplGuess;

}

#endif //_MATHUTILS_HH_
