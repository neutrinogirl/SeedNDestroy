//
// Created by zsoldos on 1/15/21.
//

#ifndef SEEDNDESTROY_INCLUDE_NLL_HH_
#define SEEDNDESTROY_INCLUDE_NLL_HH_

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TMath.h>

#include <random>

#include <SnD/ZAxis.hh>
#include <SnD/Hit.hh>

static std::string random_string(const unsigned int& nChars=32);

static double EvalLL(double nObs, double nPred);
static double EvalExtendedLL(double nObs, double nPred);
static double EvalNLL(double nObs, double nPred);

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
double UnbinnedLL(TH1D *hPDF, const std::vector<double> &vTRes);

double GetNLL(const std::vector<Hit>& vHits, TH2D* hPDF,
			  const TVector3& Pos, const double& T, const TVector3& Dir,
			  double(*fW)(const Hit&, const int&) = fWeight, const int& wPower = 1);
double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const TVector3& Pos, const double& T,
			  double(*fW)(const Hit&, const int&) = fWeight, const int& wPower = 1,
			  const bool &isUnbinned = false);

#endif //SEEDNDESTROY_INCLUDE_NLL_HH_
