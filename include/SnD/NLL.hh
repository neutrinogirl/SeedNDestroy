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
double CalculateLL(T const *hPDF, T const *hExp,
				   bool isNormalized){

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

  return Chi2;
}

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const TVector3& Pos, const double& T,
			  bool isNormalized=false);

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const std::vector<double> &x,
			  bool isNormalized=false);

double GetUNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			   const TVector3& Pos, const double& T);

double GetUNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			   const std::vector<double> &x);

#endif //SEEDNDESTROY_INCLUDE_NLL_HH_
