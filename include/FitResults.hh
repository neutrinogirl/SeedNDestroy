//
// Created by zsoldos on 1/13/21.
//

#ifndef SEEDNDESTROY_INCLUDE_FITRESULTS_HH_
#define SEEDNDESTROY_INCLUDE_FITRESULTS_HH_

#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

// General structure to hold params
typedef struct Params{
  double val=0, err=0;
  Params() = default;
  Params(double val, double err) : val(val), err(err) {}
  Params(const Params& P) = default;
  Params(std::initializer_list<Params>) {}
  void Print() const {
	std::cout << val << "+-" << err << std::endl;
  }
} Params;

// #################################### //
// #### DEDICATED TO GAUSSIAN FITS #### //
// #################################### //

namespace GausFit {

typedef struct Res {
  double Chi2=0;
  Params Norm, Mu, Sig;
  Res() = default;
  Res(double chi_2, const Params &norm, const Params &mu, const Params &sig)
	  : Chi2(chi_2), Norm(norm), Mu(mu), Sig(sig) {}
  Res(const Res& R) = default;
  Res(std::initializer_list<Res>) {}
} Res;

static Double_t FitGaus(const Double_t *x, const Double_t *par){

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

template <typename T> // TH1*
Res GetFitRes(T* h, const std::vector<double>& vBnds = {}){

  std::vector<double> bnds =
	  vBnds.empty() ? (vBnds.size() == 2 ? vBnds : std::vector<double>{ h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() } )
					: std::vector<double>{ h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() };
  TF1 fG("fG", FitGaus, bnds[0], bnds[1], 3);
  fG.SetParameter(0, h->GetMaximum());
  fG.SetParameter(1, h->GetMean());
  fG.SetParameter(2, h->GetRMS());
  TFitResultPtr fPRes = h->Fit("fG", "LQEMRNS");

  return {
	  fPRes.Get()->Chi2() / fPRes.Get()->Ndf(),
	  {fG.GetParameter(0), fG.GetParError(0)},
	  {fG.GetParameter(1), fG.GetParError(1)},
	  {fG.GetParameter(2), fG.GetParError(2)}
  };

}

}

#endif //SEEDNDESTROY_INCLUDE_FITRESULTS_HH_
