//
// Created by zsoldos on 9/15/20.
//

#ifndef _PATHFIT_HH_
#define _PATHFIT_HH_

#include <TH1D.h>
#include <TH2D.h>

#include <Hit.hh>
#include <TRandom3.h>
#include "MathUtils.hh"
#include "Recorder.hh"
#include "SnD/NLL.hh"

typedef struct DataStruct{
  std::vector<Hit> vHits;
  unsigned int wPower;
  bool isUnbinned;
  DataStruct(unsigned int w_power, bool is_unbinned)
	  : wPower(w_power), isUnbinned(is_unbinned) {}
} DataStruct;

typedef struct DataStruct1D : public DataStruct, Recorder {
  TH1D* hPDF;
  DataStruct1D(unsigned int w_power, bool is_unbinned, TH1D *h_pdf)
	  : DataStruct(w_power, is_unbinned), hPDF(h_pdf), Recorder() {}
} DataStruct1D;

const double PosScale = 1.e-1;
const double TScale   = 1.e1;

typedef struct DetParams {
  bnds *b;
  double SoL;
} DetParams;

double fPosTC(const std::vector<double> &x, std::vector<double> &grad, void *data) {
  auto d = static_cast<DetParams *>(data);

  TVector3 PosGuess(x[0] / PosScale, x[1] / PosScale, x[2] / PosScale);
  double TGuess = x[3] / TScale;
  double dWall = d->b->GetDWall(PosGuess);

  return dWall < 0.f? dWall : std::abs(TGuess) - dWall/d->SoL;
}

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<DataStruct1D*>(data);

  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0] / PosScale, x[1] / PosScale, x[2] / PosScale);
  double TGuess = x[3] / TScale;

  // Calculate NLL
  double NLL = GetNLL(d->vHits, d->hPDF, PosGuess, TGuess, fweight, d->wPower, d->isUnbinned);

  // Record
  d->iCall++;
  d->vPosGuess.emplace_back(PosGuess);
  d->vTGuess.emplace_back(TGuess);
  d->vNLL.emplace_back(NLL);

  return NLL;

}

double fPosTSmear(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<DataStruct1D*>(data);

  auto vx = Get4DSmplGuess(x, 1.);
  double NLL=0.;

  for(auto& xx:vx){
    // Create object to calculate TRes histogram
    TVector3 PosGuess(xx[0] / PosScale, xx[1] / PosScale, xx[2] / PosScale);
    double TGuess = xx[3] / TScale;
    NLL += GetNLL(d->vHits, d->hPDF, PosGuess, TGuess, fweight, d->wPower, true) / static_cast<double>(vx.size());
  }

  return NLL;

}


#endif //_PATHFIT_HH_
