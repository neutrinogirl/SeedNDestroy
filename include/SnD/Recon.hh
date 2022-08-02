//
// Created by zsoldos on 10/12/20.
//

#ifndef _RECON_HH_
#define _RECON_HH_

#include <SnD/PosT.hh>
#include <SnD/Geom.hh>
#include <SnD/Hit.hh>
#include <SnD/NLL.hh>

#include <nlopt.hpp>

typedef struct FitStruct{
  std::vector<Hit> vHits;
  TH1D *hPDF;
} FitStruct;

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data);
double fPosTC(const std::vector<double> &x, std::vector<double> &grad, void *data);

typedef struct FitResults {
  double NLL=0.f;
  PosT RecT;
  FitResults() = default;
  FitResults(double NLL, PosT RecT) : NLL(NLL), RecT(RecT) {}
  void SetTTree(TTree *t){
	t->Branch("NLL", &NLL, "NLL/D");
  }
} FitResults;

FitResults Recon(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *c, std::vector<PosT> &vSeeds);

#endif //_RECON_HH_
