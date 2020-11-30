//
// Created by zsoldos on 11/25/20.
//

#ifndef _SEEDNDESTROY_HH_
#define _SEEDNDESTROY_HH_

#include <iostream>

#include <TVector3.h>
#include <TFile.h>

#include "PathFit.hh"
#include "wRATter/include/Hit.hh"

typedef struct SnDRes {

  TVector3 Pos;
  double T;
  double Chi2;

  void Print() const {
	Pos.Print();
	std::cout << "T:" << T << "ns "
			  << "Chi2:" << Chi2 << std::endl;
  }

} SnDRes;

TVector3 GetAvgDirFromPosSeed(const TVector3& PosSeed, const std::vector<Hit>& vHits,
							  double(*fW)(const Hit&, int) = fweight, int wPower = 0){

  TVector3 v;
  for(auto& h: vHits){
	v+=h.PMTPos.Unit()*fW(h,wPower);
  }
  return v;

}

template <typename T>
static T* GetRootHisto(const char* filename, const char* histname){
  auto f = TFile::Open(filename);
  // Check if key exist
  if(!f->GetListOfKeys()->Contains(histname))
	return nullptr;
  auto hist = dynamic_cast<T *>(f->Get(histname)->Clone());
  hist->SetDirectory(nullptr);
  f->Close();
  delete f;
  return hist;
}


#endif //_SEEDNDESTROY_HH_
