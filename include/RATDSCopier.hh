//
// Created by zsoldos on 2/24/21.
//

#ifndef SEEDNDESTROY_INCLUDE_RATDSCOPIER_HH_
#define SEEDNDESTROY_INCLUDE_RATDSCOPIER_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>
#include <vector>


// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>


// ####################################### //
// #### #### ####    RAT    #### #### #### //
// ####################################### //
#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>


typedef struct{
  TFile *f;
  TTree *T;
  TTree *RunT;
} wRATStruct;

class wRATCP {

 protected:

  wRATStruct Old;
  wRATStruct New;

  RAT::DS::Root* BufDS;

 public:
  wRATCP(const std::string& oldfilename, const std::string& newfilename, bool isVerbose = true) {

	if(isVerbose)
	  std::cout << "Copying Tree" << std::endl;

	Old.f = new TFile(oldfilename.c_str());
	Old.RunT = (TTree*)Old.f->Get("runT");
	Old.T = (TTree*)Old.f->Get("T");

	New.f = new TFile(newfilename.c_str(), "RECREATE");
	New.RunT = Old.RunT->CloneTree(0);
	New.T = Old.T->CloneTree(0);

	auto RUN = new RAT::DS::Run();
	Old.RunT->SetBranchAddress("run", &RUN);
	Old.RunT->GetEntry(0);
	New.RunT->Fill();

	BufDS = new RAT::DS::Root();
	Old.T->SetBranchAddress("ds", &BufDS);

  }

  virtual ~wRATCP() {

	delete BufDS;

  }

  void PointDSToEvt(const unsigned& iEvt) const{
	Old.T->GetEntry(iEvt);
  }

  void ClearDS(const int& iEV){
	BufDS->GetEV(iEV)->GetPathFit()->SetPosition(TVector3(0.,0.,0.));
	BufDS->GetEV(iEV)->GetPathFit()->SetDirection(TVector3(0.,0.,0.));
	BufDS->GetEV(iEV)->GetPathFit()->SetTime(0);
	BufDS->GetEV(iEV)->GetPathFit()->SetTime0(0);
	BufDS->GetEV(iEV)->GetPathFit()->SetGoodness(-1);
  }

  void FillDS(const int& iEV, const std::vector<double>& x, const double& minf){
	BufDS->GetEV(iEV)->GetPathFit()->SetPosition(TVector3(x[0], x[1], x[2]));
	if(x.size()==4){
	  BufDS->GetEV(iEV)->GetPathFit()->SetTime(x[3]/100);
	  BufDS->GetEV(iEV)->GetPathFit()->SetTime0(x[3]/100);
	} else {
	  BufDS->GetEV(iEV)->GetPathFit()->SetDirection(TVector3(x[3], x[4], x[5]));
	  BufDS->GetEV(iEV)->GetPathFit()->SetTime(x[6]/100);
	  BufDS->GetEV(iEV)->GetPathFit()->SetTime0(x[6]/100);
	}
	BufDS->GetEV(iEV)->GetPathFit()->SetGoodness(minf);
  }

  void FillDS(const int& iEV, const std::vector<double>& v){
	BufDS->GetEV(iEV)->GetPathFit()->SetPosition(TVector3(v[0], v[1], v[2]));
	BufDS->GetEV(iEV)->GetPathFit()->SetTime(v[3]);
	BufDS->GetEV(iEV)->GetPathFit()->SetTime0(v[3]);
	BufDS->GetEV(iEV)->GetPathFit()->SetGoodness(v[4]);
  }

  void FillDS(const int& iEV, const TVector3& v, const double& T, const double& minf){
	BufDS->GetEV(iEV)->GetPathFit()->SetPosition(v);
	BufDS->GetEV(iEV)->GetPathFit()->SetTime(T);
	BufDS->GetEV(iEV)->GetPathFit()->SetTime0(T);
	BufDS->GetEV(iEV)->GetPathFit()->SetGoodness(minf);
  }

  void FillTree() const{
	New.T->Fill();
  }

  void Write() const{
	New.f->cd();
	New.T->Write(nullptr,TObject::kOverwrite);
	New.RunT->Write();
	Old.f->Close();
	New.f->Close();
  }

};


#endif //SEEDNDESTROY_INCLUDE_RATDSCOPIER_HH_
