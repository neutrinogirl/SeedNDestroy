//
// Created by zsoldos on 2/12/21.
//

#ifndef SEEDNDESTROY_INCLUDE_OUTPUT_HH_
#define SEEDNDESTROY_INCLUDE_OUTPUT_HH_

#include <iostream>
#include <vector>

#include <TVector3.h>
#include <TTree.h>

#include <Hit.hh>

typedef struct Vec{
  double x=0., y=0., z=0.;
  Vec() = default;
  explicit Vec(const TVector3& v) : x(v.x()), y(v.y()), z(v.z()) { };
  explicit Vec(const std::vector<double>& v) : x(v[0]), y(v[1]), z(v[2]) { };
} Vec;

#define MAXNHITS 10000

typedef struct Event{
  Vec MCPos;
  double MCT=0.;

  Vec RecPos;
  double RecT=0.;

  Vec MCDir;

  double ETrue=0.;

  int NHits=0;
  double X[MAXNHITS], Y[MAXNHITS], Z[MAXNHITS], T[MAXNHITS], Q[MAXNHITS];

  double Chi2=0.;
  int NLOPT=0;

  int iEvt=0;
  int iTrig=0;

  Event() = default;
  Event(const TVector3& mcpos, const TVector3& mcdir, const double& mct,
		const TVector3& recpos, const double& recT, const double& etrue, const int& nhits,
		const double& chi2, const int& nlopt)
	  : MCPos(mcpos), MCDir(mcdir), MCT(mct), RecPos(recpos), RecT(recT),
		ETrue(etrue), NHits(nhits),
		Chi2(chi2), NLOPT(nlopt){ };
} Event;

void SetTTree(TTree& Tree, Event& Event){

  Tree.Branch("ievt", &Event.iEvt, "ievt/I");
  Tree.Branch("itrig", &Event.iTrig, "itrig/I");

  Tree.Branch("mcx", &Event.MCPos.x, "mcx/D");
  Tree.Branch("mcy", &Event.MCPos.y, "mcy/D");
  Tree.Branch("mcz", &Event.MCPos.z, "mcz/D");
  Tree.Branch("mcT", &Event.MCT, "mcT/D");

  Tree.Branch("recx", &Event.RecPos.x, "recx/D");
  Tree.Branch("recy", &Event.RecPos.y, "recy/D");
  Tree.Branch("recz", &Event.RecPos.z, "recz/D");
  Tree.Branch("recT", &Event.RecT, "recT/D");

  Tree.Branch("mcdx", &Event.MCDir.x, "mcdx/D");
  Tree.Branch("mcdy", &Event.MCDir.y, "mcdy/D");
  Tree.Branch("mcdz", &Event.MCDir.z, "mcdz/D");

  Tree.Branch("etrue", &Event.ETrue, "etrue/D");
  Tree.Branch("nhits", &Event.NHits, "nhits/I");

  Tree.Branch("xhits", &Event.X, "xhits[nhits]/D");
  Tree.Branch("yhits", &Event.Y, "yhits[nhits]/D");
  Tree.Branch("zhits", &Event.Z, "zhits[nhits]/D");
  Tree.Branch("thits", &Event.T, "thits[nhits]/D");
  Tree.Branch("qhits", &Event.Q, "qhits[nhits]/D");

  Tree.Branch("chi2", &Event.Chi2, "chi2/D");
  Tree.Branch("nlopt", &Event.NLOPT, "nlopt/I");

}

void LoadVHits(Event& evt, const std::vector<Hit>& vHits){

  evt.NHits = vHits.size();

  for(auto i=0; i<vHits.size(); i++){

    evt.X[i] = vHits[i].PMTPos.x();
	evt.Y[i] = vHits[i].PMTPos.y();
	evt.Z[i] = vHits[i].PMTPos.z();

	evt.T[i] = vHits[i].T;
	evt.Q[i] = vHits[i].Q;

  }

}

template <typename T>
void Save2ROOT(const T* obj, const std::string& filename = "test.root"){

  TFile f(filename.c_str(), "UPDATE");
  if(!f.IsZombie()){
	obj->Write();
	f.Close();
  } else {
	std::cerr << "Can't save obj to " << filename << " because it's a ZOMBIE" << std::endl;
  }

}

#endif //SEEDNDESTROY_INCLUDE_OUTPUT_HH_
