//
// Created by zsoldos on 2/23/21.
//

#ifndef SEEDNDESTROY_INCLUDE_TRIGGERTIMEMAP_HH_
#define SEEDNDESTROY_INCLUDE_TRIGGERTIMEMAP_HH_

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

template<typename T>
class AxisGrid {

  std::vector<T> vEdges;
  std::vector<T> vCenters;

 public:

  //
  // #### #### #### CONSTRUCTORS / DESTRUCTORS #### #### #### //
  //

  AxisGrid(std::vector<T> vBnds, T width) {

	if(vBnds.empty())
	  std::cerr << "AxisGrid:: vBnds is empty" << std::endl;

	if(vBnds.size() == 1){
	  vBnds.template emplace_back(0);
	}

	if(vBnds.size() != 2)
	  std::cerr << "AxisGrid:: vBnds is > 2" << std::endl;

	// Ensure vBnds is sorted
	std::sort(vBnds.begin(), vBnds.end());
	const std::size_t nBins = std::ceil( (vBnds[1] - vBnds[0]) / width);
	vEdges.resize(nBins+1);
	vCenters.resize(nBins);
	std::iota(vEdges.begin(), vEdges.end(), 0);
	std::iota(vCenters.begin(), vCenters.end(), 0);

	std::transform(vEdges.begin(), vEdges.end(), vEdges.begin(),
				   [&](const double& val) {return vBnds[0] + val*width;});
	std::transform(vCenters.begin(), vCenters.end(), vCenters.begin(),
				   [&](const double& val) {return vBnds[0] + (val+0.5)*width;});

  }

  //
  // #### #### #### GETTERS / SETTERS #### #### #### //
  //

  const std::vector<T> &GetVEdges() const { return vEdges; }
  const std::vector<T> &GetVCenters() const { return vCenters; }

};

class TrigTimePDF{

 private:

  std::vector<int> vRho;
  std::vector<int> vZ;
  std::map< int , std::map< int , TH1D*> > mT;

  int FindNearest(const std::vector<int>& v, const int val){

	auto itLb = std::lower_bound(v.begin(), v.end(), val);

	if(itLb == v.begin()){
	  return v[0];
	}

	if(itLb != v.end()){

	  const auto idx = std::distance(v.begin(), itLb);

	  if(idx>0){

		const auto lb = val - v[idx-1];
		const auto ub = v[idx] - val;

		return lb > ub ? v[idx] : v[idx-1];

	  }

	}

	return -1;

  }

  template <typename T>
  T* GetRootHisto(TFile* f, const char* histname){
	// Check if key exist
	if(!f->GetListOfKeys()->Contains(histname))
	  return nullptr;
	auto hist = dynamic_cast<T *>(f->Get(histname)->Clone());
	hist->SetDirectory(nullptr);
	return hist;
  }


 public:

  //
  // #### #### #### CONSTRUCTORS / DESTRUCTORS #### #### #### //
  //

  TrigTimePDF() = default;

  TrigTimePDF(std::vector<int > v_rho, std::vector<int > v_z) : vRho(std::move(v_rho)), vZ(std::move(v_z)) {
	for(auto& R:vRho){
	  for(auto& Z:vZ){
		mT.insert(std::make_pair(R, std::map<int , TH1D*>()));
		mT[R].insert(std::make_pair(Z, new TH1D(Form("hR%dZ%d", R, Z),
												"",
												100, -10., 90.)));
		mT[R][Z]->SetDirectory(nullptr);
	  }
	}
  }

  //
  // #### #### #### METHODS #### #### #### //
  //

  void Fill(const TVector3& v, const double& TrigTime){

	int R = std::round(v.Perp());
	int Z = std::round(std::abs(v.z()));

	R = FindNearest(vRho, R);
	Z = FindNearest(vZ, Z);

	if(R>-1 && Z>-1)
	  mT[R][Z]->Fill(TrigTime);

  }

  void Save(){
	TTree T("TTriggerTimeMap", "Tree with trigger map bins");
	T.Branch("vRho", &vRho);
	T.Branch("vZ", &vZ);
	T.Fill();
	T.Write();
  }

  void Load(const std::string& filename){
    TFile f(filename.c_str(), "READ");

    TTree *T;
    f.GetObject("TTriggerTimeMap", T);
    std::vector<int> *pvRho=nullptr, *pvZ=nullptr;

	TBranch *bvpRho=nullptr, *bvpZ=nullptr;
    T->SetBranchAddress("vRho", &pvRho, &bvpRho);
	T->SetBranchAddress("vZ", &pvZ, &bvpZ);

	auto tentry = T->LoadTree(0);
	bvpRho->GetEntry(tentry);
	bvpZ->GetEntry(tentry);
	delete T;

	vRho = *pvRho;
	vZ = *pvZ;

	for(auto& R: vRho){
	  for(auto& Z: vZ){
	    mT.insert(std::make_pair(R, std::map<int , TH1D*>()));
		mT[R].insert(std::make_pair(Z, GetRootHisto<TH1D>(&f, Form("hR%dZ%d", R, Z))));
	  }
	}

  }

  double GetTrigTime(const TVector3& v, std::vector<double>& vBnds){
	int R = std::round(v.Perp());
	int Z = std::round(std::abs(v.z()));

	R = FindNearest(vRho, R);
	Z = FindNearest(vZ, Z);

	if(R>-1 && Z>-1){
	  vBnds = {mT[R][Z]->GetMean() - std::abs(mT[R][Z]->GetRMS()), mT[R][Z]->GetMean() + std::abs(mT[R][Z]->GetRMS())};
	  return mT[R][Z]->GetMean();
	}

	return 0.;

  }

  //
  // #### #### #### GETTERS / SETTERS #### #### #### //
  //

  const std::vector<int> &GetVRho() const { return vRho; }
  const std::vector<int> &GetVZ() const { return vZ; }
  std::map<int, std::map<int, TH1D *>> &GetMT() { return mT; }

};


#endif //SEEDNDESTROY_INCLUDE_TRIGGERTIMEMAP_HH_
