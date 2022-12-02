//
// Created by Stephane Zsoldos on 9/5/21.
//

#ifndef SND_INCLUDE_SND_GRIDWALKER_HH_
#define SND_INCLUDE_SND_GRIDWALKER_HH_

#include <iostream>
#include <vector>
#include <map>
#include <random>

#include <boost/foreach.hpp>
#include <boost/range/combine.hpp>

class Grid {
 private:
  std::vector<std::string> vAxNames;
  std::vector< std::vector<double> > vAx;
  std::vector< std::vector<double>::iterator > vIt;
  std::size_t nAx;
  bool isWalking;
  void Increment(const std::size_t &idx){
	if(++vIt[idx] == vAx[idx].end()){
	  if(idx==0){
		// Reset axis
		vIt[idx] = vAx[idx].begin();
		// Stop Walking
		isWalking=false;
	  } else {
		// Reset axis
		vIt[idx] = vAx[idx].begin();
		// Increment axis n-1
		Increment(idx-1);
	  }
	}
  }
 public:
  explicit Grid(std::vector<std::vector<double>> v_ax)
	  : vAx(std::move(v_ax)), isWalking(true), nAx(0) {
	for(auto& v:vAx){
	  vIt.emplace_back(v.begin());
	  nAx++;
	}
	for(auto i=0; i<nAx; i++)
	  vAxNames.emplace_back(Form("%d", i));
  }
  explicit Grid(const std::map< std::string, std::vector<double> > & v_m)
	  : isWalking(true), nAx(v_m.size()){
	AddAx(v_m);
  }
  void AddAx(const std::map< std::string, std::vector<double> > & v_m){
	for(const auto& m: v_m)
	  AddAx(m.first, m.second);
  }
  void AddAx(const std::string &name, const std::vector<double> &v){
	nAx++;
	vAxNames.emplace_back(name);
	vAx.emplace_back(v);
	vIt.emplace_back(vAx.back().begin());
  }
  void Print() const {
	std::for_each(vAx.begin(), vAx.end(), [](const std::vector<double> &v){
	  std::for_each(v.begin(), v.end(), [](const double &val){std::cout << val << " ";});
	  std::cout << std::endl;
	});
	std::for_each(vIt.begin(), vIt.end(), [](const std::vector<double>::iterator &it){
	  std::cout << *it << std::endl;
	});
  }
  bool Walk(){
	Increment(vAx.size()-1);
	return isWalking;
  }
  std::vector<double> Throw(std::mt19937 &gen){
	std::vector<double> vPt;
	for(const auto& v: vAx){
	  vPt.emplace_back(std::uniform_real_distribution<>(v.front(), v.back())(gen));
	}
	return vPt;
  }
  std::vector<double> GetGridPt() const {
	std::vector<double> val;
	val.reserve(vIt.size());
	for(const auto &v: vIt)
	  val.emplace_back(*v);
	return val;
  }
  void GetGridPt(std::vector<double> &val) const {
	if(val.size()!=nAx){
	  val.clear();
	  val.resize(nAx);
	}
	for(auto i=0; i<nAx;i++)
	  val[i]=*vIt[i];
  }
  std::map<std::string, double> GetMapPt() const {
	std::map<std::string, double> m;
	GetMapPt(m);
	return m;
  }
  void GetMapPt(std::map<std::string, double> &m) const {
	for(auto i=0; i<nAx; i++){
	  m[vAxNames[i]] = *vIt[i];
	}
  }
  size_t GetNAx() const {
	return nAx;
  }
};

#endif // SND_INCLUDE_SND_GRIDWALKER_HH_
