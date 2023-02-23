//
// Created by Stephane Zsoldos on 2/22/23.
//

#include "SnD/SpaceGrid.hh"

SpaceGrid::SpaceGrid(const std::vector<std::vector<double>> &vv) {
  // Init vPts and vNLL
  for(const auto& x:vv[0])
	for(const auto& y:vv[1])
	  for(const auto& z:vv[2])
		vPts.emplace_back(x,y,z);
  nPts = vv[0].size()*vv[1].size()*vv[2].size();
  vNLL.resize(nPts, 0.);
}

void SpaceGrid::Reset() {
  std::fill(vNLL.begin(), vNLL.end(), 0.f);
}

void SpaceGrid::Clear() {
  vNLL.clear();
}

void Loop(const TH1D& hPDF,
		  const std::vector<Hit>& vHits, const std::vector<TVector3>& vPts, const float &T,
		  std::vector<float>& vNLL,
		  int startIndex, int endIndex){
  //
  for (int i = startIndex; i < endIndex; i++) {
	vNLL[i] = GetUNLL(hPDF, vPts[i], T, vHits);
  }
  //
}

void SpaceGrid::Walk(const TH1D& hPDF,
					 const std::vector<Hit> &vHits, const float& T) {
  //
  int startIndex = 0;
  int endIndex = nPts;
  //
  Loop(hPDF, vHits, vPts, T, vNLL, startIndex, endIndex);
  //
}

void SpaceGrid::ParallelWalk(const TH1D& hPDF,
							 const std::vector<Hit> &vHits, const float& T) {

  //
  const int numThreads = std::thread::hardware_concurrency();
  const int numElementsPerThread = nPts / numThreads;
  //
  std::vector<std::future<void>> futures(numThreads);
  for (int i = 0; i < numThreads; i++) {
	//
	int startIndex = i * numElementsPerThread;
	int endIndex = (i + 1) * numElementsPerThread;
	//
	futures[i] = std::async(
		std::launch::async, Loop,
		std::ref(hPDF),
		std::cref(vHits), std::cref(vPts), std::cref(T),
		std::ref(vNLL),
		startIndex, endIndex
	);
  }
  for (int i = 0; i < numThreads; i++) {
	futures[i].wait();
  }
  //
}
