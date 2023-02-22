//
// Created by Stephane Zsoldos on 2/22/23.
//

#include "SnD/SpaceGrid.hh"

SpaceGrid::SpaceGrid(const std::vector<std::vector<double>> &vv) {
  // Init vPts and vNLL
  std::cout << "SpaceGrid::SpaceGrid(const std::vector<std::vector<float>> &vv)" << std::endl;
  for(const auto& x:vv[0])
	for(const auto& y:vv[1])
	  for(const auto& z:vv[2])
		vPts.emplace_back(x,y,z);
  nPts = vv[0].size()*vv[1].size()*vv[2].size();
  vNLL.resize(nPts, 0.);
  std::cout << "SpaceGrid::SpaceGrid(const std::vector<std::vector<float>> &vv) - Done" << std::endl;
}

void SpaceGrid::Reset() {
  std::fill(vNLL.begin(), vNLL.end(), 0.f);
}

void SpaceGrid::Clear() {
  vNLL.clear();
}

static float GetNLLFlat(const std::vector<double>& vPDF,
						const std::vector<Hit>& vHits, const TVector3& Pt, const float &T,
						ZHist& ZHTres){
  //
  ZHTres.Reset();
  //
  for(const auto& hit:vHits){
	ZHTres.AddData(hit.GetTRes(Pt, T));
  }
  //
  float L=0.f;
  for(int i=0; i<ZHTres.GetNumBins(); i++){
	float nObs = ZHTres.GetData()[i];
	float nPred = vPDF[i];
	if(nObs>0 && nPred>0)
	  L+=nObs*TMath::Log(nObs/nPred) + nPred-nObs;
	else
	  L+=nPred;

  }
  //
  return L;
}

//
static std::vector<double> getBinContents(TH1D* hist) {
  int numBins = hist->GetNbinsX();
  const double* binArray = hist->GetArray();
  std::vector<double> binContents(binArray, binArray + numBins + 2);
  return binContents;
}
// Get N bins, min and max from TH1D
static std::tuple<int, double, double> getNbinsMinMax(TH1D* hist) {
  int numBins = hist->GetNbinsX();
  double min = hist->GetXaxis()->GetXmin();
  double max = hist->GetXaxis()->GetXmax();
  return std::make_tuple(numBins, min, max);
}
// Remove first and last element from std::vector<double>
static std::vector<double> removeFirstLast(const std::vector<double>& v) {
  std::vector<double> v2(v.size() - 2);
  std::copy(v.begin() + 1, v.end() - 1, v2.begin());
  return v2;
}

void Loop(const std::vector<double>& vPDF,
		  const std::vector<Hit>& vHits, const std::vector<TVector3>& vPts, const float &T,
		  std::vector<float>& vNLL,
		  ZHist ZHTres,
		  int startIndex, int endIndex){
  //
  for (int i = startIndex; i < endIndex; i++) {
	vNLL[i] = GetNLLFlat(vPDF, vHits, vPts[i], T, ZHTres);
  }
  //
}

void SpaceGrid::Walk(TH1D* hPDF,
					 const std::vector<Hit> &vHits, const float& T) {

  // Init vPDF
  std::vector<double> vPDF = removeFirstLast(getBinContents(hPDF));
  auto t = getNbinsMinMax(hPDF);
  // Scale vPDF to size vHits
  std::vector<double> vNormPDF = vPDF;
  std::transform(vNormPDF.begin(), vNormPDF.end(), vNormPDF.begin(),
				 [n = vHits.size()](double x) { return x * n; });
  //
  int startIndex = 0;
  int endIndex = nPts;
  //
  ZHist ZHTres(std::get<0>(t), std::get<1>(t), std::get<2>(t));
  //
  Loop(vNormPDF, vHits, vPts, T, vNLL, ZHTres, startIndex, endIndex);

}

void SpaceGrid::ParallelWalk(TH1D* hPDF,
							 const std::vector<Hit> &vHits, const float& T) {

  // Init vPDF
  std::vector<double> vPDF = removeFirstLast(getBinContents(hPDF));
  auto t = getNbinsMinMax(hPDF);
  // Scale vPDF to size vHits
  std::vector<double> vNormPDF = vPDF;
  std::transform(vNormPDF.begin(), vNormPDF.end(), vNormPDF.begin(),
				 [n = vHits.size()](double x) { return x * n; });
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
	ZHist ZHTres(std::get<0>(t), std::get<1>(t), std::get<2>(t));
	futures[i] = std::async(
		std::launch::async, Loop,
		std::cref(vNormPDF),
		std::cref(vHits), std::cref(vPts), std::cref(T),
		std::ref(vNLL),
		ZHTres,
		startIndex, endIndex
	);
  }
  for (int i = 0; i < numThreads; i++) {
	futures[i].wait();
  }
  //
}
