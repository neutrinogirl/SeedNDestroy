//
// Created by Stephane Zsoldos on 2/21/23.
//

#ifndef SND_INCLUDE_SND_SPACEGRID_HH_
#define SND_INCLUDE_SND_SPACEGRID_HH_

// // // // // //
#include "Hit.hh"

#include <vector>
#include <future>

#include <TH1D.h>
#include <TVector3.h>

// Generate linear space between two values min/max and a constant step
// https://stackoverflow.com/questions/27028226/linear-spaced-vector-in-c
template<typename T>
std::vector<T> linspace(T min, T max, std::size_t N){
  T h = (max - min) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = min; x != xs.end(); ++x, val += h)
	*x = val;
  return xs;
}

//
class ZHist {
 private:
  std::vector<float> data;
  int numBins;
  double binSize;
  double minValue;
 public:
  ZHist() = default;
  ZHist(int bins, double min, double max) {
	Init(bins, min, max);
  }
  void Init(int bins, double min, double max) {
	numBins = bins;
	minValue = min;
	binSize = (max - min) / bins;
	data.resize(numBins, 0.f);
  }
  void Reset() {
	std::fill(data.begin(), data.end(), 0.f);
  }
  void AddData(double value) {
	if (value < minValue || value > (minValue + numBins * binSize))
	  return; // ignore values outside of range
	int bin = (int)((value - minValue) / binSize);
	data[bin]++;
  }
  void Print() const {
	for (int i = 0; i < numBins; i++) {
	  std::cout << minValue + i * binSize << " - " << minValue + (i + 1) * binSize << ": " << data[i] << std::endl;
	}
  }
  // Get data
  std::vector<float> GetData() const { return data; }
  // Get number of bins
  int GetNumBins() const { return numBins; }
  // Get bin size
  double GetBinSize() const { return binSize; }
  // Get min value
  double GetMinValue() const { return minValue; }
};

//
class SpaceGrid{
 private:
  std::size_t nPts;
 public:
  std::vector<TVector3> vPts;
  std::vector<float>    vNLL;
  //
  SpaceGrid(const std::vector< std::vector<double> >& vv);
  //
  void Reset();
  void Clear();
  //
  void Walk(TH1D* hPDF,
			const std::vector<Hit>& vHits, const float& T);
  void ParallelWalk(TH1D* hPDF,
					const std::vector<Hit>& vHits, const float& T);
  //
  std::size_t GetNPts() const { return nPts; }
};


#endif //SND_INCLUDE_SND_SPACEGRID_HH_
