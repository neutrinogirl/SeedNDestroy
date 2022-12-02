//
// Created by Stephane Zsoldos on 7/5/22.
//

#ifndef SND_INCLUDE_SND_ZAXIS_HH_
#define SND_INCLUDE_SND_ZAXIS_HH_

#include <TAxis.h>

typedef struct zAxis {
  int nBins;
  double min, max;

  zAxis() = default;

  zAxis(int n, double mmin, double mmax)
	  : nBins(n), min(mmin), max(mmax){}

  explicit zAxis(TAxis* ax)
	  : nBins(ax->GetNbins()), min(ax->GetXmin()), max(ax->GetXmax()){

  }

  void Set(TAxis* ax){
	nBins = ax->GetNbins();
	min   = ax->GetXmin();
	max   = ax->GetXmax();
  }

  std::vector<double> GetStdVec() const{
	std::vector<double> v(nBins+1);
	const double width = (max-min) / static_cast<double>(nBins);
	std::iota(v.begin(), v.end(), 0);
	std::transform(v.begin(), v.end(), v.begin(), [&](const double& val){
	  return min + (val+0.5)*width;
	});
	return v;
  }

} zAxis ;

#endif //SND_INCLUDE_SND_ZAXIS_HH_
