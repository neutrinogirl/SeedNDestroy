//
// Created by zsoldos on 10/12/20.
//

#ifndef _RECON_HH_
#define _RECON_HH_

#include <SnD/PosT.hh>

#include <nlopt.hpp>

PosT Recon(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *c, std::vector<PosT> &vSeeds){

  auto fPosT = [&](const std::vector<double> &x, std::vector<double> &grad, void *data) {
	// Create object to calculate TRes histogram
	TVector3 PosGuess(x[0], x[1], x[2]);
	double TGuess = x[3];

	// Calculate NLL
	return GetNLL(vHits, hPDF,
				  PosGuess, TGuess);

  };

  nlopt::opt opt(nlopt::LN_NELDERMEAD, 4);
  opt.set_min_objective(fPosT, nullptr);

  std::map< double, PosT > mRec;

  std::transform(
	  vSeeds.begin(), vSeeds.end(),
	  mRec.begin(),
	  [&](PosT &s) {
		double minf=std::numeric_limits<double>::max();
		int result=0;
		std::vector<double> x = s.GetStdVec();
		try{
		  result = opt_local.optimize(x, minf);
		} catch (std::exception &e) {
		  std::cout << "nlopt failed: " << e.what() << std::endl;
		}
		return std::make_pair(minf, PosT(TVector3(x[0], x[1], x[2]), x[3]));
	  }
  );

  return mRec.begin()->second;

}

#endif //_RECON_HH_
