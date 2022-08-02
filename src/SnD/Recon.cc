//
// Created by Stephane Zsoldos on 8/2/22.
//

#include <SnD/Recon.hh>

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitStruct*>(data);
  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  double TGuess = x[3];
  // Calculate NLL
  return GetNLL(d->vHits, d->hPDF,
				PosGuess, TGuess);
}

double fPosTC(const std::vector<double> &x, std::vector<double> &grad, void *data) {
  //
  auto d = static_cast<Cylinder *>(data);
  //
  TVector3 PosGuess(x[0], x[1] ,x[2]);
  double TGuess = x[3];
  double TWall = d->GetTWall(PosGuess);
  //
  return TGuess > d->GetTEdge() ? -1.f : std::min(TWall, TGuess);
}

RecT Recon(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *c, std::vector<PosT> &vSeeds){
  //
  FitStruct FS = {vHits, hPDF};
  // Create minimizer obj
  nlopt::opt opt(nlopt::LN_COBYLA, 4);
  opt.set_min_objective(fPosT, &FS);
  nlopt::result result;
  // Minimizer bounds
  std::vector<double> lb = {
	  -c->GetEdge().x(), -c->GetEdge().y(), -c->GetEdge().z(), 0.f
  };
  std::vector<double> ub = {
	  c->GetEdge().x(), c->GetEdge().y(), c->GetEdge().z(), c->GetTEdge()
  };
  // Set boundaries
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  // Set T constraints
  opt.add_inequality_constraint(fPosTC, c, 1.e-12);
  // Set stopping criteria
  opt.set_xtol_rel(1.e-18);
  opt.set_ftol_rel(1.e-18);
  // Set limits
  opt.set_maxtime(1./*sec*/);

  std::vector< RecT > vResults;

  std::transform(
	  vSeeds.begin(), vSeeds.end(),
	  std::back_inserter(vResults),
	  [&](const PosT &s) {
		double minf=std::numeric_limits<double>::max();
		std::vector<double> x = s.GetStdVec();
		try{
		  result = opt.optimize(x, minf);
		} catch (std::exception &e) {
		  std::cout << "nlopt failed: " << e.what() << std::endl;
		}
		return RecT(x[0], x[1], x[2], x[3], minf);
	  }
  );

  std::sort(vResults.begin(), vResults.end(),
		  [](const RecT &a, const RecT &b) {
			return a.NLL < b.NLL;
		  }
  );

  return vResults.front();
}

