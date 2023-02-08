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

double fPosTU(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitStruct*>(data);
  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  double TGuess = x[3];
  // Calculate NLL
  return GetUNLL(d->vHits, d->hPDF,
				 PosGuess, TGuess);
}

double fPosTPerPMT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitMapStruct*>(data);
  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  double TGuess = x[3];
  // Calculate NLL
  double NLL = 0.f;
  for (const auto& hit: d->vHits){
	if(!d->mPDF.at(hit.ID))
	  continue;
	double TRes = hit.GetTRes(PosGuess, -TGuess);
	double P_TRes = d->mPDF.at(hit.ID)->Interpolate(TRes);
	NLL += P_TRes <= 0.f ? d->vHits.size() : -TMath::Log(P_TRes/d->mPDF.at(hit.ID)->Integral());
  }
  return NLL;
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

void SetBounds(nlopt::opt &opt, Bnd *c){
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
}

void SetPosBounds(nlopt::opt &opt, Bnd *c){
  // Minimizer bounds
  std::vector<double> lb = {
	  -c->GetEdge().x(), -c->GetEdge().y(), -c->GetEdge().z(), -1.f
  };
  std::vector<double> ub = {
	  c->GetEdge().x(), c->GetEdge().y(), c->GetEdge().z(), 1.f
  };
  // Set boundaries
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
}

void SetPars(nlopt::opt &opt, Bnd *c){
  // Set stopping criteria
  opt.set_xtol_rel(1.e-18);
  opt.set_ftol_rel(1.e-18);
  // Set limits
  opt.set_maxtime(1./*sec*/);
}

std::vector<RecT> DoRecon(nlopt::opt &opt, const std::vector<PosT> &vSeeds){
  //
  nlopt::result result;
  // Create return object
  std::vector< RecT > vResults;
  // Loop over seeds and recon
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
  // Sort seeds by min NLL value
  std::sort(vResults.begin(), vResults.end(),
			[](const RecT &a, const RecT &b) {
			  return a.NLL < b.NLL;
			}
  );
  // Return results
  return vResults;
}

nlopt::algorithm GetAlgo(const int &a){
  switch (a) {
	case 0:
	  return nlopt::LN_NELDERMEAD;
	  break;
	case 1:
	  return nlopt::LN_BOBYQA;
	  break;
	case 2:
	  return nlopt::LN_COBYLA;
	  break;
	case 3:
	  return nlopt::LN_NEWUOA;
	  break;
	case 4:
	  return nlopt::LN_PRAXIS;
	  break;
	case 5:
	  return nlopt::LN_SBPLX;
	  break;
	default:
	  return nlopt::LN_NELDERMEAD;
	  break;
  }
}

RecT Recon(void* data,
		   Bnd *c,
		   std::vector<PosT> &vSeeds,
		   nlopt::algorithm alg,
		   double(*fRec)(const std::vector<double> &x, std::vector<double> &grad, void *data),
		   const std::vector<void (*)(nlopt::opt &opt, Bnd *c)>& vSetPars){
  //
  // Create minimizer obj
  nlopt::opt opt(alg, 4);
  opt.set_min_objective(fRec, data);
  //
  for(auto& fSet: vSetPars)
	fSet(opt, c);
  //
  return DoRecon(opt, vSeeds).front();
}
