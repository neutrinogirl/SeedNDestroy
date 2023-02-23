//
// Created by Stephane Zsoldos on 2/23/23.
//

#include "Algo/WOpt.hh"

double fLS(const std::vector<double> &x, std::vector<double> &grad, void *data){
  // Cast to Vector of Hits
  auto d = static_cast<std::vector<Hit>*>(data);

  // Create ZVector to handle unit conversions
  // Since 1ns ~ 1dm, we can use the same unit for both time and space
  Vector3<double> v(x[0], x[1], x[2], SpaceUnit::dm);
  return GetSum2Residuals(v.GetTVector3(), x[3], *d);
}

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitStruct*>(data);
  // Create object to calculate TRes histogram
  Vector3<double> v(x[0], x[1], x[2], SpaceUnit::dm);
  double TGuess = x[3];
  // Calculate NLL
  return GetNLL(*d->hPDF, v.GetTVector3(), TGuess, d->vHits);
}

double fPosTU(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitStruct*>(data);
  // Create object to calculate TRes histogram
  Vector3<double> v(x[0], x[1], x[2], SpaceUnit::dm);
  double TGuess = x[3];
  // Calculate NLL
  return GetUNLL(*d->hPDF, v.GetTVector3(), TGuess, d->vHits);
}

double fPosTPerPMT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  //
  auto d = static_cast<FitMapStruct*>(data);
  // Create object to calculate TRes histogram
  Vector3<double> v(x[0], x[1], x[2], SpaceUnit::dm);
  double TGuess = x[3];
  return GetUNLL(d->mPDF, v.GetTVector3(), TGuess, d->vHits);
}

double fLSC(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<Cylinder *>(data);

  Vector3<double> v(x[0], x[1], x[2], SpaceUnit::dm);
  double TWall = d->GetTWall(v.GetTVector3());

  // Inequality constraint
  return std::abs(TWall - x[3]) - 1;
}

void SetBounds(nlopt::opt &opt, Bnd *c){
  // Create ZVector to handle unit conversions
  // Since 1ns ~ 1dm, we can use the same unit for both time and space
  Vector3<double> v3(c->GetEdge().x(), c->GetEdge().y(), c->GetEdge().z(), SpaceUnit::mm);
  Vector3<double> v3dm = v3.ConvertTo(SpaceUnit::dm);
  // Minimizer bounds
  std::vector<double> lb = {
	  -v3dm.GetX(), -v3dm.GetY(), -v3dm.GetZ(), -1.f
  };
  std::vector<double> ub = {
	  v3dm.GetX(), v3dm.GetY(), v3dm.GetZ(), c->GetTEdge()+1.f
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

void SetInequalityConstraint(nlopt::opt &opt, Bnd *c){
  opt.add_inequality_constraint(fLSC, c, 1.e-8);
}

std::vector< RecT > DoRecon(nlopt::opt &opt, const std::vector< PosT > &vSeeds){
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
	  return nlopt::LN_COBYLA;
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
