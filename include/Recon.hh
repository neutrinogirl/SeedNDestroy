//
// Created by zsoldos on 10/12/20.
//

#ifndef _RECON_HH_
#define _RECON_HH_

#include <iostream>

#include <TVector3.h>

#include <nlopt.hpp>

#include "PathFit.hh"
#include "MathUtils.hh"


std::vector<double> ReconPosTime(DataStruct1D& DS, const Bnds& bnds,
								 const TVector3& PosSeed, const double& TSeed = 0.){

  const unsigned nDimf = 4;
  std::vector<double> x = {
	  PosSeed.x() * PosScale, PosSeed.y() * PosScale, PosSeed.z() * PosScale,
	  TSeed * TScale /* Get same dimensionality as space */
  };
  double minf;

  // Create minimizer obj
  nlopt::opt opt_local(nlopt::LN_SBPLX, nDimf);
  opt_local.set_min_objective(fPosT, &DS);
  // Create result obj
  nlopt::result result_local;

  // ######################################## //
  // Create fitter boundaries
  std::vector<double> lb(nDimf); lb[3] = - bnds.T[1] * TScale;
  std::vector<double> ub(nDimf); ub[3] = - bnds.T[0] * TScale;
  for(auto iDim=0; iDim<3; iDim++){
    lb[iDim] = ( - bnds.Pos[iDim]) * PosScale;
	ub[iDim] = ( + bnds.Pos[iDim]) * PosScale;
  }

  // std::cout << "[" << lb[3] << "," << ub[3] << "] " << x[3] << std::endl;

  // Set boundaries
  opt_local.set_lower_bounds(lb);
  opt_local.set_upper_bounds(ub);

  // Set stopping criteria
  opt_local.set_xtol_rel(1.e-3);

  // Set limits
  opt_local.set_maxtime(1./*sec*/);

  // Set step size
  opt_local.get_initial_step_(
	  {
		  10., 10., 10.,
		  10.
	  }
  );

  try{

	result_local = opt_local.optimize(x, minf);

  } catch (std::exception &e) {

	std::cout << "nlopt failed: " << e.what() << std::endl;

  }

  x[0] /= PosScale; x[1] /= PosScale; x[2] /= PosScale;
  x[3] /= TScale;

  x.emplace_back(minf);
  x.emplace_back(result_local);

  return x;

}

#endif //_RECON_HH_
