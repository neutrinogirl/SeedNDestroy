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


std::vector<double> ReconPosTime(DataStruct1D& DS, const bnds& b, DetParams& DP,
				 const TVector3& PosSeed, const double& TSeed,
				 const double& NLLSeed = std::numeric_limits<double>::max()){

  const unsigned nDimf = 4;
  std::vector<double> x = {
    PosSeed.x() * PosScale, PosSeed.y() * PosScale, PosSeed.z() * PosScale,
    TSeed * TScale /* Get same dimensionality as space */
  };
  double minf;

  // Create minimizer obj
  nlopt::opt opt_local(nlopt::LN_COBYLA, nDimf);
  opt_local.set_min_objective(fPosT, &DS);
  // Create result obj
  nlopt::result result;

  // ######################################## //
  // Create fitter boundaries
  std::vector<double> lb = b.GetVLB();
  std::vector<double> ub = b.GetVUB();
  for(auto iDim=0; iDim<3; iDim++){
    lb[iDim] *= PosScale;
    ub[iDim] *= PosScale;
  }
  lb[3] *= TScale;
  ub[3] *= TScale;

  // Set boundaries
  opt_local.set_lower_bounds(lb);
  opt_local.set_upper_bounds(ub);
  // Set T constraints
  opt_local.add_inequality_constraint(fPosTC, &DP, 1.e-12);
  // Set stopping criteria
  opt_local.set_xtol_rel(1.e-18);
  opt_local.set_ftol_rel(1.e-18);
  // Set limits
  opt_local.set_maxtime(2./*sec*/);
  // Set step size
  opt_local.get_initial_step_({
      10., 10., 10.,
      10.
    });

  //
  // #### AUGMENTED LAGRANGIAN METHOD
  //

  nlopt::opt opt(nlopt::AUGLAG_EQ, nDimf);
  opt.set_local_optimizer(opt_local);
  opt.set_min_objective(fPosT, &DS);

  // Set boundaries
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  // Set T constraints
  opt.add_inequality_constraint(fPosTC, &DP, 1.e-12);
  // Set stopping criteria
  opt.set_xtol_rel(1.e-12);
  opt.set_ftol_rel(1.e-12);
  // Set limits
  opt.set_maxtime(1./*sec*/);
  // Set step size
  opt.get_initial_step_({
      10., 10., 10.,
      10.
    });
  
  //
  // #### 
  //

  try{

    result = opt_local.optimize(x, minf);

  } catch (std::exception &e) {

    std::cout << "nlopt failed: " << e.what() << std::endl;

  }

  x[0] /= PosScale; x[1] /= PosScale; x[2] /= PosScale;
  x[3] /= TScale;

  x.emplace_back(minf);
  x.emplace_back(result);

  return x;

}


#endif //_RECON_HH_
