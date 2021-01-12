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

typedef struct Bnds {
  double Pos;
  double T;
} Bnds;

std::vector<double> ReconPosTime(DataStruct1D& DS, const Bnds& bnds,
								 const TVector3& PosSeed, const double& TSeed = 0.){

  const unsigned nDimf = 4;
  std::vector<double> x = {
	  PosSeed.x(), PosSeed.y(), PosSeed.z(),
	  TSeed*1.e2 /* Get same dimensionality as space */
  };
  double minf;

  // Create minimizer obj
  nlopt::opt opt_local(nlopt::LN_SBPLX, nDimf);
  opt_local.set_min_objective(fPosT, &DS);
  // Create result obj
  nlopt::result result_local;

  // ######################################## //
  // Create fitter boundaries
  std::vector<double> lb(nDimf, -bnds.Pos); lb[3] = (TSeed - bnds.T) * 1.e2;
  std::vector<double> ub(nDimf, bnds.Pos);  ub[3] = (TSeed + bnds.T) * 1.e2;

  // Set boundaries
  opt_local.set_lower_bounds(lb);
  opt_local.set_upper_bounds(ub);

  // Set stopping criteria
  opt_local.set_xtol_rel(1.e-6);

  // Set limits
  opt_local.set_maxtime(1./*sec*/);

  // Set step size
  opt_local.get_initial_step_(
	  {
		  1.e2, 1.e2, 1.e2,
		  1.e2
	  }
  );

  try{

	result_local = opt_local.optimize(x, minf);

  } catch (std::exception &e) {

	std::cout << "nlopt failed: " << e.what() << std::endl;

  }

  return x;

}

std::vector<double> ReconDir(DataStructDir& DS,
							 const TVector3& DirSeed){

  // Create minimizer obj
  const unsigned nDimfDir = 2;
  std::vector<double> xDir = {
	  DirSeed.Theta(),
	  DirSeed.Phi()
  };
  double minfDir;

  nlopt::opt opt_dir(nlopt::LN_SBPLX, nDimfDir);
  opt_dir.set_min_objective(fDir, &DS);

  // Create result obj
  nlopt::result result_dir;

  // ######################################## //
  // Create fitter boundaries
  std::vector<double> lbDir = {0., -PI};
  std::vector<double> ubDir = {PI, PI};

  // Set boundaries
  opt_dir.set_lower_bounds(lbDir);
  opt_dir.set_upper_bounds(ubDir);

  // Set stopping criteria
  opt_dir.set_xtol_rel(1.e-6);

  // Set limits
  opt_dir.set_maxtime(1./*sec*/);

  // Set step size
  opt_dir.get_initial_step_( { 0.1, 0.1 } );

  try {

	result_dir = opt_dir.optimize(xDir, minfDir);

  } catch (std::exception &e) {

	std::cout << "nlopt failed: " << e.what() << std::endl;

  }

  return xDir;

}

#endif //_RECON_HH_
