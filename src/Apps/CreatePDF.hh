//
// Created by Stephane Zsoldos on 7/3/22.
//

#ifndef SND_SRC_APPS_CREATEPDF_HH_
#define SND_SRC_APPS_CREATEPDF_HH_

#include "Argz/Args.hh"

typedef struct TAppArgs : public Args {
  TAppArgs() {
	v = {
		new bArg("-v", "--verbose"),
		new sArg("-i", "--input"),
		new sArg("-o", "--output"),
		new vfArg("-t", "--tres", {250., -5., 20.}),
		new bArg("-s", "--shift"),
		new vfArg("-ps", "--pos-shift", {0, 0, 0}),
		new bArg("-a", "--apply-trigger")
	};
  }
  TAppArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -i (--input) INPUT.root "
			  << " -o (--output) OUT. root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h (--help)\tShow this help message\n"
			  << "\t-v (--verbose)\tSet verbosity level true\n"
			  << "\t-t (--tres) nBins min max\tSet bins for TRes hist (default 250, -5., 20.)\n"
			  << "\t-s (--shift)\tShift TRes hist to 0 (default false)\n"
			  << "\t-ps (--pos-shift) x y z\tShift MC TRUE positions (default 0, 0, 0)\n"
			  << "\t-a (--apply-trigger)\tApply trigger (default false)\n"

			  << std::endl;
  }
  bool GetVerbose() const {
	return reinterpret_cast<bArg*>(v[0])->val;
  }
  const char *GetInput() const {
	return reinterpret_cast<sArg*>(v[1])->val.c_str();
  }
  const char *GetOutput() const {
	return reinterpret_cast<sArg*>(v[2])->val.c_str();
  }
  std::vector<float> GetTResBins() const {
	return reinterpret_cast<vfArg*>(v[3])->val;
  }
  bool GetShift() const {
	return reinterpret_cast<bArg*>(v[4])->val;
  }
  std::vector<float> GetPosShift() const {
	return reinterpret_cast<vfArg*>(v[5])->val;
  }
  bool GetApplyTrigger() const {
	return reinterpret_cast<bArg*>(v[6])->val;
  }
} TAppArgs;


#endif //SND_SRC_APPS_CREATEPDF_HH_
