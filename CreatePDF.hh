//
// Created by zsoldos on 1/5/21.
//

#ifndef _CREATEPDF_HH_
#define _CREATEPDF_HH_

// ####################################### //
// #### #### ####   USER    #### #### #### //
// ####################################### //
#include <Args.hh>
typedef struct CreatePDFArgs : public Args {
  CreatePDFArgs() {
	v = {
		new bArg("-v", "--verbose"),
		new sArg("-i", "--input"),
		new sArg("-o", "--output"),
		new iArg("-n", "--nevts", 0),
		new fArg("-r", "--radius"),
		new fArg("-hh", "--half-height"),
		new vfArg("-t", "--tres", {250., -5., 20.})
	};
  }
  CreatePDFArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -r (--radius) R "
			  << " -hh (--half-height) HH "
			  << " -i (--input) INPUT.root "
			  << " -o (--output) OUT. root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h (--help)\tShow this help message\n"
			  << "\t-v (--verbose)\tSet verbosity level true\n"
			  << "\t-n (--nevts) N\tSet n evts to process (default process all evts)\n"
			  << "\t-t (--tres) nBins min max\tSet bins for TRes hist (default 250, -5., 20.)\n"

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
  int GetNEvts() const {
	return reinterpret_cast<iArg*>(v[3])->val;
  }
  float GetRadius() const {
	return reinterpret_cast<fArg*>(v[4])->val;
  }
  float GetHHeight() const {
	return reinterpret_cast<fArg*>(v[5])->val;
  }
  std::vector<float> GetTResBins() const {
	return reinterpret_cast<vfArg*>(v[6])->val;
  };
} CreatePDFArgs;

#endif //_CREATEPDF_HH_
