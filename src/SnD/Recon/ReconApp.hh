//
// Created by Stephane Zsoldos on 7/3/22.
//

#ifndef SND_SRC_SND_TAPP_HH_
#define SND_SRC_SND_TAPP_HH_

#include "Argz/Args.hh"

typedef struct ReconAppArgs : public Args {
  ReconAppArgs() {
	v = {
		new bArg("-v", "--verbose"),
		new sArg("-i", "--input"),
		new sArg("-o", "--output")
	};
  }
  ReconAppArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -i (--input) INPUT.root "
			  << " -o (--output) OUT. root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h (--help)\tShow this help message\n"
			  << "\t-v (--verbose)\tSet verbosity level true\n"

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
} ReconAppArgs;


#endif //SND_SRC_SND_TAPP_HH_
