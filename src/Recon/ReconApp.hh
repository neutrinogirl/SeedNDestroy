//
// Created by Stephane Zsoldos on 7/3/22.
//

#ifndef SND_SRC_SND_TAPP_HH_
#define SND_SRC_SND_TAPP_HH_

#include "Argz/Args.hh"

typedef struct ReconAppArgs : public Args {
  ReconAppArgs() {
	v = {
		new bArg("-v",  "--verbose"),
		new sArg("-i",  "--input"),
		new sArg("-p",  "--pdf"),
		new sArg("-o",  "--output"),
		new fArg("-r",  "--radius"),
		new fArg("-hh", "--hheight"),
		new bArg("-u",  "--unbinned"),
		new sArg("-pn", "--pdf-name", "hCTVSTResPDF_TTOF_QW0"),
		new iArg("-n", "--n-evts", -1)
	};
  }
  ReconAppArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -r (--radius) R"
			  << " -hh (--half-height) HH"
			  << " -p (--pdf) PDF.root"
			  << " -i (--input) INPUT.root"
			  << " -o (--output) OUT.root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h  (--help)    \tShow this help message\n"
			  << "\t-v  (--verbose) \tSet verbosity level true\n"
			  << "\t-u  (--unbinned)\tSet unbinned TRes fit\n"
			  << "\t-pn (--pdf-name)\tSet PDF hist name\n"
			  << "\t-n  (--n-evts)  \tSet n evts to process (default all)\n"

			  << std::endl;
  }
  bool GetVerbose() const {
	return reinterpret_cast<bArg*>(v[0])->val;
  }
  const char *GetInput() const {
	return reinterpret_cast<sArg*>(v[1])->val.c_str();
  }
  const char *GetPDF() const {
	return reinterpret_cast<sArg*>(v[2])->val.c_str();
  }
  const char *GetOutput() const {
	return reinterpret_cast<sArg*>(v[3])->val.c_str();
  }
  float GetRadius() const {
	return reinterpret_cast<fArg*>(v[4])->val;
  }
  float GetHHeight() const {
	return reinterpret_cast<fArg*>(v[5])->val;
  }
  bool GetUnbinned() const {
	return reinterpret_cast<bArg*>(v[6])->val;
  }
  const char *GetPDFName() const {
	return reinterpret_cast<sArg*>(v[7])->val.c_str();
  }
  int GetNEvts() const {
	return reinterpret_cast<iArg*>(v[8])->val;
  }
} ReconAppArgs;


#endif //SND_SRC_SND_TAPP_HH_
