//
// Created by Stephane Zsoldos on 7/3/22.
//

#ifndef SND_SRC_APPS_RECON_HH_
#define SND_SRC_APPS_RECON_HH_

#include "Argz/Args.hh"

typedef struct ReconAppArgs : public Args {
  ReconAppArgs() {
	v = {
		new bArg("-v",   "--verbose"),
		new sArg("-i",   "--input"),
		new sArg("-p",   "--pdf"),
		new sArg("-o",   "--output"),
		new fArg("-r",   "--radius"),
		new fArg("-hh",  "--hheight"),
		new bArg("-u",   "--unbinned"),
		new sArg("-pn",  "--pdf-name",     "hCTVSTResPDF_TTOF_QW0"),
		new sArg("-ppn", "--pdf-pmt-name", "hCTVSTResPDF_TTOF_QW0_PMT"),
		new iArg("-n",   "--n-evts",       -1),
		new iArg("-a",   "--algo",          0),
		new iArg("-ms",  "--max-seed",     -1),
		new bArg("-m",   "--map"),
		new bArg("-vv",  "--vverbose"),
		new sArg("-mn",  "--map-name",     "MAP.root"),
		new bArg("-b",   "--binned"),
		new bArg("-pp",  "--pperpmt"),
	};
  }
  explicit ReconAppArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -r (--radius) R"
			  << " -hh (--half-height) HH"
			  << " -p (--pdf) PDF.root"
			  << " -i (--input) INPUT.root"
			  << " -o (--output) OUT.root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h   (--help)        \tShow this help message\n"
			  << "\t-v   (--verbose)     \tSet verbosity level true\n"
			  << "\t-vv  (--vverbose)    \tSet MEGA verbosity\n"
			  << "\t-b   (--binned)      \tSet binned TRes fit\n"
			  << "\t-u   (--unbinned)    \tSet unbinned TRes fit\n"
			  << "\t-pp  (--pperpmt)     \tSet PDF p PMT TRes fit (unbinned)\n"
			  << "\t-pn  (--pdf-name)    \tSet PDF hist name\n"
			  << "\t-ppn (--pdf-pmt-name)\tSet PDFpPMT hist name\n"
			  << "\t-n   (--n-evts)      \tSet n evts to process (default all)\n"
			  << "\t-ms  (--max-seed)    \tSelect max seeds to try (default all)\n"
			  << "\t-m   (--map)         \tPlot and save NLL map in MAP.root\n"
			  << "\t-mn  (--map-name)    \tSet output name map (default MAP.root)\n"
			  << "\t-a   (--algo)        \tSelect algorithm (default 0):\n"
			  << "\t                     \t [0]: nlopt::LN_NELDERMEAD\n"
			  << "\t                     \t [1]: nlopt::LN_BOBYQA\n"
			  << "\t                     \t [2]: nlopt::LN_COBYLA\n"
			  << "\t                     \t [3]: nlopt::LN_NEWUOA\n"
			  << "\t                     \t [4]: nlopt::LN_PRAXIS\n"
			  << "\t                     \t [5]: nlopt::LN_SBPLX\n"

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
  const char* GetPDFPMTName() const {
	return reinterpret_cast<sArg*>(v[8])->val.c_str();
  }
  int GetNEvts() const {
	return reinterpret_cast<iArg*>(v[9])->val;
  }
  int GetAlgo() const {
	return reinterpret_cast<iArg*>(v[10])->val;
  }
  int GetMaxSeed() const {
	return reinterpret_cast<iArg*>(v[11])->val;
  }
  bool GetMap() const {
	return reinterpret_cast<bArg*>(v[12])->val;
  }
  bool GetVVerbose() const {
	return reinterpret_cast<bArg*>(v[13])->val;
  }
  const char *GetMapName() const {
	return reinterpret_cast<sArg*>(v[14])->val.c_str();
  }
  bool GetBinned() const {
	return reinterpret_cast<bArg*>(v[15])->val;
  }
  bool GetPerPMT() const {
	return reinterpret_cast<bArg*>(v[16])->val;
  }
} ReconAppArgs;


#endif //SND_SRC_APPS_RECON_HH_
