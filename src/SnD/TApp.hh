//
// Created by Stephane Zsoldos on 7/3/22.
//

#ifndef SND_SRC_SND_TAPP_HH_
#define SND_SRC_SND_TAPP_HH_

#include <Argz/Args.hh>

typedef struct TAppArgs : public Args {
  TAppArgs() {
	v = {
		new bArg("-v", "--verbose"),
		new sArg("-i", "--input"),
		new sArg("-o", "--output"),
		new vfArg("-t", "--tres", {250., -5., 20.})
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
  };
} TAppArgs;

#include <TH2D.h>
#include <TH1D.h>

#include "TAnalysis.hh"

class Analysis : public TAnalysis {
 private:
  std::vector< std::vector<TH2D*> > vvHPDFs;
  TH1D* hNHits;
  TH1D* hN400;
 public:
  explicit Analysis(const TAppArgs &args);
  void Do(void* Data) override;
};

#include <wRATter/Wrapper.hh>

class RATReader : public TReader {
 private:
  wRAT w_rat;
  ProgressBar progress_bar_;
  bool verbose_;
 protected:
  bool GetNextEvent() override;
  bool GetNextTrigger() override;
  void* GetData() override;
  ProgressBar *GetProgressBar() override { return &progress_bar_; }
  bool GetVerbosity() override { return verbose_; }
 public:
  explicit RATReader(const TAppArgs &args);
  ~RATReader() = default;
};


#endif //SND_SRC_SND_TAPP_HH_
