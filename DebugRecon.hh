//
// Created by zsoldos on 1/5/21.
//

#ifndef _DEBUGRECON_HH_
#define _DEBUGRECON_HH_

// ####################################### //
// #### #### ####   USER    #### #### #### //
// ####################################### //
#include <Args.hh>
typedef struct ReconArgs : public Args {
  ReconArgs() {
	v = {
		new bArg("-v", "--verbose"),
		new bArg("-d", "--debug"),
		new bArg("-dd", "--ddebug"),

		new sArg("-i", "--input"),
		new sArg("-p", "--pdf"),
		new sArg("-o", "--output"),

		new iArg("-n", "--nevts", 0),

		new iArg("-w", "--weight", 0),

		new fArg("-r", "--radius"),
		new fArg("-hh", "--half-height"),

		new bArg("-u", "--unbinned"),

		new fArg("-c", "--cvg", 0.),
	};
  }
  ReconArgs(const std::vector<BaseArg *> &v) : Args(v) {}
  void ShowUsage(const std::string &name) override {
	std::cout << "Usage: " << name
			  << " <option(s)>"
			  << " -r (--radius) R "
			  << " -hh (--half-height) HH "
			  << " -i (--input) INPUT.root "
			  << " -p (--pdf) PDF.root "
			  << " -o (--output) OUT. root\n\n"

			  << "Options: [default]\n\n"

			  << "\t-h (--help)    \tShow this help message\n"
			  << "\t-v (--verbose) \tSet verbosity level true\n"
			  << "\t-d (--debug)   \tSet debug plots for each evt\n"
			  << "\t-dd (--ddebug) \tSet debug prints for each evt\n"
			  << "\t-n (--nevts)  N\tSet n evts to process [process all evts]\n"
			  << "\t-w (--weight) N\tSet TRes kernel reweighting [0]\n"
			  << "\t-u (--unbinned)\tSet unbinned TRes fit\n"
			  << "\t-c (--cvg) N   \tSet artifical cvg\n"
			  << std::endl;
  }
  bool GetVerbose() const {
	return reinterpret_cast<bArg*>(v[0])->val;
  }
  bool GetDebug() const {
	return reinterpret_cast<bArg*>(v[1])->val;
  }
  bool GetDDebug() const {
	return reinterpret_cast<bArg*>(v[2])->val;
  }
  const char *GetInput() const {
	return reinterpret_cast<sArg*>(v[3])->val.c_str();
  }
  const char *GetPDF() const {
	return reinterpret_cast<sArg*>(v[4])->val.c_str();
  }
  const char *GetOutput() const {
	return reinterpret_cast<sArg*>(v[5])->val.c_str();
  }
  int GetNEvts() const {
	return reinterpret_cast<iArg*>(v[6])->val;
  }
  int GetWeight() const {
	return reinterpret_cast<iArg*>(v[7])->val;
  }
  float GetRadius() const {
	return reinterpret_cast<fArg*>(v[8])->val;
  }
  float GetHHeight() const {
	return reinterpret_cast<fArg*>(v[9])->val;
  }
  bool GetUnbinned() const {
	return reinterpret_cast<bArg*>(v[10])->val;
  }
  float GetCvg() const {
	return reinterpret_cast<fArg*>(v[11])->val;
  }
} ReconArgs;

template <typename T>
T* GetRootHisto(const char* filename, const char* histname){
  auto f = TFile::Open(filename);
  // Check if key exist
  if(!f->GetListOfKeys()->Contains(histname))
    return nullptr;
  auto hist = dynamic_cast<T *>(f->Get(histname)->Clone());
  hist->SetDirectory(nullptr);
  f->Close();
  delete f;
  return hist;
}

std::vector<double> GetTBndsLocal(const double& TSeed,
				  const std::vector<double>& vTBnds, const double& range = 3.){
  double TMin = TSeed - range > vTBnds[0] ? TSeed - range : vTBnds[0] ;
  double TMax = TSeed + range < vTBnds[1] ? TSeed + range : vTBnds[1] ;
  return {TMin, TMax};
}

void SlimVHits(std::vector<Hit>& vHOrig, const double ScaleFactor = 0.9){

  const auto nHitsScaled = static_cast<unsigned>( std::round(vHOrig.size() * (ScaleFactor / 0.9) ) );

  auto rng = std::default_random_engine {};
  std::shuffle(vHOrig.begin(), vHOrig.end(), rng);

  vHOrig.resize(nHitsScaled);

}

//the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


void PrintDD(const std::string& name,
	     const TVector3& Pos, const double& T, const double& dWall,
	     const double& NLL = 0.){

  std::cout << "#### #### #### " << name << " #### #### ####" << std::endl;
  Pos.Print();
  std::cout << RED << T << "ns" << RESET << std::endl;
  std::cout << BLUE << dWall << "mm" << RESET << std::endl;
  if(NLL>0.)
    std::cout << GREEN << "NLL=" << NLL << RESET << std::endl;
  std::cout << "#### #### #### -------- #### #### ####" << std::endl;

}


#endif //_DEBUGRECON_HH_
