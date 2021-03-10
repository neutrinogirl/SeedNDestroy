//
// Created by zsoldos on 1/5/21.
//

#ifndef _READFILE_HH_
#define _READFILE_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>
#include <utility>

// ####################################### //
// #### #### ####   BOOST   #### #### #### //
// ####################################### //
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TApplication.h>

// ####################################### //
// #### #### ####   USER    #### #### #### //
// ####################################### //
#include <Utils.hh>
#include <MathUtils.hh>

std::vector<std::string> GetFilesInDir(const boost::filesystem::path& dir, const std::string& ext = ".root"){
  std::vector<std::string> vPaths;

  if(boost::filesystem::exists(dir) && boost::filesystem::is_directory(dir)){
    for(const auto& entry: boost::filesystem::recursive_directory_iterator(dir)){
	  if (boost::filesystem::is_regular_file(entry) && entry.path().extension() == ext){
		vPaths.push_back(dir.string() + "/" + entry.path().filename().string());
	  }
    }

  }

  return vPaths;

}

typedef struct Args{
  bool isVerbose = false;
  std::vector<std::string> filename;
  std::string outname = "PDF.root";
  unsigned int nEvts = 0;
  TVector3 bnds = TVector3(10.e3, 10.e3, 10.e3);
  bool isBox = false;
} Args;

static void ShowUsage(const std::string& name){

  std::cerr << "Usage: " << name << " <option(s)> -i (--input) IN.root -o (--output) OUT.root" << std::endl
			<< "Options:\n"

			<< "\t-h\tShow this help message\n"
			<< "\t-v\tSet verbose mode (display progress bar)\n"
			<< "\t-n <nEvts>\tSet nEvts to process\n"
			<< "\t-b <XX YY ZZ>\tSet boundaries for box geom (in mm) \n"
			<< "\t-c <R H>\tSet boundaries for cylinder geom (in mm) \n"

			<< "\t--dir\tRead all .root files in directory \n"


			<< std::endl;

}

static void ProcessArgs(TApplication &theApp,
						Args &args) {

  // Reading user input parameters
  if (theApp.Argc() < 2) {
	ShowUsage(theApp.Argv(0));
	exit(0);
  }

  for (int i = 1; i < theApp.Argc(); i++) {
	std::string arg = theApp.Argv(i);
	if ((arg == "-h") || (arg == "--help")) {
	  ShowUsage(theApp.Argv(0));
	  exit(0);
	} else if (boost::iequals(arg, "-v")) {
	  args.isVerbose=true;
	} else if (boost::iequals(arg, "-n")) {
	  args.nEvts=std::stoi(theApp.Argv(++i));
	} else if (boost::iequals(arg, "-b")) {
      args.bnds.SetX(std::stod(theApp.Argv(++i)));
      args.bnds.SetY(std::stod(theApp.Argv(++i)));
      args.bnds.SetZ(std::stod(theApp.Argv(++i)));
      std::cout << "Setting box geom boundaries" << std::endl;
      args.bnds.Print();
      args.isBox = true;
	} else if (boost::iequals(arg, "-c")) {
	  const double R = std::stod(theApp.Argv(++i));
	  const double Z = std::stod(theApp.Argv(++i));
	  args.bnds.SetX(1);
	  args.bnds.SetY(0);
	  args.bnds.SetZ(0);
	  std::cout << "Setting cylinder geom boundaries" << std::endl;
	  args.bnds.SetPerp(R);
	  args.bnds.SetZ(Z);
	  args.bnds.Print();
	  args.isBox = false;

	} else if (boost::iequals(arg,"-i") || boost::iequals(arg,"--input")) {
	  args.filename.emplace_back(theApp.Argv(++i));
	} else if (boost::iequals(arg,"--dir")) {
	  auto v = GetFilesInDir(theApp.Argv(++i));
	  std::merge(args.filename.begin(), args.filename.end(), v.begin(), v.end(), std::back_inserter(args.filename));
	} else if (boost::iequals(arg,"-o") || boost::iequals(arg,"--output")) {
	  args.outname=theApp.Argv(++i);

	} else {
	  std::cout << "Unkown parameter" << std::endl;
	  continue;
	}

  }


  auto RemoveEmptyFile = [](std::vector<std::string>& v){
	for(auto itFile=v.begin(); itFile!=v.end(); itFile++){
	  if(!IsFileExist((*itFile).c_str())){
		v.erase(itFile);
	  }
	}
  };

  RemoveEmptyFile(args.filename);

  if(args.filename.empty()){
	std::cerr << "ERROR: No input file provided!" << std::endl;
	exit(EXIT_FAILURE);
  }

}

#endif //_READFILE_HH_
