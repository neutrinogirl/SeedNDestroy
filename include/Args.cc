//
// Created by Stephane Zsoldos on 11/14/21.
//

#include <Args.hh>

void Args::ProcessArgs(const TApplication &theApp) {

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
	}

	if(std::none_of(v.begin(), v.end(),
				   [&arg, &theApp, &i](BaseArg *ArgPtr){
					 return (*ArgPtr)(arg, theApp, i);
				   }))
	  std::cerr << "Unkown parameter: " << arg << std::endl;

  }

}
