//
// Created by zsoldos on 1/5/21.
//

#include <iostream>

#include "ReadFile.hh"
#include <Wrapper.hh>
#include <ProgressBar.hpp>

int main(int argc, char *argv[]){

  // ######################################## //
  // Create TApp
  TApplication theApp("App", &argc, argv);


  // ######################################## //
  // Parse arguments
  // Simple struct containing filename and verbosity level
  Args args;
  ProcessArgs(theApp, args);
  const bool isVerbose = args.isVerbose;
  const std::string input = args.filename;


  // ######################################## //
  // Create wrapper object
  wRAT w_rat(input);
  const unsigned long int nEvts = w_rat.GetNEvts();


  // ######################################## //
  // Loop and get vector of NHits
  ProgressBar progress_bar(nEvts, 70);
  for(auto iEvt=0; iEvt<nEvts; iEvt++){

    // Record the tick
    ++progress_bar;

    // Point to evt
    w_rat.SetEvt(iEvt);

    // Get number of trigger associated with an event
    // i.e, number of EV inside the rat DS
	auto nTriggers = w_rat.GetNTriggers();

	for(auto iTrigger=0; iTrigger<nTriggers; iTrigger++){

	  // Get vector of hits
	  std::vector<Hit> vHits = w_rat.GetVHits(iTrigger);

	  // DO STUFF
	  std::cout << vHits.size() << std::endl;
	  // ...

	}

	if(isVerbose)
	  progress_bar.display();

  }

  if(isVerbose)
	progress_bar.done();


  return EXIT_SUCCESS;
}