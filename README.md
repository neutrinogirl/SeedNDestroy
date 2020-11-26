# SeedNDestroy
PathFitter adapted to WbLS reconstruction.


## Set up dependencies
Depends on RAT, ROOT and BOOST. Have $RATROOT, $ROOTSYS env variable defined. 
- First step is to compile the RAT wrapper (wRATter) and the hit structure
```bash
cd wRATter
source SetROOTEnv.sh
cmake .
make
```
