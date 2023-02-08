# SnD

A general framework for building and deploying applications looping through a collection of events.
The goal is to provide to the user a way to decouple the I/O from the analysis. 

The file processing is done in a streaming fashion, so that the user can process a large amount of data without having to load it all in memory.
The include `<Templates/TReader.hh>` expresses the architecture to follow for the I/O:
```c++
class TReader {
 protected:
  virtual bool GetNextEvent() = 0;
  virtual bool GetNextTrigger() = 0;
  virtual void *GetData() = 0;
  virtual ProgressBar *GetProgressBar() = 0;
  virtual bool GetVerbosity() = 0;
 public:
  static volatile sig_atomic_t gSignalStatus;
  static void MyHandler(int sig){
	TReader::gSignalStatus = sig;
  };
  virtual void Read(TAnalysis *Ana){
	signal(SIGINT, MyHandler);
	while(this->GetNextEvent()){
	  if(TReader::gSignalStatus == SIGINT)
		break;
	  ++(*GetProgressBar());
	  while(this->GetNextTrigger()){
		Ana->Do(this->GetData());
	  }
	  if(GetVerbosity())
		GetProgressBar()->display();
	}
	GetProgressBar()->done();
  }
};
```
Afterwords, the user can implement the `TReader` class for the specific I/O, and the `TAnalysis` class for the analysis.
The include `<Templates/TAnalysis.hh>` is straightforward and in the only requirements to implement any analysis.
```c++
class TAnalysis {
 protected:
 public:
  virtual void Do(void *Data) = 0;
```

## Full example

Check `src/Readers/NTuple.hh` a flat ROOT tree reader

## Dependencies

`SnD` is based on C++11 for compatibility with Universities HPCs.
There is no specific dependencies to derive the templates and implement your own. 
However, the current apps implemented uses the following libraries:
```bash
ROOT  6.24.26
BOOST 1.70.00
```
If you are using the [RAT]() package, `SnD` is compatible with the tag version `theia/v1.0`.

## Environment
Is you want to use `RAT` with `SnD`:
```bash
$> export ROOT_INCLUDE_PATH=${RATROOT}/include
```
Be sure to have `ROOT` and `BOOST` in your cmake prefix path (it should be for any global installation).
For `RAT`, add when you do cmake:
```bash
-DCMAKE_PREFIX_PATH=/path/to/rat
```

## Build
```bash
$> git clone --recurse-submodules git@github.com:P3tru/SeedNDestroy.git
$> cd SeedNDestroy
$> mkdir -p build
$> cd build
$> cmake ../
$> make
```

## Usage
Verbosity level is controlled with flag -v

### Create PDF
```bash
$> ./TApp -v -i IN.root -o PDF.root
```

### Recon
```bash
$> ./ReconApp -v -r <R> -hh <HH> -p PDF.root -i IN.root -o OUT.root
```