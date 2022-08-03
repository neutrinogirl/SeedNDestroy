# SeedNDestroy

## Dependencies
```bash
ROOT  6.24.26
BOOST 1.70.00
NLOPT 2.06.02
RAT   theia/v1.0
```

## Environment
```bash
$> export ROOT_INCLUDE_PATH=${RATROOT}/include
```

## Build
```bash
$> git clone --recurse-submodules -b v2 git@github.com:P3tru/SeedNDestroy.git
$> cd SeedNDestroy
$> mkdir -p build
$> cd build
$> cmake -DNLOPT_INCLUDE_DIRS=/path/to/nlopt/include -DNLOPT_LIB=/path/to/nlopt/lib ../
$> make
```

## Usage
### Create PDF


### Recon