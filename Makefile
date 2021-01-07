# Basic Makefile

### Compilers
CC  = gcc
CXX = g++

DEBUG_LEVEL    = -g
EXTRA_CCFLAGS  = -W -Wall -std=c++11
CPPFLAGS       = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS        = $(CPPFLAGS)

RM = rm -f
MV = mv

### ROOT
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

### RAT
RATLIBS  := -L$(RATROOT)/lib -lRATEvent

### BOOST
BOOSTCFLAGS := -I/data/snoplus/home/zsoldos/.local/boost-1.71.0
BOOSTLIBS   := -L/data/snoplus/home/zsoldos/.local/boost-1.71.0/lib -lboost_system -lboost_filesystem

### NLOPT
NLOPTCFLAGS := -I/data/snoplus/home/zsoldos/.local/nlopt-2.6.2-install/include
NLOPTLIBS   := -L/data/snoplus/home/zsoldos/.local/nlopt-2.6.2-install/lib -lnlopt -lm

CPPFLAGS  += -Iinclude -IwRATter/include $(ROOTCFLAGS) -I$(RATROOT)/include
CPPFLAGS  +=  $(BOOSTCFLAGS)
CPPFLAGS  +=  $(NLOPTCFLAGS)

EXTRALIBS  = $(ROOTLIBS)
EXTRALIBS += $(RATLIBS)
EXTRALIBS += $(BOOSTLIBS)
EXTRALIBS += $(NLOPTLIBS)

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(subst .cc,.o,$(SRCS))

.PHONY: all clean 
.DEFAULT_GOAL = CreatePDF

help:
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: CreatePDF

CreatePDF: 
	$(CXX) $(CPPFLAGS) -o CreatePDF CreatePDF.cc $(OBJS) $(EXTRALIBS)

SeedNDestroy: 
	$(CXX) $(CPPFLAGS) -o SeedNDestroy SeedNDestroy.cc $(OBJS) $(EXTRALIBS)

clean:
	$(RM) $(OBJS) CreatePDF SeedNDestroy
