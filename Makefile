# Basic Makefile

### Compilers
CC  = gcc
CXX = g++

### Machine name
NAME = $(shell uname -n)

DEBUG_LEVEL    = -g
EXTRA_CCFLAGS  = -std=c++11 -Wall
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
ifeq ($(NAME), nino)
BOOSTCFLAGS := -I/data/snoplus/home/zsoldos/.local/boost-1.71.0
BOOSTLIBS   := -L/data/snoplus/home/zsoldos/.local/boost-1.71.0/lib -lboost_system -lboost_filesystem
else
BOOSTCFLAGS := -I/usr/include/boost
BOOSTLIBS   := -lboost_system -lboost_filesystem
endif

### NLOPT
ifeq ($(NAME), nino)
NLOPTCFLAGS := -I/data/snoplus/home/zsoldos/.local/nlopt-2.6.2-install/include
NLOPTLIBS   := -L/data/snoplus/home/zsoldos/.local/nlopt-2.6.2-install/lib -lnlopt -lm
else
# NLOPTCFLAGS := -I/usr/include # No need to add it, g++ will find it
NLOPTLIBS   := -lnlopt -lm
endif

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

DebugRecon: 
	$(CXX) $(CPPFLAGS) -o DebugRecon DebugRecon.cc $(OBJS) $(EXTRALIBS)

clean:
	$(RM) $(OBJS) CreatePDF SeedNDestroy
