ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX            = g++

CXXFLAGS       = -fPIC -Wall -O3 -g
CXXFLAGS       += $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))

GLIBS          = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))

INCLUDEDIR       = ./include/
SRCDIR           = ./src/
CXX	         += -I$(INCLUDEDIR) -I.
OUTOBJ	         = ./obj/

CC_FILES := $(wildcard src/*.C)
HH_FILES := $(wildcard include/*.hh)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.C=.o)))

all: octuplet.x

octuplet.x:  $(SRCDIR)octuplet_ana.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o octuplet.x $(GLIBS) $ $<
	touch octuplet.x
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -rf *.dSYM
