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

all: RunMMAnalysisTemplate.x RunMMAnalysis.x RunMMAnalysis_v2.x MakeClusterTree.x RunAlignment.x HitEffAnalysis.x ResolutionAnalysis.x TPAnalysis.x GBTAnalysis.x

RunMMAnalysisTemplate.x:  $(SRCDIR)RunMMAnalysisTemplate.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o RunMMAnalysisTemplate.x $(GLIBS) $ $<
	touch RunMMAnalysisTemplate.x

RunMMAnalysis.x:  $(SRCDIR)RunMMAnalysis.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o RunMMAnalysis.x $(GLIBS) $ $<
	touch RunMMAnalysis.x

RunMMAnalysis_v2.x:  $(SRCDIR)RunMMAnalysis_v2.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o RunMMAnalysis_v2.x $(GLIBS) $ $<
	touch RunMMAnalysis_v2.x

HitEffAnalysis.x:  $(SRCDIR)HitEffAnalysis.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o HitEffAnalysis.x $(GLIBS) $ $<
	touch HitEffAnalysis.x

ResolutionAnalysis.x:  $(SRCDIR)ResolutionAnalysis.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o ResolutionAnalysis.x $(GLIBS) $ $<
	touch ResolutionAnalysis.x

TPAnalysis.x:  $(SRCDIR)TPAnalysis.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o TPAnalysis.x $(GLIBS) $ $<
	touch TPAnalysis.x

GBTAnalysis.x:  $(SRCDIR)GBTAnalysis.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o GBTAnalysis.x $(GLIBS) $ $<
	touch GBTAnalysis.x

MakeClusterTree.x:  $(SRCDIR)MakeClusterTree.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o MakeClusterTree.x $(GLIBS) $ $<
	touch MakeClusterTree.x

RunAlignment.x:  $(SRCDIR)RunAlignment.C $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o RunAlignment.x $(GLIBS) $ $<
	touch RunAlignment.x
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -rf *.dSYM
