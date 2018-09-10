CXX = g++
#CXXFLAGS = -O0 -g -std=c++11
CXXFLAGS = -O2 -Wall -march=native -std=c++11
INCLUDE = 
TLIB = -lm

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

.SUFFIXES: .cpp

.o :
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) 

#-----File Dependencies----------------------

SRC = $(SRC1) $(SRC2)

# SRC1 = chromosome.cpp dsmga2.cpp fastcounting.cpp global.cpp main.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp mkp.cpp
SRC1 = chromosome.cpp CVRP_Evaluator.cpp dice.cpp fssp_evaluator.cpp Gamma.cpp global.cpp hybrid.cpp LOP.cpp lop_evaluator.cpp main.cpp MallowsModel.cpp mkpath.cpp model_EHBSA.cpp model_Hybrid.cpp model_NHBSA.cpp model_order.cpp Model_PlackettLuce.cpp mt19937ar.cpp myrand.cpp PBP.cpp PermutationTools.cpp PFSP.cpp PlackettLuce.cpp PlackettLuceInference.cpp QAP.cpp qap_evaluator.cpp RankingModel.cpp SortingTools.cpp tinyxml2.cpp Tools.cpp TSP.cpp tsp_evaluator.cpp 

SRC2 = EDA.cc GeneralizedMallowsModel.cc Individual.cc Population.cc Problem.cc Seed.cc Solution.cc Timer.cc

# SRC2 = chromosome.cpp dsmga2.cpp fastcounting.cpp global.cpp sweep.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp mkp.cpp

# SRC3 = genZobrist.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

OBJ1 = $(addsuffix .o, $(basename $(SRC)))
# OBJ2 = $(addsuffix .o, $(basename $(SRC2)))
# OBJ3 = $(addsuffix .o, $(basename $(SRC3)))

all: DSMGA2


DSMGA2: $(OBJ1)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ1)

# sweep: $(OBJ2) DSMGA2
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ2)
# 
# genZobrist: $(OBJ3)
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ3)

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ) DSMGA2

# DO NOT DELETE

# chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
# chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h mkp.h
# dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
# dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h
# dsmga2.o: statistics.h trimatrix.h fastcounting.h mkp.h
# fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
# fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h
# global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
# global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h mkp.h
# main.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
# main.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h zkey.h
# main.o: sat.h trimatrix.h fastcounting.h mkp.h
# myrand.o: myrand.h mt19937ar.h
# spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
# spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h mkp.h
# nk-wa.o: nk-wa.h
# sat.o: sat.h
# mkp.o: mkp.h
# chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
# chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h mkp.h
# dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
# dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h
# dsmga2.o: statistics.h trimatrix.h fastcounting.h mkp.h
# fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
# fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h mkp.h
# global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
# global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h mkp.h
# sweep.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
# sweep.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h
# sweep.o: zkey.h sat.h trimatrix.h fastcounting.h mkp.h
# myrand.o: myrand.h mt19937ar.h
# spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
# spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h mkp.h
# nk-wa.o: nk-wa.h
# sat.o: sat.h
# mkp.o: mkp.h
