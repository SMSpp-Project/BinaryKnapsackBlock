##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BinaryKnapsackBlock / DPBinaryKnapsackSolver                 #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by BKBkBlock  #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(SMS++INC) = the -I$( core SMS++ directory )                    #
#           $(SMS++OBJ) = the libSMS++ library itself                        #
#           $(BKBkSDR)  = the directory where the source is                  #
#                                                                            #
#   The external RECORD and COMBO solvers are found through $(RECORD_ROOT)   #
#   (a checkout of RECORD) and $(COMBO_ROOT) (the directory with combo.c     #
#   and combo.h), set in extlib/makefile-default-paths-* (overridable via    #
#   extlib/makefile-paths): RECORDBinaryKnapsackSolver is built if RECORD    #
#   is there, COMBOBinaryKnapsackSolver if COMBO is, each from a copy of     #
#   the source made in obj/ [see CMakeLists.txt].                            #
#                                                                            #
#   Output: $(BKBkOBJ)  = the final object(s) / library                      #
#           $(BKBkH)    = the .h files to include                            #
#           $(BKBkINC)  = the -I$( source directory )                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BKBkOBJ = $(BKBkSDR)/obj/BinaryKnapsackBlock.o \
          $(BKBkSDR)/obj/BinaryKnapsackSolver.o \
          $(BKBkSDR)/obj/DPBinaryKnapsackSolver.o \
          $(BKBkSDR)/obj/ParallelDPBinaryKnapsackSolver.o \
          $(BKBkSDR)/obj/CoreDPBinaryKnapsackSolver.o \
          $(BKBkSDR)/obj/GreedyRelaxationBinaryKnapsackSolver.o \
          $(BKBkSDR)/obj/IncrementalGreedyRelaxationBinaryKnapsackSolver.o

BKBkINC = -I$(BKBkSDR)/include

BKBkH   = $(BKBkSDR)/include/BinaryKnapsackBlock.h \
          $(BKBkSDR)/include/BinaryKnapsackSolver.h \
          $(BKBkSDR)/include/DPBinaryKnapsackSolver.h \
          $(BKBkSDR)/include/ParallelDPBinaryKnapsackSolver.h \
          $(BKBkSDR)/include/CoreDPBinaryKnapsackSolver.h \
          $(BKBkSDR)/include/GreedyRelaxationBinaryKnapsackSolver.h \
          $(BKBkSDR)/include/IncrementalGreedyRelaxationBinaryKnapsackSolver.h

# the external solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# non-empty if and only if the source is there
BKBkRECORD := $(wildcard $(subst $\",,$(RECORD_ROOT))/source/RECORD.cpp)
BKBkCOMBO := $(wildcard $(subst $\",,$(COMBO_ROOT))/combo.c)

ifneq ($(BKBkRECORD),)
BKBkOBJ += $(BKBkSDR)/obj/RECORDBridge.o \
           $(BKBkSDR)/obj/RECORDBinaryKnapsackSolver.o
BKBkH += $(BKBkSDR)/include/RECORDBinaryKnapsackSolver.h
endif

ifneq ($(BKBkCOMBO),)
BKBkOBJ += $(BKBkSDR)/obj/COMBOGuard.o \
           $(BKBkSDR)/obj/COMBOBridge.o \
           $(BKBkSDR)/obj/COMBOBinaryKnapsackSolver.o
BKBkH += $(BKBkSDR)/include/COMBOBinaryKnapsackSolver.h
endif

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BKBkOBJ) $(BKBkSDR)/*~
	rm -rf $(BKBkSDR)/obj/record-shim $(BKBkSDR)/obj/combo-shim

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BKBkSDR)/obj/BinaryKnapsackBlock.o: $(BKBkSDR)/src/BinaryKnapsackBlock.cpp \
	$(BKBkSDR)/include/BinaryKnapsackBlock.h $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/BinaryKnapsackBlock.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/BinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/BinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/BinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/DPBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/DPBinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)  
	$(CC) -c $(BKBkSDR)/src/DPBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/ParallelDPBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/ParallelDPBinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/ParallelDPBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/CoreDPBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/CoreDPBinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/CoreDPBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/GreedyRelaxationBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/GreedyRelaxationBinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/GreedyRelaxationBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/IncrementalGreedyRelaxationBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp \
	$(BKBkH) $(SMS++OBJ)
	$(CC) -c \
	$(BKBkSDR)/src/IncrementalGreedyRelaxationBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

# the copy of RECORD without its main(), compiled alone in RECORDBridge.cpp
$(BKBkSDR)/obj/record-shim/RECORD.hpp: $(BKBkRECORD)
	mkdir -p $(BKBkSDR)/obj/record-shim
	sed '/^int main (/,$$d' $< > $@

$(BKBkSDR)/obj/RECORDBridge.o: $(BKBkSDR)/src/RECORDBridge.cpp \
	$(BKBkSDR)/include/RECORDBridge.h $(BKBkSDR)/obj/record-shim/RECORD.hpp
	$(CC) -c $(BKBkSDR)/src/RECORDBridge.cpp -o $@ \
	-I$(BKBkSDR)/obj/record-shim $(BKBkINC) $(SW) -fopenmp-simd

$(BKBkSDR)/obj/RECORDBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/RECORDBinaryKnapsackSolver.cpp $(BKBkH) \
	$(BKBkSDR)/include/RECORDBridge.h $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/RECORDBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

# the copy of COMBO with its products widened to 128-bit integers, compiled
# as C alone in COMBOGuard.c
$(BKBkSDR)/obj/combo-shim/combo.c: $(BKBkCOMBO)
	mkdir -p $(BKBkSDR)/obj/combo-shim
	sed 's/typedef[[:space:]]*double[[:space:]]*prod;/typedef __int128 prod;/' \
	$< > $@
	cp $(dir $(BKBkCOMBO))combo.h $(BKBkSDR)/obj/combo-shim/

$(BKBkSDR)/obj/COMBOGuard.o: $(BKBkSDR)/src/COMBOGuard.c \
	$(BKBkSDR)/obj/combo-shim/combo.c
	$(CC) -x c -c $(BKBkSDR)/src/COMBOGuard.c -o $@ \
	-I$(BKBkSDR)/obj/combo-shim -O3 -DNDEBUG -w

$(BKBkSDR)/obj/COMBOBridge.o: $(BKBkSDR)/src/COMBOBridge.cpp \
	$(BKBkSDR)/include/COMBOBridge.h $(BKBkSDR)/obj/combo-shim/combo.c
	$(CC) -c $(BKBkSDR)/src/COMBOBridge.cpp -o $@ \
	-I$(BKBkSDR)/obj/combo-shim $(BKBkINC) $(SW)

$(BKBkSDR)/obj/COMBOBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/COMBOBinaryKnapsackSolver.cpp $(BKBkH) \
	$(BKBkSDR)/include/COMBOBridge.h $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/COMBOBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
