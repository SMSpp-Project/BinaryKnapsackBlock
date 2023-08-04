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
#   Output: $(BKBkOBJ)  = the final object(s) / library                      #
#           $(BKBkH)    = the .h files to include                            #
#           $(BKBkINC)  = the -I$( source directory )                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                            Federica Di Pasquale                            #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BKBkOBJ = $(BKBkSDR)/obj/BinaryKnapsackBlock.o \
          $(BKBkSDR)/obj/DPBinaryKnapsackSolver.o  

BKBkINC = -I$(BKBkSDR)/include

BKBkH   = $(BKBkSDR)/include/BinaryKnapsackBlock.h \
          $(BKBkSDR)/include/DPBinaryKnapsackSolver.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BKBkOBJ) $(BKBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BKBkSDR)/obj/BinaryKnapsackBlock.o: $(BKBkSDR)/src/BinaryKnapsackBlock.cpp \
	$(BKBkSDR)/include/BinaryKnapsackBlock.h $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)/src/BinaryKnapsackBlock.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)/obj/DPBinaryKnapsackSolver.o: \
	$(BKBkSDR)/src/DPBinaryKnapsackSolver.cpp $(BKBkH) $(SMS++OBJ)  
	$(CC) -c $(BKBkSDR)/src/DPBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
