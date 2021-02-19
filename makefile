##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#                              Antonio Frangioni                             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

BKBkOBJ = $(BKBkSDR)obj/BinaryKnapsackBlock.o\
 $(BKBkSDR)obj/DPBinaryKnapsackSolver.o 

BKBkINC = -I$(BKBkSDR)/include

BKBkH   = $(BKBkSDR)include/BinaryKnapsackBlock.h\
 $(BKBkSDR)include/DPBinaryKnapsackSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BKBkOBJ) $(BKBkSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BKBkSDR)obj/BinaryKnapsackBlock.o: $(BKBkSDR)src/BinaryKnapsackBlock.cpp \
	$(BKBkSDR)include/BinaryKnapsackBlock.h $(SMS++OBJ)
	$(CC) -c $(BKBkSDR)src/BinaryKnapsackBlock.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

$(BKBkSDR)obj/DPBinaryKnapsackSolver.o:\
 $(BKBkSDR)src/DPBinaryKnapsackSolver.cpp $(BKBkH) \
	$(SMS++OBJ)  
	$(CC) -c $(BKBkSDR)src/DPBinaryKnapsackSolver.cpp -o $@ \
	$(BKBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
