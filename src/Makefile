CC=icpc
FLAGS=-std=c++11 -m64 -O2 -qopenmp -xCORE-AVX512# sandybridge
DIR=-I /users/lwang70/Eigen/eigen
 
OBJ=main.o functions.o bitstr.o introsort.o sparse_solver.o dense_algo.o
HEADERS=all.h 

#Rules....
%.o:%.cpp
	$(CC) -c $(FLAGS) $(DIR) -o $@ $< 

sfci: $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(DIR) -o sfci 

$(OBJ):$(HEADERS)

clean:
	rm -f *.o sfci
