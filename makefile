FLAGS = -o -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

start: a.out
%.o: main%.o
	mpic++ $(FLAGS) $< -o -O3 $@
%.o: %.cpp 
	mpic++ -c $(FLAGS) $< -o $@

main.o: main.cpp
matrix.o: matrix.cpp matrix.h

a.out: main.o matrix.o
	mpic++ main.o matrix.o -o a.out

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.exe
