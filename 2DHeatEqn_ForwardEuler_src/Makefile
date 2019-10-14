2DHeatEquationForwardEuler: 2DHeat_ForwEuler.h AdjointObjective.c AdjointUpdateFunction.c Objective.c UpdateFunction.c
	gcc -O3 -fPIC -c AdjointObjective.c      -o AdjointObjective.o
	gcc -O3 -fPIC -c AdjointUpdateFunction.c -o AdjointUpdateFunction.o
	gcc -O3 -fPIC -c Objective.c             -o Objective.o
	gcc -O3 -fPIC -c UpdateFunction.c        -o UpdateFunction.o
	gcc -shared -O3 AdjointObjective.o AdjointUpdateFunction.o Objective.o UpdateFunction.o -o lib2DHeatEquationForwardEuler.so

clean:
	rm -rf *.o

distclean:
	rm -rf *.so
