HEADERS = double2.h, double3.h, d2func.h, d3func.h, Enfunc.h, Gradfunc.h, Misc.h, Globalvars.h

default: SCsim

SCsim: 
	gcc -O3 d2func.c d3func.c Globalvars.c Enfunc.c Gradfunc.c Misc.c main.c -o SCsim -lm
