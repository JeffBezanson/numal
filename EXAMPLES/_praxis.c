#include <stdio.h>

real_t f(int n, real_t x[])
{
	real_t temp;

	temp=x[2]-x[1]*x[1];
	return temp*temp*100.0+(1.0-x[1])*(1.0-x[1]);
}

void main ()
{
	void praxis(int, real_t [], real_t (*)(int, real_t[]),
					real_t [], real_t []);
	real_t x[3],in[10],out[7];

	in[0]=1.0e-6;  in[1]=in[2]=1.0e-6;  in[5]=250.0;
	in[6]=in[7]=in[8]=in[9]=1.0;
	x[1] = -1.2;  x[2]=1.0;
	praxis(2,x,f,in,out);
	if (out[1] == 0.0) printf("Normal Termination\n\n");
	printf("Minimum is  %e\nFor x is  %e  %e\n"
			"The initial function value was  %e\n"
			"The number of function evaluations needed was %4.0f\n"
			"The number of line searches was %4.0f\n"
			"The step size in the last iteration step was  %e\n",
			out[2],x[1],x[2],out[3],out[4],out[5],out[6]);
}

