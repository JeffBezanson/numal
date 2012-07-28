#include <stdio.h>

real_t rosenbrock(int n, real_t x[], real_t g[])
{
	real_t temp;

	temp=x[2]-x[1]*x[1];
	g[1]=(-temp*400.0+2.0)*x[1]-2.0;
	g[2]=temp*200.0;
	return temp*temp*100.0+(1.0-x[1])*(1.0-x[1]);
}

void main ()
{
	real_t rnk1min(int, real_t [], real_t [], real_t [],
						real_t (*)(int, real_t[], real_t[]),
						real_t [], real_t []);
	real_t flemin(int, real_t [], real_t [], real_t [],
					real_t (*)(int, real_t[], real_t[]), real_t [], real_t []);
	int again;
	real_t f,x[3],g[3],h[4],in[9],out[5];

	in[0]=1.0e-6;  in[1]=in[2]=1.0e-5;  in[3]=1.0e-4;  in[4]=1.0e-5;
	in[5] = -10.0;  in[6]=1.0;  in[7]=100.0;  in[8]=0.01;
	x[1] = -1.2;  x[2]=1.0;
	again=1;
	f=rnk1min(2,x,g,h,rosenbrock,in,out);
	while (1) {
		printf("\nLeast value:  %e\nx:  %e  %e\n"
			"Gradient:  %e  %e\n"
			"Metric:     %e  %e\n                         %e\n"
			"OUT: %e\n     %e\n     %e\n     %e\n     %e\n",
			f,x[1],x[2],g[1],g[2],h[1],h[2],h[3],out[0],out[1],out[2],
			out[3],out[4]);
		if (again) {
			x[1] = -1.2;  x[2]=1.0;
			again=0;
			f=flemin(2,x,g,h,rosenbrock,in,out);
		} else
			break;
	}
}

