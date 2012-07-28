#include <math.h>
#include <stdio.h>
void main ()
{
	void gsswtssym(int, real_t [], real_t [], real_t []);
	real_t pi,zer[3],w[4],c[5];

	pi=4.0*atan(1.0);
	c[1]=0.5;	c[2]=0.25;	c[3]=0.25;	c[4]=0.25;
	zer[1]=cos(0.9*pi);
	zer[2]=cos(0.7*pi);
	gsswtssym(5,zer,c,w);
	printf("Results:\n %7.3f %7.3f %7.3f %7.3f %7.3f",
				w[1]*pi,w[2]*pi,w[3]*pi,w[2]*pi,w[1]*pi);
}

