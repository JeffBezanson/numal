#include <stdio.h>
void main ()
{
	void allortpol(int, real_t, real_t [], real_t [], real_t []);
	int i;
	real_t b[5],c[5],p[6];

	b[0]=1.0;
	for (i=1; i<=4; i++) {
		b[i]=2*i+1;
		c[i]=i*i;
	}
	allortpol(5,0.0,b,c,p);
	printf("ALLORTPOL delivers:  %-6.1f%-6.1f%-6.1f%-6.1f%-6.1f%-6.1f",
			p[0],p[1],p[2],p[3],p[4],p[5]);
}

