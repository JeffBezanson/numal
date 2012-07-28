#include "../real.h"
real_t recipgamma(real_t x, real_t *odd, real_t *even)
{
	int i;
	real_t alfa,beta,x2,b[13];

	b[1] = -0.283876542276024;   b[2]  = -0.076852840844786;
	b[3] =  0.001706305071096;   b[4]  =  0.001271927136655;
	b[5] =  0.000076309597586;   b[6]  = -0.000004971736704;
	b[7] = -0.000000865920800;   b[8]  = -0.000000033126120;
	b[9] =  0.000000001745136;   b[10] =  0.000000000242310;
	b[11]=  0.000000000009161;   b[12] = -0.000000000000170;
	x2=x*x*8.0;
	alfa = -0.000000000000001;
	beta=0.0;
	for (i=12; i>=2; i -= 2) {
		beta = -(alfa*2.0+beta);
		alfa = -beta*x2-alfa+b[i];
	}
	*even=(beta/2.0+alfa)*x2-alfa+0.921870293650453;
	alfa = -0.000000000000034;
	beta=0.0;
	for (i=11; i>=1; i -= 2) {
		beta = -(alfa*2.0+beta);
		alfa = -beta*x2-alfa+b[i];
	}
	*odd=(alfa+beta)*2.0;
	return (*odd)*x+(*even);
}
