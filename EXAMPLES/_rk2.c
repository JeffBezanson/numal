#define real_t double

#include <math.h>
#include <stdio.h>

real_t fxyz(real_t x, real_t y, real_t z)
{
	return 10.0*(1.0-y*y)*z-y;
}

void main ()
{
	void rk2(real_t *, real_t, real_t, real_t *, real_t, real_t *, real_t,
			real_t (*)(real_t, real_t, real_t), real_t [], real_t [], int);
	int i,fi;
	real_t x,y,z,e[5],d[6],
			b[4]={9.32386578, 18.86305405, 28.40224162, 37.94142918};

	e[1]=e[2]=e[3]=e[4]=e[5]=1.0e-8;
	printf("RK2 delivers :\n");
	for (i=0; i<=3; i++) {
		fi=(b[i] < 10.0);
		rk2(&x,0.0,b[i],&y,2.0,&z,0.0,fxyz,e,d,fi);
		printf("  X = %+e  Y = %+e    DY/DX = %+e\n",x,y,z);
	}
}

