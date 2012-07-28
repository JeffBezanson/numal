#include <math.h>
#include <stdio.h>

void der(int r, real_t y[], real_t *delta)
{
	real_t y1,y2;

	y1=y[1];
	y2=y[2];
	y[1] = -1000.0*y1*(y1+y2-1.999987);
	y[2] = -2500.0*y2*(y1+y2-2.0);
}

void jac(int r, real_t **j, real_t y[], real_t *delta)
{
	real_t y1,y2;

	y1=y[1];
	y2=y[2];
	j[1][1] = 1999.987-1000.0*(2.0*y1+y2);
	j[1][2] = -1000.0*y1;
	j[2][1] = -2500.0*y2;
	j[2][2] = 2500.0*(2.0-y1-2.0*y2);
}

void outp(real_t x, real_t xe, int r, real_t y[], real_t delta,
			int n, int jev, int lu)
{
	real_t ye1,ye2;

	if (x == 50.0) {
		ye1=0.5976546988;
		ye2=1.4023434075;
		printf(" X = %2.0f    N = %4d    JEV = %3d    LU = %4d\n"
				" Y1 = %e    REL.ERR. = %5.2e\n"
				" Y2 = %e    REL.ERR. = %5.2e\n",
				x,n,jev,lu,y[1],fabs((y[1]-ye1)/ye1),
				y[2],fabs((y[2]-ye2)/ye2));
	}
}

void main ()
{
	void gms(real_t *, real_t, int, real_t [], real_t, real_t,
			real_t, real_t *, void (*)(int, real_t [], real_t *),
			void (*)(int, real_t **, real_t [], real_t *),
			real_t, real_t, int *, int *, int *, int, int,
			void (*)(real_t, real_t, int, real_t [], real_t,
						int, int, int));
	int n,jev,lu;
	real_t x,y[3],delta;

	printf("The results with GMS are:\n");
	y[1]=y[2]=1.0;
	x=0.0;
	delta=0.0;
	gms(&x,50.0,2,y,0.01,0.001,0.5,&delta,der,jac,1.0e-5,1.0e-5,
			&n,&jev,&lu,0,0,outp);
}

