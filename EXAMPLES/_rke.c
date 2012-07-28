#include <math.h>
#include <stdio.h>

void rhs(int n, real_t t, real_t y[])
{
	real_t xx,yy,zz;

	xx=y[1];
	yy=y[2];
	zz=y[3];
	y[1]=yy-zz;
	y[2]=xx*xx+2.0*yy+4.0*t;
	y[3]=xx*(xx+5.0)+2.0*zz+4.0*t;
}

void info(int n, real_t t, real_t te, real_t y[], real_t data[])
{
	real_t et,t2,aex,aey,aez,rex,rey,rez;

	if (t == te) {
		et=exp(t);
		t2=2.0*t;
		rex = -et*sin(t2);
		aex=rex-y[1];
		rex=fabs(aex/rex);
		rey=et*et*(8.0+2.0*t2-sin(2.0*t2))/8.0-t2-1.0;
		rez=et*(sin(t2)+2.0*cos(t2))+rey;
		aey=rey-y[2];
		rey=fabs(aey/rey);
		aez=rez-y[3];
		rez=fabs(aez/rez);
		printf("\nT = %2.0f\n"
			"Relative and absolute errors in x, y and z:\n"
			"  RE(X)    RE(Y)    RE(Z)    AE(X)    AE(Y)    AE(Z)\n"
			" %7.2e  %7.2e  %7.2e  %7.2e  %7.2e  %7.2e\n"
			"Number of integration steps performed : %3.0f\n"
			"Number of integration steps skipped   : %3.0f\n"
			"Number of integration steps rejected  : %3.0f\n",
			t,rex,rey,rez,fabs(aex),fabs(aey),fabs(aez),
			data[4],data[6],data[5]);
	}
}

void main ()
{
	void rke(real_t *, real_t *, int, real_t [],
				void (*)(int, real_t, real_t[]), real_t [], int,
				void (*)(int, real_t, real_t, real_t [], real_t []));
	real_t t,te,y[4],data[7];

	te=1.0;
	while (1) {
		y[1]=y[2]=0.0;
		y[3]=2.0;
		t=0.0;
		data[1]=data[2]=1.0e-5;
		rke(&t,&te,3,y,rhs,data,1,info);
		if (te != 1.0) break;
		te = -1.0;
	}
}

