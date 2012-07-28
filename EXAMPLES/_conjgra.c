#include <stdio.h>

void a(real_t x[], real_t b[])
{
	int i;

	b[0]=2.0*x[0]-x[1];
	for (i=1; i<=11; i++) b[i] = -x[i-1]+2.0*x[i]-x[i+1];
	b[12]=2.0*x[12]-x[11];
}

int moveon(int iter, real_t norm)
{
	return (iter<20 && norm>1.0e-10);
}

void main ()
{
	void conjgrad(void (*)(real_t [], real_t []), real_t [], real_t [],
					int, int, int (*)(int, real_t), int *, real_t *);
	int it,i;
	real_t no,x[13],b[13];

	for (i=0; i<=12; i++) x[i]=b[i]=0.0;
	b[0]=1.0;  b[12]=4.0;
	conjgrad(a,x,b,0,12,moveon,&it,&no);
	printf("Delivers:\nIT = %3d,   NO = %e\n"
			"\n      X            R\n\n",it,no);
	for (i=0; i<=12; i++)
		printf(" %+e  %+e\n",x[i],b[i]);
}

