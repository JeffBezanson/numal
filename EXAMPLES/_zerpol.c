#include <stdio.h>
void main ()
{
	int zerpol(int, real_t [], real_t [], real_t [], real_t [], real_t []);
	void bounds(int, real_t [], real_t [], real_t [], real_t,
					real_t, real_t [], real_t [], real_t []);
	int i,j;
	real_t a[8],d[8],re[8],im[8],em[5],recentre[8],imcentre[8],bound[8];

	a[7]=1.0;   a[6]=a[5] = -3.0;  a[4]=25.0;  a[3] = -46.0;
	a[2]=38.0;  a[1] = -12.0;      a[0]=0.0;
	em[0]=1.0e-6;  em[1]=40.0;
	i=zerpol(7,a,em,re,im,d);
	printf("Coefficients of polynomial:\n  ");
	for (j=7; j>=0; j--) printf("%5.0f",a[j]);
	printf("\n\nNumber of not found zeros: %2d\n"
			"Fail indication:          %3.0f\n"
			"Number of new starts:     %3.0f\n"
			"Number of iterations:     %3.0f\n\nZeros:\n",
			i,em[2],em[3],em[4]);
	for (j=i+1; j<=7; j++)
		if (im[j] == 0.0)
			printf(" %12.6e\n",re[j]);
		else
			printf(" %12.6e  %12.6e\n",re[j],im[j]);
	if (i == 0) {
		bounds(7,a,re,im,0.0,0.0,recentre,imcentre,bound);
		printf("\nReal and imaginary part of centre + radius\n");
		for (j=1; j<=7; j++)
			printf("  %12.6e  %12.6e  %12.6e\n",
					recentre[j],imcentre[j],bound[j]);
	}
}

