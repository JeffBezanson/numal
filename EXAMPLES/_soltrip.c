#include <stdio.h>
void main ()
{
	void dectripiv(real_t [], real_t [], real_t [], int,
						real_t [], real_t [], int []);
	void soltripiv(real_t [], real_t [], real_t [], int,
						real_t [], int [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i,piv[30];
	real_t d[31],sub[31],super[31],aid[31],b1[31],b2[31],aux[6];

	for (i=1; i<=30; i++) {
		sub[i]=i*2;
		super[i]=i;
		d[i]=i+10;
		b1[i]=b2[i]=0.0;
	}
	b1[1]=1.0;  b1[2]=12.0;  b1[3]=4.0;
	b2[2]=2.0;  b2[3]=13.0;  b2[4]=6.0;
	aux[2]=1.0e-6;
	dectripiv(sub,d,super,30,aid,aux,piv);
	soltripiv(sub,d,super,30,aid,piv,b1);
	soltripiv(sub,d,super,30,aid,piv,b2);
	b1[2]--;  b2[3]--;
	printf("AUX[3] and AUX[5]:  %e   %e\n"
			"Error in b1:  %e\nError in b2:  %e\n",aux[3],aux[5],
			sqrt(vecvec(1,30,0,b1,b1)),sqrt(vecvec(1,30,0,b2,b2)));
}

