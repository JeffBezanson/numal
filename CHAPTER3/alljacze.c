#include "../real.h"


void alljaczer(int n, real_t alfa, real_t beta, real_t zer[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void allzerortpol(int, real_t [], real_t [], real_t [], real_t []);
	int i,m;
	real_t sum,min,*a,*b,em[6],gamma,zeri;

	if (alfa == beta) {
		a=allocate_real_vector(0,n/2);
		b=allocate_real_vector(0,n/2);
		m=n/2;
		if (n != 2*m) {
			gamma=0.5;
			zer[m+1]=0.0;
		} else
			gamma = -0.5;
		min=0.25-alfa*alfa;
		sum=alfa+gamma+2.0;
		a[0]=(gamma-alfa)/sum;
		a[1]=min/sum/(sum+2.0);
		b[1]=4.0*(1.0+alfa)*(1.0+gamma)/sum/sum/(sum+1.0);
		for (i=2; i<=m-1; i++) {
			sum=i+i+alfa+gamma;
			a[i]=min/sum/(sum+2.0);
			sum *= sum;
			b[i]=4.0*i*(i+alfa+gamma)*(i+alfa)*(i+gamma)/sum/(sum-1.0);
		}
		em[0]=FLT_MIN;
		em[2]=FLT_EPSILON;
		em[4]=6*m;
		allzerortpol(m,a,b,zer,em);
		for (i=1; i<=m; i++) {
			zer[i] = zeri = -sqrt((1.0+zer[i])/2.0);
			zer[n+1-i] = -zeri;
		}
		free_real_vector(a,0);
		free_real_vector(b,0);
	} else {
		a=allocate_real_vector(0,n);
		b=allocate_real_vector(0,n);
		min=(beta-alfa)*(beta+alfa);
		sum=alfa+beta+2.0;
		b[0]=0.0;
		a[0]=(beta-alfa)/sum;
		a[1]=min/sum/(sum+2.0);
		b[1]=4.0*(1.0+alfa)*(1.0+beta)/sum/sum/(sum+1.0);
		for (i=2; i<=n-1; i++) {
			sum=i+i+alfa+beta;
			a[i]=min/sum/(sum+2.0);
			sum *= sum;
			b[i]=4.0*i*(i+alfa+beta)*(i+alfa)*(i+beta)/(sum-1.0)/sum;
		}
		em[0]=FLT_MIN;
		em[2]=FLT_EPSILON;
		em[4]=6*n;
		allzerortpol(n,a,b,zer,em);
		free_real_vector(a,0);
		free_real_vector(b,0);
	}
}
