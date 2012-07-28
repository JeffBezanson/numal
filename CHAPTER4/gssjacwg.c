#include "../real.h"


void gssjacwghts(int n, real_t alfa, real_t beta, real_t x[], real_t w[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t gamma(real_t);
	void alljaczer(int, real_t, real_t, real_t []);
	int i,j,m;
	real_t r0,r1,r2,s,h0,alfa2,xi,*a,*b,min,sum,alfabeta,temp;

	if (alfa == beta) {
		b=allocate_real_vector(1,n-1);
		alljaczer(n,alfa,alfa,x);
		alfa2=2.0*alfa;
		temp=gamma(1.0+alfa);
		h0=pow(2.0,alfa2+1.0)*temp*temp/gamma(alfa2+2.0);
		b[1]=1.0/sqrt(3.0+alfa2);
		m=n-n/2;
		for (i=2; i<=n-1; i++)
			b[i]=sqrt(i*(i+alfa2)/(4.0*(i+alfa)*(i+alfa)-1.0));
		for (i=1; i<=m; i++) {
			xi=fabs(x[i]);
			r0=1.0;
			r1=xi/b[1];
			s=1.0+r1*r1;
			for (j=2; j<=n-1; j++) {
				r2=(xi*r1-b[j-1]*r0)/b[j];
				r0=r1;
				r1=r2;
				s += r2*r2;
			}
			w[i]=w[n+1-i]=h0/s;
		}
		free_real_vector(b,1);
	} else {
		a=allocate_real_vector(0,n);
		b=allocate_real_vector(0,n);
		alfabeta=alfa+beta;
		min=(beta-alfa)*alfabeta;
		b[0]=0.0;
		sum=alfabeta+2.0;
		a[0]=(beta-alfa)/sum;
		a[1]=min/sum/(sum+2.0);
		b[1]=2.0*sqrt((1.0+alfa)*(1.0+beta)/(sum+1.0))/sum;
		for (i=2; i<=n-1; i++) {
			sum=i+i+alfabeta;
			a[i]=min/sum/(sum+2.0);
			b[i]=(2.0/sum)*sqrt(i*(sum-i)*(i+alfa)*(i+beta)/(sum*sum-1.0));
		}
		h0=pow(2.0,alfabeta+1.0)*gamma(1.0+alfa)*gamma(1.0+beta)/
					gamma(2.0+alfabeta);
		alljaczer(n,alfa,beta,x);
		for (i=1; i<=n; i++) {
			xi=x[i];
			r0=1.0;
			r1=(xi-a[0])/b[1];
			sum=1.0+r1*r1;
			for (j=2; j<=n-1; j++) {
				r2=((xi-a[j-1])*r1-b[j-1]*r0)/b[j];
				sum += r2*r2;
				r0=r1;
				r1=r2;
			}
			w[i]=h0/sum;
		}
		free_real_vector(a,0);
		free_real_vector(b,0);
	}
}
