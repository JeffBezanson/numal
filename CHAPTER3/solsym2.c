#include "../real.h"
void solsym2(real_t **a, int n, real_t b[], int p[], real_t detaux[])
{
	real_t matvec(int, int, int, real_t **, real_t []);
	void elmvecrow(int, int, int, real_t [], real_t **, real_t);
	int i,ii,k,ip1,pi,pii;
	real_t det,temp,save;

	i=1;
	while (i < n) {
		ip1=i+1;
		pi=p[i];
		save=b[pi];
		if (p[ip1] > 0) {
			b[pi]=b[i];
			b[i]=save/a[i][i];
			elmvecrow(ip1,n,i,b,a,save);
			i=ip1;
		} else {
			temp=b[i];
			b[pi]=b[ip1];
			det=detaux[ip1];
			b[i]=(temp*a[ip1][ip1]-save*a[i][ip1])/det;
			b[ip1]=(save*a[i][i]-temp*a[i][ip1])/det;
			elmvecrow(i+2,n,i,b,a,temp);
			elmvecrow(i+2,n,ip1,b,a,save);
			i += 2;
		}
	}
	if (i == n) {
		b[i] /= a[i][i];
		i=n-1;
	} else
		i=n-2;
	while (i > 0) {
		if (p[i] == 0)
			ii=i-1;
		else
			ii=i;
		for (k=ii; k<=i; k++) {
			save=b[k];
			save += matvec(i+1,n,k,a,b);
			b[k]=save;
		}
		pii=p[ii];
		b[i]=b[pii];
		b[pii]=save;
		i=ii-1;
	}
}
