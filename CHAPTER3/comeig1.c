#include "../real.h"
int comeig1(real_t **a, int n, real_t em[], real_t re[], real_t im[],
				real_t **vec)
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	void eqilbr(real_t **, int, real_t [], real_t [], int []);
	void tfmreahes(real_t **, int, real_t [], int []);
	void bakreahes2(real_t **, int, int, int, int [], real_t **);
	void baklbr(int, int, int, real_t [], int [], real_t **);
	void reaveches(real_t **, int, real_t, real_t [], real_t []);
	void comscl(real_t **, int, int, int, real_t []);
	int comvalqri(real_t **, int, real_t [], real_t [], real_t []);
	void comveches(real_t **, int, real_t, real_t,
						real_t [], real_t [], real_t []);
	int i,j,k,pj,itt,again,*ind,*ind0;
	real_t x,y,max,neps,**ab,*d,*u,*v,temp1,temp2;

	ind=allocate_integer_vector(1,n);
	ind0=allocate_integer_vector(1,n);
	d=allocate_real_vector(1,n);
	u=allocate_real_vector(1,n);
	v=allocate_real_vector(1,n);
	ab=allocate_real_matrix(1,n,1,n);

	eqilbr(a,n,em,d,ind0);
	tfmreahes(a,n,em,ind);
	for (i=1; i<=n; i++)
		for (j=((i == 1) ? 1 : i-1); j<=n; j++) ab[i][j]=a[i][j];
	k=comvalqri(ab,n,em,re,im);
	neps=em[0]*em[1];
	max=0.0;
	itt=0;
	for (i=k+1; i<=n; i++) {
		x=re[i];
		y=im[i];
		pj=0;
		again=1;
		do {
			for (j=k+1; j<=i-1; j++) {
				temp1=x-re[j];
				temp2=y-im[j];
				if (temp1*temp1+temp2*temp2 <= neps*neps) {
					if (pj == j)
						neps=em[2]*em[1];
					else
						pj=j;
					x += 2.0*neps;
					again = (!again);
					break;
				}
			}
			again = (!again);
		} while (again);
		re[i]=x;
		for (i=1; i<=n; i++)
			for (j=((i == 1) ? 1 : i-1); j<=n; j++) ab[i][j]=a[i][j];
		if (y != 0.0) {
			comveches(ab,n,re[i],im[i],em,u,v);
			for (j=1; j<=n; j++) vec[j][i]=u[j];
			i++;
			re[i]=x;
		} else
			reaveches(ab,n,x,em,v);
		for (j=1; j<=n; j++) vec[j][i]=v[j];
		if (em[7] > max) max=em[7];
		if (itt < em[9]) itt=em[9];
	}
	em[7]=max;
	em[9]=itt;
	bakreahes2(a,n,k+1,n,ind,vec);
	baklbr(n,k+1,n,d,ind0,vec);
	comscl(vec,n,k+1,n,im);
	free_integer_vector(ind,1);
	free_integer_vector(ind0,1);
	free_real_vector(d,1);
	free_real_vector(u,1);
	free_real_vector(v,1);
	free_real_matrix(ab,1,n,1);
	return k;
}
