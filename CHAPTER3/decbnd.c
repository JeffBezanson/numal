#include "../real.h"


void decbnd(real_t a[], int n, int lw, int rw, real_t aux[],
				real_t m[], int p[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	void ichvec(int, int, int, real_t []);
	int i,j,k,kk,kk1,pk,mk,ik,lw1,f,q,w,w1,w2,nrw,iw,sdet;
	real_t r,s,eps,min,*v;

	v=allocate_real_vector(1,n);
	f=lw;
	w1=lw+rw;
	w=w1+1;
	w2=w-2;
	iw=0;
	sdet=1;
	nrw=n-rw;
	lw1=lw+1;
	q=lw-1;
	for (i=2; i<=lw; i++) {
		q--;
		iw += w1;
		for (j=iw-q; j<=iw; j++) a[j]=0.0;
	}
	iw = -w2;
	q = -lw;
	for (i=1; i<=n; i++) {
		iw += w;
		if (i <= lw1) iw--;
		q += w;
		if (i > nrw) q--;
		v[i]=sqrt(vecvec(iw,q,0,a,a));
	}
	eps=aux[2];
	min=1.0;
	kk = -w1;
	mk = -lw;
	if (f > nrw) w2 += nrw-f;
	for (k=1; k<=n; k++) {
		if (f < n) f++;
		ik = kk += w;
		mk += lw;
		s=fabs(a[kk])/v[k];
		pk=k;
		kk1=kk+1;
		for (i=k+1; i<=f; i++) {
			ik += w1;
			m[mk+i-k]=r=a[ik];
			a[ik]=0.0;
			r=fabs(r)/v[i];
			if (r > s) {
				s=r;
				pk=i;
			}
		}
		if (s < min) min=s;
		if (s < eps) {
			aux[3]=k-1;
			aux[5]=s;
			aux[1]=sdet;
			free_real_vector(v,1);
			return;
		}
		if (k+w2 >= n) w2--;
		p[k]=pk;
		if (pk != k) {
			v[pk]=v[k];
			pk -= k;
			ichvec(kk1,kk1+w2,pk*w1,a);
			sdet = -sdet;
			r=m[mk+pk];
			m[mk+pk]=a[kk];
			a[kk]=r;
		} else
			r=a[kk];
		if (r < 0.0) sdet = -sdet;
		iw=kk1;
		lw1=f-k+mk;
		for (i=mk+1; i<=lw1; i++) {
			s = m[i] /= r;
			iw += w1;
			elmvec(iw,iw+w2,kk1-iw,a,a,-s);
		}
	}
	aux[3]=n;
	aux[5]=min;
	aux[1]=sdet;
	free_real_vector(v,1);
}
