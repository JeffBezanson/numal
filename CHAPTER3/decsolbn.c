#include "../real.h"


void decsolbnd(real_t a[], int n, int lw, int rw, real_t aux[], real_t b[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	void ichvec(int, int, int, real_t []);
	int i,j,k,kk,kk1,pk,ik,lw1,f,q,w,w1,w2,iw,nrw,shift,sdet;
	real_t r,s,eps,min,*m,*v;

	m=allocate_real_vector(0,lw);
	v=allocate_real_vector(1,n);
	f=lw;
	sdet=1;
	w1=lw+rw;
	w=w1+1;
	w2=w-2;
	iw=0;
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
	if (f > nrw) w2 += nrw-f;
	for (k=1; k<=n; k++) {
		if (f < n) f++;
		ik = kk += w;
		s=fabs(a[kk])/v[k];
		pk=k;
		kk1=kk+1;
		for (i=k+1; i<=f; i++) {
			ik += w1;
			m[i-k]=r=a[ik];
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
			free_real_vector(m,0);
			free_real_vector(v,1);
			return;
		}
		if (k+w2 >= n) w2--;
		if (pk != k) {
			v[pk]=v[k];
			pk -= k;
			ichvec(kk1,kk1+w2,pk*w1,a);
			sdet = -sdet;
			r=b[k];
			b[k]=b[pk+k];
			b[pk+k]=r;
			r=m[pk];
			m[pk]=a[kk];
			a[kk]=r;
		} else
			r=a[kk];
		iw=kk1;
		lw1=f-k;
		if (r < 0.0) sdet = -sdet;
		for (i=1; i<=lw1; i++) {
			s = m[i] /= r;
			iw += w1;
			elmvec(iw,iw+w2,kk1-iw,a,a,-s);
			b[k+i] -= b[k]*s;
		}
	}
	aux[3]=n;
	aux[5]=min;
	kk=(n+1)*w-w1;
	w2 = -1;
	shift=n*w1;
	for (k=n; k>=1; k--) {
		kk -= w;
		shift -= w1;
		if (w2 < w1) w2++;
		b[k]=(b[k]-vecvec(k+1,k+w2,shift,b,a))/a[kk];
	}
	aux[1]=sdet;
	free_real_vector(m,0);
	free_real_vector(v,1);
}
