#include "../real.h"
void hshcomhes(real_t **ar, real_t **ai, int n, real_t em[], real_t b[],
					real_t tr[], real_t ti[], real_t del[])
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void elmrowcol(int, int, int, int, real_t **, real_t **, real_t);
	void hshcomprd(int, int, int, int, int, real_t **,
						real_t **, real_t **, real_t **, real_t);
	void comcolcst(int, int, int, real_t **, real_t **, real_t, real_t);
	void comrowcst(int, int, int, real_t **, real_t **, real_t, real_t);
	void carpol(real_t, real_t, real_t *, real_t *, real_t *);
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);
	int hshcomcol(int, int, int, real_t **, real_t **, real_t,
						real_t *, real_t *, real_t *, real_t *);
	int r,rm1,i,nm1;
	real_t tol,t,xr,xi;

	nm1=n-1;
	t=em[0]*em[1];
	tol=t*t;
	rm1=1;
	for (r=2; r<=nm1; r++) {
		if (hshcomcol(r,n,rm1,ar,ai,tol,&(b[rm1]),&(tr[r]),&(ti[r]),&t)) {
			for (i=1; i<=n; i++) {
				xr=(matmat(r,n,i,rm1,ai,ai)-matmat(r,n,i,rm1,ar,ar))/t;
				xi=(-matmat(r,n,i,rm1,ar,ai)-matmat(r,n,i,rm1,ai,ar))/t;
				elmrowcol(r,n,i,rm1,ar,ar,xr);
				elmrowcol(r,n,i,rm1,ar,ai,xi);
				elmrowcol(r,n,i,rm1,ai,ar,xi);
				elmrowcol(r,n,i,rm1,ai,ai,-xr);
			}
			hshcomprd(r,n,r,n,rm1,ar,ai,ar,ai,t);
		}
		del[rm1]=t;
		rm1=r;
	}
	if (n > 1) carpol(ar[n][nm1],ai[n][nm1],&(b[nm1]),&(tr[n]),&(ti[n]));
	rm1=1;
	tr[1]=1.0;
	ti[1]=0.0;
	for (r=2; r<=n; r++) {
		commul(tr[rm1],ti[rm1],tr[r],ti[r],&(tr[r]),&(ti[r]));
		comcolcst(1,rm1,r,ar,ai,tr[r],ti[r]);
		comrowcst(r+1,n,r,ar,ai,tr[r],-ti[r]);
		rm1=r;
	}
}
