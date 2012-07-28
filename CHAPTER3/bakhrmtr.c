#include "../real.h"
void bakhrmtri(real_t **a, int n, int n1, int n2, real_t **vecr,
					real_t **veci, real_t tr[], real_t ti[])
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	void elmcolrow(int, int, int, int, real_t **, real_t **, real_t);
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);
	void comrowcst(int, int, int, real_t **, real_t **, real_t, real_t);
	int i,j,r,rm1;
	real_t c,s,t,qr,qi;

	for (i=1; i<=n; i++)
		for (j=n1; j<=n2; j++) veci[i][j]=0.0;
	c=1.0;
	s=0.0;
	for (j=n-1; j>=1; j--) {
		commul(c,s,tr[j],ti[j],&c,&s);
		comrowcst(n1,n2,j,vecr,veci,c,s);
	}
	rm1=2;
	for (r=3; r<=n; r++) {
		t=a[r][r];
		if (t > 0.0)
			for (j=n1; j<=n2; j++) {
				qr=(tammat(1,rm1,r,j,a,vecr)-matmat(1,rm1,r,j,a,veci))/t;
				qi=(tammat(1,rm1,r,j,a,veci)+matmat(1,rm1,r,j,a,vecr))/t;
				elmcol(1,rm1,j,r,vecr,a,-qr);
				elmcolrow(1,rm1,j,r,vecr,a,-qi);
				elmcolrow(1,rm1,j,r,veci,a,qr);
				elmcol(1,rm1,j,r,veci,a,-qi);
			}
		rm1=r;
	}
}
