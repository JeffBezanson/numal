#include "../real.h"


real_t airyzeros(int n, int d, real_t zai[], real_t vai[])
{
	void airy(real_t,real_t *,real_t *,real_t *,real_t *,real_t *,int);
	int a,found,i;
	real_t c,e,r,zaj,zak,vaj,daj,kaj,zz;

	a=((d == 0) || (d == 2));
	r = (d == 0 || d == 3) ? -1.17809724509617 : -3.53429173528852;
	airy(0.0,&zaj,&vaj,&daj,&kaj,&zz,1);
	for (i=1; i<=n; i++) {
		r += 4.71238898038469;
		zz=r*r;
		zaj = (i == 1 && d == 1) ? -1.01879297 :
				((i == 1 && d == 2) ? -1.17371322 :
				pow(r,0.666666666666667)*
				(a ? -(1.0+(5.0/48.0-(5.0/36.0-(77125.0/82944.0-
				(108056875.0/6967296.0-(162375596875.0/334430208.0)/
				zz)/zz)/zz)/zz)/zz) : -(1.0-(7.0/48.0-(35.0/288.0-
				(181223.0/207360.0-(18683371.0/1244160.0-
				(91145884361.0/191102976.0)/zz)/zz)/zz)/zz)/zz)));
		if (d <= 1.0)
			airy(zaj,&vaj,&daj,&c,&e,&zz,0);
		else
			airy(zaj,&c,&e,&vaj,&daj,&zz,0);
		found=(fabs(a ? vaj : daj) < 1.0e-12);
		while (!found) {
			if (a) {
				kaj=vaj/daj;
				zak=zaj-kaj*(1.0+zaj*kaj*kaj);
			} else {
				kaj=daj/(zaj*vaj);
				zak=zaj-kaj*(1.0+kaj*(kaj*zaj+1.0/zaj));
			}
			if (d <= 1)
				airy(zak,&vaj,&daj,&c,&e,&zz,0);
			else
				airy(zak,&c,&e,&vaj,&daj,&zz,0);
			found=(fabs(zak-zaj) < 1.0e-14*fabs(zak) ||
					fabs(a ? vaj : daj) < 1.0e-12);
			zaj=zak;
		}
		vai[i]=(a ? daj : vaj);
		zai[i]=zaj;
	}
	return zai[n];
}
