#include "../real.h"


void sndremez(int n, int m, int s[], real_t g[], real_t em[])
{
	real_t infnrmvec(int, int, int *, real_t []);
	int s0,sn,sjp1,i,j,k,up,low,nm1;
	real_t max,msjp1,hi,hj,he,abse,h,temp1,temp2;

	s0=sjp1=s[0];
	he=em[0];
	low=s0+1;
	max=msjp1=abse=fabs(he);
	nm1=n-1;
	for (j=0; j<=nm1; j++) {
		up=s[j+1]-1;
		h=infnrmvec(low,up,&i,g);
		if (h > max) max=h;
		if (h > abse)
			if (he*g[i] > 0.0) {
				s[j] = (msjp1 < h) ? i : sjp1;
				sjp1=s[j+1];
				msjp1=abse;
			} else {
				s[j]=sjp1;
				sjp1=i;
				msjp1=h;
			}
		else {
			s[j]=sjp1;
			sjp1=s[j+1];
			msjp1=abse;
		}
		he = -he;
		low=up+2;
	}
	sn=s[n];
	s[n]=sjp1;
	hi=infnrmvec(0,s0-1,&i,g);
	hj=infnrmvec(sn+1,m,&j,g);
	if (j > m) j=m;
	if (hi > hj) {
		if (hi > max) max=hi;
		temp1 = (g[i] == 0.0) ? 0.0 : ((g[i] > 0.0) ? 1.0 : -1.0);
		temp2 = (g[s[0]]==0.0) ? 0.0 : ((g[s[0]]>0.0) ? 1.0 : -1.0);
		if (temp1 == temp2) {
			if (hi > fabs(g[s[0]])) {
				s[0]=i;
				if (g[j]/g[s[n]] > 1.0) s[n]=j;
			}
		}
		else {
			if (hi > fabs(g[s[n]])) {
				s[n] = (g[j]/g[s[nm1]] > 1.0) ? j : s[nm1];
				for (k=nm1; k>=1; k--) s[k]=s[k-1];
				s[0]=i;
			}
		}
	} else {
		if (hj > max) max=hj;
		temp1 = (g[j] == 0.0) ? 0.0 : ((g[j] > 0.0) ? 1.0 : -1.0);
		temp2 = (g[s[n]]==0.0) ? 0.0 : ((g[s[n]]>0.0) ? 1.0 : -1.0);
		if (temp1 == temp2) {
			if (hj > fabs(g[s[n]])) {
				s[n]=j;
				if (g[i]/g[s[0]] > 1.0) s[0]=i;
			}
		} else
			if (hj > fabs(g[s[0]])) {
				s[0] = (g[i]/g[s[1]] > 1.0) ? i : s[1];
				for (k=1; k<=nm1; k++) s[k]=s[k+1];
				s[n]=j;
			}
	}
	em[1]=max;
}
