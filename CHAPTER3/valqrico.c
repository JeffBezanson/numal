#include "../real.h"


int valqricom(real_t **a1, real_t **a2, real_t b[], int n, real_t em[],
					real_t val1[], real_t val2[])
{
	void comcolcst(int, int, int, real_t **, real_t **, real_t, real_t);
	void rotcomcol(int, int, int, int, real_t **, real_t **,
						real_t, real_t, real_t);
	void rotcomrow(int, int, int, int, real_t **, real_t **,
						real_t, real_t, real_t);
	void comkwd(real_t, real_t, real_t, real_t,
					real_t *, real_t *, real_t *, real_t *);
	int nm1,i,i1,q,q1,max,count;
	real_t r,z1,z2,dd1,dd2,cc,g1,g2,k1,k2,hc,a1nn,a2nn,aij1,aij2,
			ai1i,kappa,nui,mui1,mui2,muim11,muim12,nuim1,tol;

	tol=em[1]*em[2];
	max=em[4];
	count=0;
	r=0.0;
	if (n > 1) hc=b[n-1];
	do {
		nm1=n-1;
		i=n;
		do {
			q=i;
			i--;
		} while ((i >= 1) ? (fabs(b[i]) > tol) : 0);
		if (q > 1)
			if (fabs(b[q-1]) > r) r=fabs(b[q-1]);
		if (q == n) {
			val1[n]=a1[n][n];
			val2[n]=a2[n][n];
			n=nm1;
			if (n > 1) hc=b[n-1];
		} else {
			dd1=a1[n][n];
			dd2=a2[n][n];
			cc=b[nm1];
			comkwd((a1[nm1][nm1]-dd1)/2.0,(a2[nm1][nm1]-dd2)/2.0,
					cc*a1[nm1][n],cc*a2[nm1][n],&g1,&g2,&k1,&k2);
			if (q == nm1) {
				val1[nm1]=g1+dd1;
				val2[nm1]=g2+dd2;
				val1[n]=k1+dd1;
				val2[n]=k2+dd2;
				n -= 2;
				if (n > 1) hc=b[n-1];
			} else {
				count++;
				if (count > max) break;
				z1=k1+dd1;
				z2=k2+dd2;
				if (fabs(cc) > fabs(hc)) z1 += fabs(cc);
				hc=cc/2.0;
				i=q1=q+1;
				aij1=a1[q][q]-z1;
				aij2=a2[q][q]-z2;
				ai1i=b[q];
				kappa=sqrt(aij1*aij1+aij2*aij2+ai1i*ai1i);
				mui1=aij1/kappa;
				mui2=aij2/kappa;
				nui=ai1i/kappa;
				a1[q][q]=kappa;
				a2[q][q]=0.0;
				a1[q1][q1] -= z1;
				a2[q1][q1] -= z2;
				rotcomrow(q1,n,q,q1,a1,a2,mui1,mui2,nui);
				rotcomcol(q,q,q,q1,a1,a2,mui1,-mui2,-nui);
				a1[q][q] += z1;
				a2[q][q] += z2;
				for (i1=q1+1; i1<=n; i1++) {
					aij1=a1[i][i];
					aij2=a2[i][i];
					ai1i=b[i];
					kappa=sqrt(aij1*aij1+aij2*aij2+ai1i*ai1i);
					muim11=mui1;
					muim12=mui2;
					nuim1=nui;
					mui1=aij1/kappa;
					mui2=aij2/kappa;
					nui=ai1i/kappa;
					a1[i1][i1] -= z1;
					a2[i1][i1] -= z2;
					rotcomrow(i1,n,i,i1,a1,a2,mui1,mui2,nui);
					a1[i][i]=muim11*kappa;
					a2[i][i] = -muim12*kappa;
					b[i-1]=nuim1*kappa;
					rotcomcol(q,i,i,i1,a1,a2,mui1,-mui2,-nui);
					a1[i][i] += z1;
					a2[i][i] += z2;
					i=i1;
				}
				aij1=a1[n][n];
				aij2=a2[n][n];
				kappa=sqrt(aij1*aij1+aij2*aij2);
				if ((kappa < tol) ? 1 : (aij2*aij2 <= em[0]*aij1*aij1)) {
					b[nm1]=nui*aij1;
					a1[n][n]=aij1*mui1+z1;
					a2[n][n] = -aij1*mui2+z2;
				} else {
					b[nm1]=nui*kappa;
					a1nn=mui1*kappa;
					a2nn = -mui2*kappa;
					mui1=aij1/kappa;
					mui2=aij2/kappa;
					comcolcst(q,nm1,n,a1,a2,mui1,mui2);
					a1[n][n]=mui1*a1nn-mui2*a2nn+z1;
					a2[n][n]=mui1*a2nn+mui2*a1nn+z2;
				}
			}
		}
	} while (n > 0);
	em[3]=r;
	em[5]=count;
	return n;
}

