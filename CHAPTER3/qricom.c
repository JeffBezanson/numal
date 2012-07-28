#include "../real.h"


int qricom(real_t **a1, real_t **a2, real_t b[], int n, real_t em[],
				real_t val1[], real_t val2[], real_t **vec1, real_t **vec2)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void comkwd(real_t, real_t, real_t, real_t,
					real_t *, real_t *, real_t *, real_t *);
	void rotcomrow(int, int, int, int, real_t **, real_t **,
						real_t, real_t, real_t);
	void rotcomcol(int, int, int, int, real_t **, real_t **,
						real_t, real_t, real_t);
	void comcolcst(int, int, int, real_t **, real_t **, real_t, real_t);
	void comrowcst(int, int, int, real_t **, real_t **, real_t, real_t);
	real_t matvec(int, int, int, real_t **, real_t []);
	void commatvec(int, int, int, real_t **, real_t **,
						real_t [], real_t [], real_t *, real_t *);
	void comdiv(real_t, real_t, real_t, real_t, real_t *, real_t *);
	int m,nm1,i,i1,j,q,q1,max,count;
	real_t r,z1,z2,dd1,dd2,cc,p1,p2,t1,t2,delta1,delta2,mv1,mv2,h,h1,
			h2,g1,g2,k1,k2,hc,aij12,aij22,a1nn,a2nn,aij1,aij2,ai1i,
			kappa,nui,mui1,mui2,muim11,muim12,nuim1,tol,machtol,
			*tf1,*tf2;

	tf1=allocate_real_vector(1,n);
	tf2=allocate_real_vector(1,n);
	tol=em[1]*em[2];
	machtol=em[0]*em[1];
	max=em[4];
	count=0;
	r=0.0;
	m=n;
	if (n > 1) hc=b[n-1];
	for (i=1; i<=n; i++) {
		vec1[i][i]=1.0;
		vec2[i][i]=0.0;
		for (j=i+1; j<=n; j++)
			vec1[i][j]=vec1[j][i]=vec2[i][j]=vec2[j][i]=0.0;
	}
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
			p1=(a1[nm1][nm1]-dd1)*0.5;
			p2=(a2[nm1][nm1]-dd2)*0.5;
			comkwd(p1,p2,cc*a1[nm1][n],cc*a2[nm1][n],&g1,&g2,&k1,&k2);
			if (q == nm1) {
				a1[n][n]=val1[n]=g1+dd1;
				a2[n][n]=val2[n]=g2+dd2;
				a1[q][q]=val1[q]=k1+dd1;
				a2[q][q]=val2[q]=k2+dd2;
				kappa=sqrt(k1*k1+k2*k2+cc*cc);
				nui=cc/kappa;
				mui1=k1/kappa;
				mui2=k2/kappa;
				aij1=a1[q][n];
				aij2=a2[q][n];
				h1=mui1*mui1-mui2*mui2;
				h2=2.0*mui1*mui2;
				h = -nui*2.0;
				a1[q][n]=h*(p1*mui1+p2*mui2)-nui*nui*cc+aij1*h1+aij2*h2;
				a2[q][n]=h*(p2*mui1-p1*mui2)+aij2*h1-aij1*h2;
				rotcomrow(q+2,m,q,n,a1,a2,mui1,mui2,nui);
				rotcomcol(1,q-1,q,n,a1,a2,mui1,-mui2,-nui);
				rotcomcol(1,m,q,n,vec1,vec2,mui1,-mui2,-nui);
				n -= 2;
				if (n > 1) hc=b[n-1];
				b[q]=0.0;
			} else {
				count++;
				if (count > max) {
					em[3]=r;
					em[5]=count;
					free_real_vector(tf1,1);
					free_real_vector(tf2,1);
					return n;
				}
				z1=k1+dd1;
				z2=k2+dd2;
				if (fabs(cc) > fabs(hc)) z1 += fabs(cc);
				hc=cc/2.0;
				q1=q+1;
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
				rotcomrow(q1,m,q,q1,a1,a2,mui1,mui2,nui);
				rotcomcol(1,q,q,q1,a1,a2,mui1,-mui2,-nui);
				a1[q][q] += z1;
				a2[q][q] += z2;
				rotcomcol(1,m,q,q1,vec1,vec2,mui1,-mui2,-nui);
				for (i=q1; i<=nm1; i++) {
					i1=i+1;
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
					rotcomrow(i1,m,i,i1,a1,a2,mui1,mui2,nui);
					a1[i][i]=muim11*kappa;
					a2[i][i] = -muim12*kappa;
					b[i-1]=nuim1*kappa;
					rotcomcol(1,i,i,i1,a1,a2,mui1,-mui2,-nui);
					a1[i][i] += z1;
					a2[i][i] += z2;
					rotcomcol(1,m,i,i1,vec1,vec2,mui1,-mui2,-nui);
				}
				aij1=a1[n][n];
				aij2=a2[n][n];
				aij12=aij1*aij1;
				aij22=aij2*aij2;
				kappa=sqrt(aij12+aij22);
				if ((kappa < tol) ? 1 : (aij22 <= em[0]*aij12)) {
					b[nm1]=nui*aij1;
					a1[n][n]=aij1*mui1+z1;
					a2[n][n] = -aij1*mui2+z2;
				} else {
					b[nm1]=nui*kappa;
					a1nn=mui1*kappa;
					a2nn = -mui2*kappa;
					mui1=aij1/kappa;
					mui2=aij2/kappa;
					comcolcst(1,nm1,n,a1,a2,mui1,mui2);
					comcolcst(1,nm1,n,vec1,vec2,mui1,mui2);
					comrowcst(n+1,m,n,a1,a2,mui1,-mui2);
					comcolcst(n,m,n,vec1,vec2,mui1,mui2);
					a1[n][n]=mui1*a1nn-mui2*a2nn+z1;
					a2[n][n]=mui1*a2nn+mui2*a1nn+z2;
				}
			}
		}
	} while (n > 0);
	for (j=m; j>=2; j--) {
		tf1[j]=1.0;
		tf2[j]=0.0;
		t1=a1[j][j];
		t2=a2[j][j];
		for (i=j-1; i>=1; i--) {
			delta1=t1-a1[i][i];
			delta2=t2-a2[i][i];
			commatvec(i+1,j,i,a1,a2,tf1,tf2,&mv1,&mv2);
			if (fabs(delta1) < machtol && fabs(delta2) < machtol) {
				tf1[i]=mv1/machtol;
				tf2[i]=mv2/machtol;
			} else
				comdiv(mv1,mv2,delta1,delta2,&tf1[i],&tf2[i]);
		}
		for (i=1; i<=m; i++)
			commatvec(1,j,i,vec1,vec2,tf1,tf2,&vec1[i][j],&vec2[i][j]);
	}
	em[3]=r;
	em[5]=count;
	free_real_vector(tf1,1);
	free_real_vector(tf2,1);
	return n;
}

