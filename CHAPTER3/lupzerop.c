#include "../real.h"


void lupzerortpol(int n, int m, real_t b[], real_t c[], real_t zer[],
						real_t em[])
{
	real_t infnrmvec(int, int, int *, real_t []);
	void dupvec(int, int, int, real_t [], real_t []);
	int i,posdef,j,k,t,converge;
	real_t nrm,dlam,eps,delta,e,ep,err,p,q,qp,r,s,tot;

	nrm=fabs(b[0]);
	for (i=1; i<=n-2; i++)
		if (c[i]+fabs(b[i]) > nrm) nrm=c[i]+fabs(b[i]);
	if (n > 1)
		nrm = (nrm+1 >= c[n-1]+fabs(b[n-1])) ? nrm+1.0 :
					(c[n-1]+fabs(b[n-1]));
	em[1]=nrm;
	for (i=n; i>=1; i--) b[i]=b[i-1];
	for (i=n; i>=2; i--) c[i]=c[i-1];
	posdef = (em[6] == 1.0);
	dlam=em[2];
	eps=em[0];
	c[1]=err=q=s=0.0;
	tot=b[1];
	for (i=n; i>=1; i--) {
		p=q;
		q=sqrt(c[i]);
		e=b[i]-p-q;
		if (e < tot) tot=e;
	}
	if (posdef && (tot < 0.0))
		tot=0.0;
	else
		for(i=1; i<=n; i++) b[i] -= tot;
	t=0;
	for (k=1; k<=m; k++) {
		converge=0;
		/* next qr transformation */
		do {
			t++;
			tot += s;
			delta=b[n]-s;
			i=n;
			e=fabs(eps*tot);
			if (dlam < e) dlam=e;
			if (delta <= dlam) {
				converge=1;
				break;
			}
			e=c[n]/delta;
			qp=delta+e;
			p=1.0;
			for (i=n-1; i>=k; i--) {
				q=b[i]-s-e;
				r=q/qp;
				p=p*r+1.0;
				ep=e*r;
				b[i+1]=qp+ep;
				delta=q-ep;
				if (delta <= dlam) {
					converge=1;
					break;
				}
				e=c[i]/q;
				qp=delta+e;
				c[i+1]=qp*ep;
			}
			if (converge) break;
			b[k]=qp;
			s=qp/p;
		} while (tot+s > tot);	/* end of qr transformation */
		if (!converge) {
			/* irregular end of iteration,
				deflate minimum diagonal element */
			s=0.0;
			i=k;
			delta=qp;
			for (j=k+1; j<=n; j++)
				if (b[j] < delta) {
					i=j;
					delta=b[j];
				}
		}
		/* convergence */
		if (i < n) c[i+1]=c[i]*e/qp;
		for (j=i-1; j>=k; j--) {
			b[j+1]=b[j]-s;
			c[j+1]=c[j];
		}
		b[k]=tot;
		c[k] = err += fabs(delta);
	}
	em[5]=t;
	em[3]=infnrmvec(1,m,&t,c);
	dupvec(1,m,0,zer,b);
}
