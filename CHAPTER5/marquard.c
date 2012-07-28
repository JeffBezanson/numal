#include "../real.h"


void marquardt(int m, int n, real_t par[], real_t g[], real_t **v,
					int (*funct)(int, int, real_t[], real_t[]),
					void (*jacobian)(int, int, real_t[], real_t[], real_t **),
					real_t in[], real_t out[])
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	void mulcol(int, int, int, int, real_t **, real_t **, real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	real_t mattam(int, int, int, int, real_t **, real_t **);
	int qrisngvaldec(real_t **, int, int, real_t [], real_t **, real_t []);
	int maxfe,fe,it,i,j,err,emergency;
	real_t vv,ww,w,mu,res,fpar,fparpres,lambda,lambdamin,p,pw,reltolres,
			abstolres,em[8],*val,*b,*bb,*parpres,**jac,temp;

	val=allocate_real_vector(1,n);
	b=allocate_real_vector(1,n);
	bb=allocate_real_vector(1,n);
	parpres=allocate_real_vector(1,n);
	jac=allocate_real_matrix(1,m,1,n);

	vv=10.0;
	w=0.5;
	mu=0.01;
	ww = (in[6] < 1.0e-7) ? 1.0e-8 : 1.0e-1*in[6];
	em[0]=em[2]=em[6]=in[0];
	em[4]=10*n;
	reltolres=in[3];
	abstolres=in[4]*in[4];
	maxfe=in[5];
	err=0;
	fe=it=1;
	p=fpar=res=0.0;
	pw = -log(ww*in[0])/2.30;
	if (!(*funct)(m,n,par,g)) {
		err=3;
		out[4]=fe;
		out[5]=it-1;
		out[1]=err;
		free_real_vector(val,1);
		free_real_vector(b,1);
		free_real_vector(bb,1);
		free_real_vector(parpres,1);
		free_real_matrix(jac,1,m,1);
		return;
	}
	fpar=vecvec(1,m,0,g,g);
	out[3]=sqrt(fpar);
	emergency=0;
	it=1;
	do {
		(*jacobian)(m,n,par,g,jac);
		i=qrisngvaldec(jac,m,n,val,v,em);
		if (it == 1)
			lambda=in[6]*vecvec(1,n,0,val,val);
		else
			if (p == 0.0) lambda *= w;
		for (i=1; i<=n; i++) b[i]=val[i]*tamvec(1,m,i,jac,g);
		while (1) {
			for (i=1; i<=n; i++) bb[i]=b[i]/(val[i]*val[i]+lambda);
			for (i=1; i<=n; i++) parpres[i]=par[i]-matvec(1,n,i,v,bb);
			fe++;
			if (fe >= maxfe)
				err=1;
			else
				if (!(*funct)(m,n,parpres,g)) err=2;
			if (err != 0) {
				emergency=1;
				break;
			}
			fparpres=vecvec(1,m,0,g,g);
			res=fpar-fparpres;
			if (res < mu*vecvec(1,n,0,b,bb)) {
				p += 1.0;
				lambda *= vv;
				if (p == 1.0) {
					lambdamin=ww*vecvec(1,n,0,val,val);
					if (lambda < lambdamin) lambda=lambdamin;
				}
				if (p >= pw) {
					err=4;
					emergency=1;
					break;
				}
			} else {
				dupvec(1,n,0,par,parpres);
				fpar=fparpres;
				break;
			}
		}
		if (emergency) break;
		it++;
	} while (fpar > abstolres && res > reltolres*fpar+abstolres);
	for (i=1; i<=n; i++) mulcol(1,n,i,i,jac,v,1.0/(val[i]+in[0]));
	for (i=1; i<=n; i++)
		for (j=1; j<=i; j++) v[i][j]=v[j][i]=mattam(1,n,i,j,jac,jac);
	lambda=lambdamin=val[1];
	for (i=2; i<=n; i++)
		if (val[i] > lambda)
			lambda=val[i];
		else
			if (val[i] < lambdamin) lambdamin=val[i];
	temp=lambda/(lambdamin+in[0]);
	out[7]=temp*temp;
	out[2]=sqrt(fpar);
	out[6]=sqrt(res+fpar)-out[2];
	out[4]=fe;
	out[5]=it-1;
	out[1]=err;
	free_real_vector(val,1);
	free_real_vector(b,1);
	free_real_vector(bb,1);
	free_real_vector(parpres,1);
	free_real_matrix(jac,1,m,1);
}
