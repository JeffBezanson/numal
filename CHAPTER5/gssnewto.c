#include "../real.h"


void gssnewton(int m, int n, real_t par[], real_t rv[], real_t **jjinv,
					int (*funct)(int, int, real_t[], real_t[]),
					void (*jacobian)(int, int, real_t[], real_t[], real_t **),
					real_t in[], real_t out[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void dupvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	void lsqortdec(real_t **, int, int, real_t [], real_t [], int []);
	void lsqsol(real_t **, int, int, real_t [], int [], real_t []);
	void lsqinv(real_t **, int, real_t [], int []);
	int i,j,inr,mit,text,it,itmax,inrmax,tim,feval,fevalmax,conv,
			testthf,dampingon,*ci,fail;
	real_t rho,res1,res2,rn,reltolpar,abstolpar,abstolres,stap,normx,
			**jac,*pr,*aid,*sol,*fu2,aux[6];

	ci=allocate_integer_vector(1,n);
	pr=allocate_real_vector(1,n);
	aid=allocate_real_vector(1,n);
	sol=allocate_real_vector(1,n);
	fu2=allocate_real_vector(1,m);
	jac=allocate_real_matrix(1,m+1,1,n);

	itmax=fevalmax=in[5];
	aux[2]=n*in[0];
	tim=in[7];
	reltolpar=in[1]*in[1];
	abstolpar=in[2]*in[2];
	abstolres=in[4]*in[4];
	inrmax=in[6];
	dupvec(1,n,0,pr,par);
	if (m < n)
		for (i=1; i<=n; i++) jac[m+1][i]=0.0;
	text=4;
	mit=0;
	testthf=1;
	res2=stap=out[5]=out[6]=out[7]=0.0;
	(*funct)(m,n,par,fu2);
	rn=vecvec(1,m,0,fu2,fu2);
	out[3]=sqrt(rn);
	feval=1;
	dampingon=0;
	fail=0;
	it=1;
	do {
		out[5]=it;
		(*jacobian)(m,n,par,fu2,jac);
		if (!testthf) {
			text=7;
			fail=1;
			break;
		}
		lsqortdec(jac,m,n,aux,aid,ci);
		if (aux[3] != n) {
			text=5;
			fail=1;
			break;
		}
		lsqsol(jac,m,n,aid,ci,fu2);
		dupvec(1,n,0,sol,fu2);
		stap=vecvec(1,n,0,sol,sol);
		rho=2.0;
		normx=vecvec(1,n,0,par,par);
		if (stap > reltolpar*normx+abstolpar || it == 1 && stap > 0.0) {
			inr=0;
			do {
				rho /= 2.0;
				if (inr > 0) {
					res1=res2;
					dupvec(1,m,0,rv,fu2);
					dampingon = inr > 1;
				}
				for (i=1; i<=n; i++) pr[i]=par[i]-sol[i]*rho;
				feval++;
				if (!(*funct)(m,n,pr,fu2)) {
					text=6;
					fail=1;
					break;
				}
				res2=vecvec(1,m,0,fu2,fu2);
				conv = inr >= inrmax;
				inr++;
			} while ((inr == 1) ? (dampingon || res2 >= rn) :
						(!conv && (rn <= res1 || res2 < res1)));
			if (fail) break;
			if (conv) {
				mit++;
				if (mit < tim) conv=0;
			} else
				mit=0;
			if (inr > 1) {
				rho *= 2.0;
				elmvec(1,n,0,par,sol,-rho);
				rn=res1;
				if (inr > 2) out[7]=it;
			} else {
				dupvec(1,n,0,par,pr);
				rn=res2;
				dupvec(1,m,0,rv,fu2);
			}
			if (rn <= abstolres) {
				text=1;
				itmax=it;
			} else
				if (conv && inrmax > 0) {
					text=3;
					itmax=it;
				} else
					dupvec(1,m,0,fu2,rv);
		} else {
			text=2;
			rho=1.0;
			itmax=it;
		}
		it++;
	} while (it <= itmax && feval < fevalmax);
	if (!fail) {
		lsqinv(jac,n,aid,ci);
		for (i=1; i<=n; i++) {
			jjinv[i][i]=jac[i][i];
			for (j=i+1; j<=n; j++) jjinv[i][j]=jjinv[j][i]=jac[i][j];
		}
	}
	out[6]=sqrt(stap)*rho;
	out[2]=sqrt(rn);
	out[4]=feval;
	out[1]=text;
	out[8]=aux[3];
	out[9]=aux[5];
	free_integer_vector(ci,1);
	free_real_vector(pr,1);
	free_real_vector(aid,1);
	free_real_vector(sol,1);
	free_real_vector(fu2,1);
	free_real_matrix(jac,1,m+1,1);
}
