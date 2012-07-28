#include "../real.h"


void minmaxpol(int n, int m, real_t y[], real_t fy[], real_t co[],
					real_t em[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	real_t pol(int, real_t, real_t []);
	void newton(int, real_t [], real_t []);
	void newgrn(int, real_t [], real_t []);
	void ini(int, int, int []);
	void sndremez(int, int, int [], real_t [], real_t []);
	int np1,k,pomk,count,cnt,j,mi,*s,sjm1,sj,s0,up;
	real_t e,abse,abseh,*x,*b,*coef,*g;

	s=allocate_integer_vector(0,n+1);
	x=allocate_real_vector(0,n+1);
	b=allocate_real_vector(0,n+1);
	coef=allocate_real_vector(0,n+1);
	g=allocate_real_vector(0,m);

	np1=n+1;
	ini(np1,m,s);
	mi=em[2];
	abse=0.0;
	count=1;
	do {
		pomk=1;
		for (k=0; k<=np1; k++) {
			x[k]=y[s[k]];
			coef[k]=fy[s[k]];
			b[k]=pomk;
			pomk = -pomk;
		}
		newton(np1,x,coef);
		newton(np1,x,b);
		em[0]=e=coef[np1]/b[np1];
		elmvec(0,n,0,coef,b,-e);
		newgrn(n,x,coef);
		s0=sjm1=s[0];
		g[s0]=e;
		for (j=1; j<=np1; j++) {
			sj=s[j];
			up=sj-1;
			for (k=sjm1+1; k<=up; k++) g[k]=fy[k]-pol(n,y[k],coef);
			g[sj] = e = -e;
			sjm1=sj;
		}
		for (k=s0-1; k>=0; k--) g[k]=fy[k]-pol(n,y[k],coef);
		for (k=sj+1; k<=m; k++) g[k]=fy[k]-pol(n,y[k],coef);
		sndremez(np1,m,s,g,em);
		abseh=abse;
		abse=fabs(e);
		cnt=count;
		count++;
	} while (count <= mi && abse > abseh);
	em[2]=mi;
	em[3]=cnt;
	dupvec(0,n,0,co,coef);
	free_integer_vector(s,0);
	free_real_vector(x,0);
	free_real_vector(b,0);
	free_real_vector(coef,0);
	free_real_vector(g,0);
}
