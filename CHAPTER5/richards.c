#include "../real.h"


void richardson(real_t **u, int lj, int uj, int ll, int ul,
			int inap, void (*residual)(int, int, int, int, real_t **),
			real_t a, real_t b, int *n, real_t discr[], int *k,
			real_t *rateconv, real_t *domeigval,
			void (*out)(real_t **, int, int, int, int, int *, real_t [],
							int, real_t, real_t))
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	int j,l;
	real_t x,y,z,y0,c,d,alfa,omega,omega0,eigmax,eigeucl,euclres,maxres,
			rcmax,rceucl,maxres0,euclres0,**v,**res,auxres0,auxv,auxu,
			auxres,eucluv,maxuv;

	v=allocate_real_matrix(lj,uj,ll,ul);
	res=allocate_real_matrix(lj,uj,ll,ul);

	alfa=2.0;
	omega=4.0/(b+a);
	y0=(b+a)/(b-a);
	x=0.5*(b+a);
	y=(b-a)*(b-a)/16.0;
	z=4.0*y0*y0;
	c=a*b;
	c=sqrt(c);
	d=sqrt(a)+sqrt(b);
	d=d*d;
	if (!inap)
		for (j=lj; j<=uj; j++)
			for (l=ll; l<=ul; l++) u[j][l]=1.0;
	*k=0;
	for (j=lj; j<=uj; j++)
		for (l=ll; l<=ul; l++) res[j][l]=u[j][l];
	(*residual)(lj,uj,ll,ul,res);
	omega0=2.0/(b+a);
	maxres0=euclres0=0.0;
	for (j=lj; j<=uj; j++)
		for (l=ll; l<=ul; l++) {
			auxres0=res[j][l];
			v[j][l]=u[j][l]-omega0*auxres0;
			auxres0=fabs(auxres0);
			maxres0 = (maxres0 < auxres0) ? auxres0 : maxres0;
			euclres0 += auxres0*auxres0;
		}
	euclres0=sqrt(euclres0);
	discr[1]=euclres0;
	discr[2]=maxres0;
	(*out)(u,lj,uj,ll,ul,n,discr,*k,*rateconv,*domeigval);
	while (*k < *n) {
		(*k)++;
		/* calculate parameters alfa and omega for each iteration */
		alfa=z/(z-alfa);
		omega=1.0/(x-omega*y);
		/* iteration */
		eucluv=euclres=maxuv=maxres=0.0;
		for (j=lj; j<=uj; j++)
			for (l=ll; l<=ul; l++) res[j][l]=v[j][l];
		(*residual)(lj,uj,ll,ul,res);
		for (j=lj; j<=uj; j++)
			for (l=ll; l<=ul; l++) {
				auxv=u[j][l];
				auxu=v[j][l];
				auxres=res[j][l];
				auxv=alfa*auxu-omega*auxres+(1.0-alfa)*auxv;
				v[j][l]=auxv;
				u[j][l]=auxu;
				auxu=fabs(auxu-auxv);
				auxres=fabs(auxres);
				maxuv = (maxuv < auxu) ? auxu : maxuv;
				maxres = (maxres < auxres) ? auxres : maxres;
				eucluv += auxu*auxu;
				euclres += auxres*auxres;
			}
		eucluv=sqrt(eucluv);
		euclres=sqrt(euclres);
		discr[1]=euclres;
		discr[2]=maxres;
		maxuv=maxres/maxuv;
		eucluv=euclres/eucluv;
		eigmax=maxuv*(c-maxuv)/(0.25*d-maxuv);
		eigeucl=eucluv*(c-eucluv)/(0.25*d-eucluv);
		*domeigval=0.5*(eigmax+eigeucl);
		rceucl = -log(euclres/euclres0)/(*k);
		rcmax = -log(maxres/maxres0)/(*k);
		*rateconv=0.5*(rceucl+rcmax);
		(*out)(u,lj,uj,ll,ul,n,discr,*k,*rateconv,*domeigval);
	}
	free_real_matrix(v,lj,uj,ll);
	free_real_matrix(res,lj,uj,ll);
}
