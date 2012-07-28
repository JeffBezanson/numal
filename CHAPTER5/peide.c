#include "../real.h"


int *allocate_integer_vector(int, int);
real_t *allocate_real_vector(int, int);
real_t **allocate_real_matrix(int, int, int, int);
void free_integer_vector(int *, int);
void free_real_vector(real_t *, int);
void free_real_matrix(real_t **, int, int, int);
int peidefunct(int nrow, int ncol, real_t par[], real_t res[],
			   int n, int m, int nobs, int *nbp, int first, int *sec,
			   int *max, int *nis, real_t eps1, int weight, int bp[],
			   real_t save[], real_t ymax[], real_t y[], real_t **yp,
			   real_t **fy, real_t **fp, int cobs[], real_t tobs[],
			   real_t obs[], real_t in[], real_t aux[], int clean,
			   int (*deriv)(int,int,real_t [],real_t [],real_t,real_t []),
			   int (*jacdfdy)(int,int,real_t [],real_t [],real_t,real_t **),
			   int (*jacdfdp)(int,int,real_t [],real_t [],real_t,real_t **),
			   void (*callystart)(int,int,real_t [],real_t [],real_t[]),
			   void (*monitor)(int,int,int,real_t [],real_t [],int,int));
void inivec(int, int, real_t [], real_t);
void inimat(int, int, int, int, real_t **, real_t);
void mulvec(int, int, int, real_t [], real_t [], real_t);
void mulrow(int, int, int, int, real_t **, real_t **, real_t);
void dupvec(int, int, int, real_t [], real_t []);
void dupmat(int, int, int, int, real_t **, real_t **);
real_t vecvec(int, int, int, real_t [], real_t []);
real_t matvec(int, int, int, real_t **, real_t []);
void elmvec(int, int, int, real_t [], real_t [], real_t);
void sol(real_t **, int, int [], real_t []);
void dec(real_t **, int, real_t [], int []);
void mulcol(int, int, int, int, real_t **, real_t **, real_t);
real_t tamvec(int, int, int, real_t **, real_t []);
real_t mattam(int, int, int, int, real_t **, real_t **);
int qrisngvaldec(real_t **,int,int,real_t [],real_t **,real_t []);

void peide(int n, int m, int nobs, int *nbp, real_t par[],
		real_t res[], int bp[], real_t **jtjinv,
		real_t in[], real_t out[],
		int (*deriv)(int,int,real_t [],real_t [],real_t,real_t []),
		int (*jacdfdy)(int,int,real_t [],real_t [],real_t,real_t **),
		int (*jacdfdp)(int,int,real_t [],real_t [],real_t,real_t **),
		void (*callystart)(int,int,real_t [],real_t [],real_t[]),
		void (*data)(int,real_t [],real_t [],int[]),
		void (*monitor)(int,int,int,real_t [],real_t [],int,int))
{
	int i,j,weight,ncol,nrow,away,max,nfe,nis,*cobs,
			first,sec,clean,nbpold,maxfe,fe,it,err,emergency;
	real_t eps1,res1,in3,in4,fac3,fac4,aux[4],*obs,*save,*tobs,
			**yp,*ymax,*y,**fy,**fp,w,**aid,temp,
			vv,ww,w2,mu,res2,fpar,fparpres,lambda,lambdamin,p,pw,
			reltolres,abstolres,em[8],*val,*b,*bb,*parpres,**jaco;
	static real_t save1[35]={1.0, 1.0, 9.0, 4.0, 0.0, 2.0/3.0, 1.0,
			1.0/3.0, 36.0, 20.25, 1.0, 6.0/11.0, 1.0, 6.0/11.0,
			1.0/11.0, 84.028, 53.778, 0.25, 0.48, 1.0, 0.7, 0.2,
			0.02, 156.25, 108.51, 0.027778, 120.0/274.0, 1.0,
			225.0/274.0, 85.0/274.0, 15.0/274.0, 1.0/274.0, 0.0,
			187.69, 0.0047361};

	nbpold=(*nbp);
	cobs=allocate_integer_vector(1,nobs);
	obs=allocate_real_vector(1,nobs);
	save=allocate_real_vector(-38,6*n);
	tobs=allocate_real_vector(0,nobs);
	ymax=allocate_real_vector(1,n);
	y=allocate_real_vector(1,6*n*(nbpold+m+1));
	yp=allocate_real_matrix(1,nbpold+nobs,1,nbpold+m);
	fy=allocate_real_matrix(1,n,1,n);
	fp=allocate_real_matrix(1,n,1,m+nbpold);
	aid=allocate_real_matrix(1,m+nbpold,1,m+nbpold);

	for (i=0; i<=34; i++) save[-38+i]=save1[i];
	(*data)(nobs,tobs,obs,cobs);
	weight=1;
	first=sec=0;
	clean=(*nbp > 0);
	aux[2]=FLT_EPSILON;
	eps1=1.0e10;
	out[1]=0.0;
	bp[0]=max=0;
	/* smooth integration without break-points */
	if (!peidefunct(nobs,m,par,res,
			n,m,nobs,nbp,first,&sec,&max,&nis,eps1,weight,bp,
			save,ymax,y,yp,fy,fp,cobs,tobs,obs,in,aux,clean,deriv,
			jacdfdy,jacdfdp,callystart,monitor)) goto Escape;
	res1=sqrt(vecvec(1,nobs,0,res,res));
	nfe=1;
	if (in[5] == 1.0) {
		out[1]=1.0;
		goto Escape;
	}
	if (clean) {
		first=1;
		clean=0;
		fac3=sqrt(sqrt(in[3]/res1));
		fac4=sqrt(sqrt(in[4]/res1));
		eps1=res1*fac4;
		if (!peidefunct(nobs,m,par,res,
				n,m,nobs,nbp,first,&sec,&max,&nis,eps1,weight,bp,
				save,ymax,y,yp,fy,fp,cobs,tobs,obs,in,aux,clean,deriv,
				jacdfdy,jacdfdp,callystart,monitor)) goto Escape;
		first=0;
	} else
		nfe=0;
	ncol=m+(*nbp);
	nrow=nobs+(*nbp);
	sec=1;
	in3=in[3];
	in4=in[4];
	in[3]=res1;
	weight=away=0;
	out[4]=out[5]=w=0.0;
	temp=sqrt(weight)+1.0;
	weight=temp*temp;
	while (weight != 16 && *nbp > 0) {
		if (away == 0 && w != 0.0) {
			/* if no break-points were omitted then one function
				function evaluation is saved */
			w=weight/w;
			for (i=nobs+1; i<=nrow; i++) {
				for (j=1; j<=ncol; j++) yp[i][j] *= w;
				res[i] *= w;
			}
			sec=1;
			nfe--;
		}
		in[3] *= fac3*weight;
		in[4]=eps1;
		(*monitor)(2,ncol,nrow,par,res,weight,nis);
		/* marquardt's method */
		val=allocate_real_vector(1,ncol);
		b=allocate_real_vector(1,ncol);
		bb=allocate_real_vector(1,ncol);
		parpres=allocate_real_vector(1,ncol);
		jaco=allocate_real_matrix(1,nrow,1,ncol);
		vv=10.0;
		w2=0.5;
		mu=0.01;
		ww = (in[6] < 1.0e-7) ? 1.0e-8 : 1.0e-1*in[6];
		em[0]=em[2]=em[6]=in[0];
		em[4]=10*ncol;
		reltolres=in[3];
		abstolres=in[4]*in[4];
		maxfe=in[5];
		err=0;
		fe=it=1;
		p=fpar=res2=0.0;
		pw = -log(ww*in[0])/2.30;
		if (!peidefunct(nrow,ncol,par,res,
					n,m,nobs,nbp,first,&sec,&max,&nis,eps1,
					weight,bp,save,ymax,y,yp,fy,fp,cobs,tobs,obs,
					in,aux,clean,deriv,jacdfdy,jacdfdp,
					callystart,monitor))
			err=3;
		else {
			fpar=vecvec(1,nrow,0,res,res);
			out[3]=sqrt(fpar);
			emergency=0;
			it=1;
			do {
				dupmat(1,nrow,1,ncol,jaco,yp);
				i=qrisngvaldec(jaco,nrow,ncol,val,aid,em);
				if (it == 1)
					lambda=in[6]*vecvec(1,ncol,0,val,val);
				else
					if (p == 0.0) lambda *= w2;
				for (i=1; i<=ncol; i++)
					b[i]=val[i]*tamvec(1,nrow,i,jaco,res);
				while (1) {
					for (i=1; i<=ncol; i++)
						bb[i]=b[i]/(val[i]*val[i]+lambda);
					for (i=1; i<=ncol; i++)
						parpres[i]=par[i]-matvec(1,ncol,i,aid,bb);
					fe++;
					if (fe >= maxfe)
						err=1;
					else
						if (!peidefunct(nrow,ncol,parpres,res,
								n,m,nobs,nbp,first,&sec,&max,&nis,
								eps1,weight,bp,save,ymax,y,yp,fy,fp,
								cobs,tobs,obs,in,aux,clean,deriv,
								jacdfdy,jacdfdp,callystart,monitor))
							err=2;
					if (err != 0) {
						emergency=1;
						break;
					}
					fparpres=vecvec(1,nrow,0,res,res);
					res2=fpar-fparpres;
					if (res2 < mu*vecvec(1,ncol,0,b,bb)) {
						p += 1.0;
						lambda *= vv;
						if (p == 1.0) {
							lambdamin=ww*vecvec(1,ncol,0,val,val);
							if (lambda < lambdamin) lambda=lambdamin;
						}
						if (p >= pw) {
							err=4;
							emergency=1;
							break;
						}
					} else {
						dupvec(1,ncol,0,par,parpres);
						fpar=fparpres;
						break;
					}
				}
				if (emergency) break;
				it++;
			} while (fpar>abstolres &&
							res2>reltolres*fpar+abstolres);
			for (i=1; i<=ncol; i++)
				mulcol(1,ncol,i,i,jaco,aid,1.0/(val[i]+in[0]));
			for (i=1; i<=ncol; i++)
				for (j=1; j<=i; j++)
					aid[i][j]=aid[j][i]=mattam(1,ncol,i,j,jaco,jaco);
			lambda=lambdamin=val[1];
			for (i=2; i<=ncol; i++)
				if (val[i] > lambda)
					lambda=val[i];
				else
					if (val[i] < lambdamin) lambdamin=val[i];
			temp=lambda/(lambdamin+in[0]);
			out[7]=temp*temp;
			out[2]=sqrt(fpar);
			out[6]=sqrt(res2+fpar)-out[2];
		}
		out[4]=fe;
		out[5]=it-1;
		out[1]=err;
		free_real_vector(val,1);
		free_real_vector(b,1);
		free_real_vector(bb,1);
		free_real_vector(parpres,1);
		free_real_matrix(jaco,1,nrow,1);
		if (out[1] > 0.0) goto Escape;
		/* the relative starting value of lambda is adjusted
			to the last value of lambda used */
		away=out[4]-out[5]-1.0;
		in[6] *= pow(5.0,away)*pow(2.0,away-out[5]);
		nfe += out[4];
		w=weight;
		temp=sqrt(weight)+1.0;
		eps1=temp*temp*in[4]*fac4;
		away=0;
		/* omit useless break-points */
		for (j=1; j<=(*nbp); j++)
			if (fabs(obs[bp[j]]+res[bp[j]]-par[j+m]) < eps1) {
				(*nbp)--;
				for (i=j; i<=(*nbp); i++) bp[i]=bp[i+1];
				dupvec(j+m,(*nbp)+m,1,par,par);
				j--;
				away++;
				bp[*nbp+1]=0;
			}
		ncol -= away;
		nrow -= away;
		temp=sqrt(weight)+1.0;
		weight=temp*temp;
	}
	in[3]=in3;
	in[4]=in4;
	*nbp=0;
	weight=1;
	(*monitor)(2,m,nobs,par,res,weight,nis);
	/* marquardt's method */
	val=allocate_real_vector(1,m);
	b=allocate_real_vector(1,m);
	bb=allocate_real_vector(1,m);
	parpres=allocate_real_vector(1,m);
	jaco=allocate_real_matrix(1,nobs,1,m);
	vv=10.0;
	w2=0.5;
	mu=0.01;
	ww = (in[6] < 1.0e-7) ? 1.0e-8 : 1.0e-1*in[6];
	em[0]=em[2]=em[6]=in[0];
	em[4]=10*m;
	reltolres=in[3];
	abstolres=in[4]*in[4];
	maxfe=in[5];
	err=0;
	fe=it=1;
	p=fpar=res2=0.0;
	pw = -log(ww*in[0])/2.30;
	if (!peidefunct(nobs,m,par,res,
				n,m,nobs,nbp,first,&sec,&max,&nis,eps1,weight,bp,
				save,ymax,y,yp,fy,fp,cobs,tobs,obs,in,aux,clean,
				deriv,jacdfdy,jacdfdp,callystart,monitor))
		err=3;
	else {
		fpar=vecvec(1,nobs,0,res,res);
		out[3]=sqrt(fpar);
		emergency=0;
		it=1;
		do {
			dupmat(1,nobs,1,m,jaco,yp);
			i=qrisngvaldec(jaco,nobs,m,val,jtjinv,em);
			if (it == 1)
				lambda=in[6]*vecvec(1,m,0,val,val);
			else
				if (p == 0.0) lambda *= w2;
			for (i=1; i<=m; i++)
				b[i]=val[i]*tamvec(1,nobs,i,jaco,res);
			while (1) {
				for (i=1; i<=m; i++)
					bb[i]=b[i]/(val[i]*val[i]+lambda);
				for (i=1; i<=m; i++)
					parpres[i]=par[i]-matvec(1,m,i,jtjinv,bb);
				fe++;
				if (fe >= maxfe)
					err=1;
				else
					if (!peidefunct(nobs,m,parpres,res,
							n,m,nobs,nbp,first,&sec,&max,&nis,eps1,
							weight,bp,save,ymax,y,yp,fy,fp,cobs,tobs,
							obs,in,aux,clean,deriv,jacdfdy,jacdfdp,
							callystart,monitor))
						err=2;
				if (err != 0) {
					emergency=1;
					break;
				}
				fparpres=vecvec(1,nobs,0,res,res);
				res2=fpar-fparpres;
				if (res2 < mu*vecvec(1,m,0,b,bb)) {
					p += 1.0;
					lambda *= vv;
					if (p == 1.0) {
						lambdamin=ww*vecvec(1,m,0,val,val);
						if (lambda < lambdamin) lambda=lambdamin;
					}
					if (p >= pw) {
						err=4;
						emergency=1;
						break;
					}
				} else {
					dupvec(1,m,0,par,parpres);
					fpar=fparpres;
					break;
				}
			}
			if (emergency) break;
			it++;
		} while (fpar>abstolres && res2>reltolres*fpar+abstolres);
		for (i=1; i<=m; i++)
			mulcol(1,m,i,i,jaco,jtjinv,1.0/(val[i]+in[0]));
		for (i=1; i<=m; i++)
			for (j=1; j<=i; j++)
				jtjinv[i][j]=jtjinv[j][i]=mattam(1,m,i,j,jaco,jaco);
		lambda=lambdamin=val[1];
		for (i=2; i<=m; i++)
			if (val[i] > lambda)
				lambda=val[i];
			else
				if (val[i] < lambdamin) lambdamin=val[i];
		temp=lambda/(lambdamin+in[0]);
		out[7]=temp*temp;
		out[2]=sqrt(fpar);
		out[6]=sqrt(res2+fpar)-out[2];
	}
	out[4]=fe;
	out[5]=it-1;
	out[1]=err;
	free_real_vector(val,1);
	free_real_vector(b,1);
	free_real_vector(bb,1);
	free_real_vector(parpres,1);
	free_real_matrix(jaco,1,nobs,1);
	nfe += out[4];

	Escape:
	if (out[1] == 3.0)
		out[1]=2.0;
	else
		if (out[1] == 4.0) out[1]=6.0;
	if (save[-3] != 0.0) out[1]=save[-3];
	out[3]=res1;
	out[4]=nfe;
	out[5]=max;
	free_integer_vector(cobs,1);
	free_real_vector(obs,1);
	free_real_vector(save,-38);
	free_real_vector(tobs,0);
	free_real_vector(ymax,1);
	free_real_vector(y,1);
	free_real_matrix(yp,1,nbpold+nobs,1);
	free_real_matrix(fy,1,n,1);
	free_real_matrix(fp,1,n,1);
	free_real_matrix(aid,1,m+nbpold,1);
}

int peidefunct(int nrow, int ncol, real_t par[], real_t res[],
		int n, int m, int nobs, int *nbp, int first, int *sec,
		int *max, int *nis, real_t eps1, int weight, int bp[],
		real_t save[], real_t ymax[], real_t y[], real_t **yp,
		real_t **fy, real_t **fp, int cobs[], real_t tobs[],
		real_t obs[], real_t in[], real_t aux[], int clean,
		int (*deriv)(int,int,real_t [],real_t [],real_t,real_t []),
		int (*jacdfdy)(int,int,real_t [],real_t [],real_t,real_t **),
		int (*jacdfdp)(int,int,real_t [],real_t [],real_t,real_t **),
		void (*callystart)(int,int,real_t [],real_t [],real_t[]),
		void (*monitor)(int,int,int,real_t [],real_t [],int,int))
{
	/* this function is internally used by PEIDE */

	void peidereset(int, int, real_t, real_t, real_t, real_t, real_t [],
				real_t [], real_t *, real_t *, real_t *, int *);
	void peideorder(int, int, real_t, real_t [], real_t [],
			real_t *, real_t *, real_t *, real_t *, real_t *, int *);
	void peidestep(int, int, int, real_t, real_t, real_t, real_t,
			real_t [], real_t [], real_t [], real_t [], int *, real_t *);
	real_t peideinterpol(int, int, int, real_t, real_t []);
	int l,k,knew,fails,same,kpold,n6,nnpar,j5n,cobsii,*p,evaluate,
			evaluated,decompose,conv,extra,npar,i,j,jj,ii;
	real_t xold,hold,a0,tolup,tol,toldwn,tolconv,h,ch,chnew,error,
			dfi,tobsdif,a[6],*delta,*lastdelta,*df,*y0,**jacob,xend,
			hmax,hmin,eps,s,aa,x,t,c;

	p=allocate_integer_vector(1,n);
	delta=allocate_real_vector(1,n);
	lastdelta=allocate_real_vector(1,n);
	df=allocate_real_vector(1,n);
	y0=allocate_real_vector(1,n);
	jacob=allocate_real_matrix(1,n,1,n);

	if (*sec) {
		*sec=0;
		goto Finish;
	}
	xend=tobs[nobs];
	eps=in[2];
	npar=m;
	extra=(*nis)=0;
	ii=1;
	jj = (*nbp == 0) ? 0 : 1;
	n6=n*6;
	inivec(-3,-1,save,0.0);
	inivec(n6+1,(6+m)*n,y,0.0);
	inimat(1,nobs+(*nbp),1,m+(*nbp),yp,0.0);
	t=tobs[1];
	x=tobs[0];
	(*callystart)(n,m,par,y,ymax);
	hmax=tobs[1]-tobs[0];
	hmin=hmax*in[1];
	/* evaluate jacobian */
	evaluate=0;
	decompose=evaluated=1;
	if (!(*jacdfdy)(n,m,par,y,x,fy)) {
		save[-3]=4.0;
		goto Finish;
	}
	nnpar=n*npar;

	Newstart:
	k=1;
	kpold=0;
	same=2;
	peideorder(n,k,eps,a,save,&tol,&tolup,&toldwn,&tolconv,
					&a0,&decompose);
	if (!(*deriv)(n,m,par,y,x,df)) {
		save[-3]=3.0;
		goto Finish;
	}
	s=FLT_MIN;
	for (i=1; i<=n; i++) {
		aa=matvec(1,n,i,fy,df)/ymax[i];
		s += aa*aa;
	}
	h=sqrt(2.0*eps/sqrt(s));
	if (h > hmax)
		h=hmax;
	else
		if (h < hmin) h=hmin;
	xold=x;
	hold=h;
	ch=1.0;
	for (i=1; i<=n; i++) {
		save[i]=y[i];
		save[n+i]=y[n+i]=df[i]*h;
	}
	fails=0;
	while (x < xend) {
		if (x+h <= xend)
			x += h;
		else {
			h=xend-x;
			x=xend;
			ch=h/hold;
			c=1.0;
			for (j=n; j<=k*n; j += n) {
				c *= ch;
				for (i=j+1; i<=j+n; i++) y[i] *= c;
			}
			same = (same < 3) ? 3 : same+1;
		}
		/* prediction */
		for (l=1; l<=n; l++) {
			for (i=l; i<=(k-1)*n+l; i += n)
				for (j=(k-1)*n+l; j>=i; j -= n) y[j] += y[j+n];
			delta[l]=0.0;
		}
		evaluated=0;
		/* correction and estimation local error */
		for (l=1; l<=3; l++) {
			if (!(*deriv)(n,m,par,y,x,df)) {
				save[-3]=3;
				goto Finish;
			}
			for (i=1; i<=n; i++) df[i]=df[i]*h-y[n+i];
			if (evaluate) {
				/* evaluate jacobian */
				evaluate=0;
				decompose=evaluated=1;
				if (!(*jacdfdy)(n,m,par,y,x,fy)) {
					save[-3]=4.0;
					goto Finish;
				}
			}
			if (decompose) {
				/* decompose jacobian */
				decompose=0;
				c = -a0*h;
				for (j=1; j<=n; j++) {
					for (i=1; i<=n; i++) jacob[i][j]=fy[i][j]*c;
					jacob[j][j] += 1.0;
				}
				dec(jacob,n,aux,p);
			}
			sol(jacob,n,p,df);
			conv=1;
			for (i=1; i<=n; i++) {
				dfi=df[i];
				y[i] += a0*dfi;
				y[n+i] += dfi;
				delta[i] += dfi;
				conv=(conv && (fabs(dfi) < tolconv*ymax[i]));
			}
			if (conv) {
				s=FLT_MIN;
				for (i=1; i<=n; i++) {
					aa=delta[i]/ymax[i];
					s += aa*aa;
				}
				error=s;
				break;
			}
		}
		/* acceptance or rejection */
		if (!conv) {
			if (!evaluated)
				evaluate=1;
			else {
				ch /= 4.0;
				if (h < 4.0*hmin) {
					save[-1] += 10.0;
					hmin /= 10.0;
					if (save[-1] > 40.0) goto Finish;
				}
			}
			peidereset(n,k,hmin,hmax,hold,xold,y,save,&ch,&x,
							&h,&decompose);
		} else if (error > tol) {
			fails++;
			if (h > 1.1*hmin) {
				if (fails > 2) {
					peidereset(n,k,hmin,hmax,hold,xold,y,save,&ch,&x,
								&h,&decompose);
					goto Newstart;
				} else {
					/* calculate step and order */
					peidestep(n,k,fails,tolup,toldwn,tol,error,delta,
								lastdelta,y,ymax,&knew,&chnew);
					if (knew != k) {
						k=knew;
						peideorder(n,k,eps,a,save,&tol,&tolup,
									&toldwn,&tolconv,&a0,&decompose);
					}
					ch *= chnew;
					peidereset(n,k,hmin,hmax,hold,xold,y,save,&ch,&x,
								&h,&decompose);
				}
			} else {
				if (k == 1) {
					/* violate eps criterion */
					save[-2] += 1.0;
					same=4;
					goto Errortestok;
				}
				k=1;
				peidereset(n,k,hmin,hmax,hold,xold,y,save,&ch,&x,
							&h,&decompose);
				peideorder(n,k,eps,a,save,&tol,&tolup,
							&toldwn,&tolconv,&a0,&decompose);
				same=2;
			}
		} else {
			Errortestok:
			fails=0;
			for (i=1; i<=n; i++) {
				c=delta[i];
				for (l=2; l<=k; l++) y[l*n+i] += a[l]*c;
				if (fabs(y[i]) > ymax[i]) ymax[i]=fabs(y[i]);
			}
			same--;
			if (same == 1)
				dupvec(1,n,0,lastdelta,delta);
			else if (same == 0) {
				/* calculate step and order */
				peidestep(n,k,fails,tolup,toldwn,tol,error,delta,
							lastdelta,y,ymax,&knew,&chnew);
				if (chnew > 1.1) {
					if (k != knew) {
						if (knew > k)
							mulvec(knew*n+1,knew*n+n,-knew*n,y,delta,
									a[k]/knew);
						k=knew;
						peideorder(n,k,eps,a,save,&tol,&tolup,
									&toldwn,&tolconv,&a0,&decompose);
					}
					same=k+1;
					if (chnew*h > hmax) chnew=hmax/h;
					h *= chnew;
					c=1.0;
					for (j=n; j<=k*n; j += n) {
						c *= chnew;
						mulvec(j+1,j+n,0,y,y,c);
					}
					decompose=1;
				} else
					same=10;
			}
			(*nis)++;
			/* start of an integration step of yp */
			if (clean) {
				hold=h;
				xold=x;
				kpold=k;
				ch=1.0;
				dupvec(1,k*n+n,0,save,y);
			} else {
				if (h != hold) {
					ch=h/hold;
					c=1.0;
					for (j=n6+nnpar; j<=kpold*nnpar+n6; j += nnpar) {
						c *= ch;
						for (i=j+1; i<=j+nnpar; i++) y[i] *= c;
					}
					hold=h;
				}
				if (k > kpold)
					inivec(n6+k*nnpar+1,n6+k*nnpar+nnpar,y,0.0);
				xold=x;
				kpold=k;
				ch=1.0;
				dupvec(1,k*n+n,0,save,y);
				/* evaluate jacobian */
				evaluate=0;
				decompose=evaluated=1;
				if (!(*jacdfdy)(n,m,par,y,x,fy)) {
					save[-3]=4.0;
					goto Finish;
				}
				/* decompose jacobian */
				decompose=0;
				c = -a0*h;
				for (j=1; j<=n; j++) {
					for (i=1; i<=n; i++) jacob[i][j]=fy[i][j]*c;
					jacob[j][j] += 1.0;
				}
				dec(jacob,n,aux,p);
				if (!(*jacdfdp)(n,m,par,y,x,fp)) {
					save[-3]=5.0;
					goto Finish;
				}
				if (npar > m) inimat(1,n,m+1,npar,fp,0.0);
				/* prediction */
				for (l=0; l<=k-1; l++)
					for (j=k-1; j>=l; j--)
						elmvec(j*nnpar+n6+1,j*nnpar+n6+nnpar,nnpar,
									y,y,1.0);
				/* correction */
				for (j=1; j<=npar; j++) {
					j5n=(j+5)*n;
					dupvec(1,n,j5n,y0,y);
					for (i=1; i<=n; i++)
						df[i]=h*(fp[i][j]+matvec(1,n,i,fy,y0))-
									y[nnpar+j5n+i];
					sol(jacob,n,p,df);
					for (l=0; l<=k; l++) {
						i=l*nnpar+j5n;
						elmvec(i+1,i+n,-i,y,df,a[l]);
					}
				}
			}
			while (x >= t) {
				/* calculate a row of the jacobian matrix and an
					element of the residual vector */
				tobsdif=(tobs[ii]-x)/h;
				cobsii=cobs[ii];
				res[ii]=peideinterpol(cobsii,n,k,tobsdif,y)-obs[ii];
				if (!clean) {
					for (i=1; i<=npar; i++)
						yp[ii][i]=peideinterpol(cobsii+(i+5)*n,nnpar,k,
														tobsdif,y);
					/* introducing break-points */
					if (bp[jj] != ii) {
					} else if (first && fabs(res[ii]) < eps1) {
						(*nbp)--;
						for (i=jj; i<=(*nbp); i++) bp[i]=bp[i+1];
						bp[*nbp+1]=0;
					} else {
						extra++;
						if (first) par[m+jj]=obs[ii];
						/* introducing a jacobian row and a residual
							vector element for continuity requirements */
						yp[nobs+jj][m+jj] = -weight;
						mulrow(1,npar,nobs+jj,ii,yp,yp,weight);
						res[nobs+jj]=weight*(res[ii]+obs[ii]-par[m+jj]);
					}
				}
				if (ii == nobs)
					goto Finish;
				else {
					t=tobs[ii+1];
					if (bp[jj] == ii && jj < *nbp) jj++;
					hmax=t-tobs[ii];
					hmin=hmax*in[1];
					ii++;
				}
			}
			/* break-points introduce new initial values for y & yp */
			if (extra > 0) {
				for (i=1; i<=n; i++) {
					y[i]=peideinterpol(i,n,k,tobsdif,y);
					for (j=1; j<=npar; j++)
						y[i+(j+5)*n]=peideinterpol(i+(j+5)*n,nnpar,
															k,tobsdif,y);
				}
				for (l=1; l<=extra; l++) {
					cobsii=cobs[bp[npar-m+l]];
					y[cobsii]=par[npar+l];
					for (i=1; i<=npar+extra; i++) y[cobsii+(5+i)*n]=0.0;
					inivec(1+nnpar+(l+5)*n,nnpar+(l+6)*n,y,0.0);
					y[cobsii+(5+npar+l)*n]=1.0;
				}
				npar += extra;
				extra=0;
				x=tobs[ii-1];
				/* evaluate jacobian */
				evaluate=0;
				decompose=evaluated=1;
				if (!(*jacdfdy)(n,m,par,y,x,fy)) {
					save[-3]=4.0;
					goto Finish;
				}
				nnpar=n*npar;
				goto Newstart;
			}
		}
	}
	Finish:
	if (save[-2] > *max) *max=save[-2];
	if (!first) (*monitor)(1,ncol,nrow,par,res,weight,*nis);
	free_integer_vector(p,1);
	free_real_vector(delta,1);
	free_real_vector(lastdelta,1);
	free_real_vector(df,1);
	free_real_vector(y0,1);
	free_real_matrix(jacob,1,n,1);
	return (save[-1] <= 40.0 && save[-3] == 0.0);
}

void peidereset(int n, int k, real_t hmin, real_t hmax, real_t hold,
		real_t xold, real_t y[], real_t save[], real_t *ch, real_t *x,
		real_t *h, int *decompose)
{
	/* this function is internally used by PEIDEFUNCT of PEIDE */

	int i,j;
	real_t c;

	if (*ch < hmin/hold)
		*ch = hmin/hold;
	else
		if (*ch > hmax/hold) *ch = hmax/hold;
	*x = xold;
	*h = hold*(*ch);
	c=1.0;
	for (j=0; j<=k*n; j += n) {
		for (i=1; i<=n; i++) y[j+i]=save[j+i]*c;
		c *= (*ch);
	}
	*decompose = 1;
}

void peideorder(int n, int k, real_t eps, real_t a[], real_t save[],
		real_t *tol, real_t *tolup, real_t *toldwn, real_t *tolconv,
		real_t *a0, int *decompose)
{
	/* this function is internally used by PEIDEFUNCT of PEIDE */

	int i,j;
	real_t c;

	c=eps*eps;
	j=((k-1)*(k+8))/2-38;
	for (i=0; i<=k; i++) a[i]=save[i+j];
	j += k+1;
	*tolup = c*save[j];
	*tol = c*save[j+1];
	*toldwn = c*save[j+2];
	*tolconv = eps/(2*n*(k+2));
	*a0 = a[0];
	*decompose = 1;
}

void peidestep(int n, int k, int fails, real_t tolup, real_t toldwn,
		real_t tol, real_t error, real_t delta[], real_t lastdelta[],
		real_t y[], real_t ymax[], int *knew, real_t *chnew)
{
	/* this function is internally used by PEIDEFUNCT of PEIDE */

	int i;
	real_t a1,a2,a3,aa,s;

	if (k <= 1)
		a1=0.0;
	else {
		s=FLT_MIN;
		for (i=1; i<=n; i++) {
			aa=y[k*n+i]/ymax[i];
			s += aa*aa;
		}
		a1=0.75*pow(toldwn/s,0.5/k);
	}
	a2=0.80*pow(tol/error,0.5/(k+1));
	if (k >= 5 || fails != 0)
		a3=0.0;
	else {
		s=FLT_MIN;
		for (i=1; i<=n; i++) {
			aa=(delta[i]-lastdelta[i])/ymax[i];
			s += aa*aa;
		}
		a3=0.70*pow(tolup/s,0.5/(k+2));
	}
	if (a1 > a2 && a1 > a3) {
		*knew = k-1;
		*chnew = a1;
	} else if (a2 > a3) {
		*knew = k;
		*chnew = a2;
	} else {
		*knew = k+1;
		*chnew = a3;
	}
}

real_t peideinterpol(int startindex, int jump, int k, real_t tobsdif,
							real_t y[])
{
	/* this function is internally used by PEIDEFUNCT of PEIDE */

	int i;
	real_t s,r;

	s=y[startindex];
	r=tobsdif;
	for (i=1; i<=k; i++) {
		startindex += jump;
		s += y[startindex]*r;
		r *= tobsdif;
	}
	return s;
}

