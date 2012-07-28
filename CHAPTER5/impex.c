#include "../real.h"


int *allocate_integer_vector(int, int);
real_t *allocate_real_vector(int, int);
real_t **allocate_real_matrix(int, int, int, int);
void free_integer_vector(int *, int);
void free_real_vector(real_t *, int);
void free_real_matrix(real_t **, int, int, int);
void inivec(int, int, real_t [], real_t);
void inimat(int, int, int, int, real_t **, real_t);
void mulvec(int, int, int, real_t [], real_t [], real_t);
void mulrow(int, int, int, int, real_t **, real_t **, real_t);
void dupvec(int, int, int, real_t [], real_t []);
void duprowvec(int, int, int, real_t **, real_t []);
void dupmat(int, int, int, int, real_t **, real_t **);
real_t vecvec(int, int, int, real_t [], real_t []);
real_t matvec(int, int, int, real_t **, real_t []);
real_t matmat(int, int, int, int, real_t **, real_t **);
void elmvec(int, int, int, real_t [], real_t [], real_t);
void elmrow(int, int, int, int, real_t **, real_t **, real_t);
void dec(real_t **, int, real_t [], int []);
void sol(real_t **, int, int [], real_t []);
int impexrecomp(real_t **, real_t, real_t, real_t [], int [],
				int, int (*)(real_t, real_t [], real_t **, int),
				void (*)(real_t, real_t [], real_t [], int));
int impexlargestep(int, real_t [], real_t, real_t *,
				   real_t *, real_t *, real_t [], real_t [], real_t [],
				   real_t, real_t, real_t [], real_t [], real_t [],
				   real_t [], real_t [], real_t [], int [], int [],
				   real_t [], real_t **, real_t **, real_t,
				   void (*)(real_t, real_t [], real_t [], int),
				   int (*)(real_t, real_t [], real_t **, int));
void impexbackdiff(int, real_t [], real_t [], real_t [],
				   real_t [], real_t [], real_t [], real_t [],
				   real_t [], real_t **, real_t **);

void impex(int n, real_t t0, real_t tend, real_t y0[],
			void (*deriv)(real_t, real_t [], real_t [], int),
			int (*available)(real_t, real_t [], real_t **, int),
			real_t h0, real_t hmax, int presch, real_t eps,
			real_t weights[],
			void (*update)(real_t [], real_t [], int), int *fail,
			void (*control)(real_t *, real_t, real_t, real_t, real_t **,
							real_t [], int, real_t))
{
	int i,k,eci,*ps1,*ps2,start,two,halv;
	real_t t,t1,t2,t3,tp,h,h2,hnew,alf,lq,*y,*z,*s1,*s2,*s3,*u1,
			*u3,*w1,*w2,*w3,*ehr,**r,**rf,**a1,**a2,err[4],
			alf1,c0,c1,c2,c3,**kof,e[5],d[5],b0,b1,b2,b3,w,sl1,sn,lr;

	ps1=allocate_integer_vector(1,n);
	ps2=allocate_integer_vector(1,n);
	y=allocate_real_vector(1,n);
	z=allocate_real_vector(1,n);
	s1=allocate_real_vector(1,n);
	s2=allocate_real_vector(1,n);
	s3=allocate_real_vector(1,n);
	u1=allocate_real_vector(1,n);
	u3=allocate_real_vector(1,n);
	w1=allocate_real_vector(1,n);
	w2=allocate_real_vector(1,n);
	w3=allocate_real_vector(1,n);
	ehr=allocate_real_vector(1,n);
	r=allocate_real_matrix(1,5,1,n);
	rf=allocate_real_matrix(1,5,1,n);
	a1=allocate_real_matrix(1,n,1,n);
	a2=allocate_real_matrix(1,n,1,n);
	kof=allocate_real_matrix(2,4,2,4);

	if (presch)
		h=h0;
	else {
		if (h0 > hmax)
			h=hmax;
		else
			h=h0;
		if (h > (tend-t0)/4.0) h=(tend-t0)/4.0;
	}
	hnew=h;
	alf=0.0;
	t=tp=t0;
	inivec(1,3,err,0.0);
	inivec(1,n,ehr,0.0);
	duprowvec(1,n,1,r,y0);
	(*control)(&tp,t,h,hnew,r,err,n,tend);
	Init:
	/* initialization */
	h2=hnew;
	h=h2/2.0;
	dupvec(1,n,0,s1,y0);
	dupvec(1,n,0,s2,y0);
	dupvec(1,n,0,s3,y0);
	dupvec(1,n,0,w1,y0);
	duprowvec(1,n,1,r,y0);
	inivec(1,n,u1,0.0);
	inivec(1,n,w2,0.0);
	inimat(2,5,1,n,r,0.0);
	inimat(1,5,1,n,rf,0.0);
	t=t1=t0;
	t2=t0-2.0*h-1.0e-6;
	t3=2.0*t2+1.0;
	if (impexrecomp(a1,h,t,s1,ps1,n,available,deriv)) goto Miss;
	if (impexrecomp(a2,h2,t,w1,ps2,n,available,deriv)) goto Miss;
	start=1;
	for (eci=0; eci<=3; eci++) {
		/* one large step */
		if (impexlargestep(n,y,t,&t1,&t2,&t3,s1,s2,s3,h,h2,
							z,u1,u3,w1,w2,w3,ps1,ps2,weights,
							a1,a2,eps,deriv,available)) goto Miss;
		t += h2;
		if (eci > 0) {
			/* backward differences */
			impexbackdiff(n,u1,u3,w1,w2,w3,s1,s2,s3,r,rf);
			(*update)(weights,s2,n);
		}
	}
	eci=4;
	Mstp:
	if (hnew != h2) {
		eci=1;
		/* change of information */
		c1=hnew/h2;
		c2=c1*c1;
		c3=c2*c1;
		kof[2][2]=c1;
		kof[2][3]=(c1-c2)/2.0;
		kof[2][4]=c3/6.0-c2/2.0+c1/3.0;
		kof[3][3]=c2;
		kof[3][4]=c2-c3;
		kof[4][4]=c3;
		for (i=1; i<=n; i++) u1[i]=r[2][i]+r[3][i]/2.0+r[4][i]/3.0;
		alf1=matvec(1,n,1,rf,u1)/vecvec(1,n,0,u1,u1);
		alf=(alf+alf1)*c1;
		for (i=1; i<=n; i++) {
			e[1]=rf[1][i]-alf1*u1[i];
			e[2]=rf[2][i]-alf1*2.0*r[3][i];
			e[3]=rf[3][i]-alf1*4.0*r[4][i];
			e[4]=rf[4][i];
			d[1]=r[1][i];
			rf[1][i] = e[1] *= c2;
			for (k=2; k<=4; k++) {
				r[k][i]=d[k]=matmat(k,4,k,i,kof,r);
				rf[k][i]=e[k]=c2*matvec(k,4,k,kof,e);
			}
			s1[i]=d[1]+e[1];
			w1[i]=d[1]+4.0*e[1];
			s2[i]=s1[i]-(d[2]+e[2]/2.0);
			s3[i]=s2[i]-(d[2]+e[2])+(d[3]+e[3]/2.0);
		}
		t3=t-hnew;
		t2=t-hnew/2.0;
		t1=t;
		h2=hnew;
		h=h2/2.0;
		err[1]=0.0;
		if (halv) {
			for (i=1; i<=n; i++) ps2[i]=ps1[i];
			dupmat(1,n,1,n,a2,a1);
		}
		if (two) {
			for (i=1; i<=n; i++) ps1[i]=ps2[i];
			dupmat(1,n,1,n,a1,a2);
		} else
			if (impexrecomp(a1,hnew/2.0,t,s1,ps1,
								n,available,deriv)) goto Miss;
		if (!halv)
			if (impexrecomp(a2,hnew,t,w1,ps2,
							n,available,deriv)) goto Miss;
		/* one large step */
		if (impexlargestep(n,y,t,&t1,&t2,&t3,s1,s2,s3,h,h2,
							z,u1,u3,w1,w2,w3,ps1,ps2,weights,
							a1,a2,eps,deriv,available)) goto Miss;
		t += h2;
		eci=2;
	}
	/* one large step */
	if (impexlargestep(n,y,t,&t1,&t2,&t3,s1,s2,s3,h,h2,
							z,u1,u3,w1,w2,w3,ps1,ps2,weights,
							a1,a2,eps,deriv,available)) goto Miss;
	/* backward differences */
	impexbackdiff(n,u1,u3,w1,w2,w3,s1,s2,s3,r,rf);
	(*update)(weights,s2,n);
	/* error estimates */
	c0=c1=c2=c3=0.0;
	for (i=1; i<=n; i++) {
		w=weights[i]*weights[i];
		b0=rf[4][i]/36.0;
		c0 += b0*b0*w;
		lr=fabs(b0);
		b1=rf[1][i]+alf*r[2][i];
		c1 += b1*b1*w;
		b2=rf[3][i];
		c2 += b2*b2*w;
		sl1=fabs(rf[1][i]-rf[2][i]);
		sn = (sl1 < 1.0e-10) ? 1.0 : fabs(rf[1][i]-r[4][i]/6.0)/sl1;
		if (sn > 1.0) sn=1.0;
		if (start) {
			sn *= sn*sn*sn;
			lr *= 4.0;
		}
		ehr[i]=b3=sn*ehr[i]+lr;
		c3 += b3*b3*w;
	}
	b0=err[1];
	err[1]=b1=sqrt(c0);
	err[2]=sqrt(c1);
	err[3]=sqrt(c3)+sqrt(c2)/2.0;
	lq=eps/((b0 < b1) ? b1 : b0);
	if (b0 < b1 && lq >= 80.0) lq=10.0;
	if (eci < 4 && lq > 80.0) lq=20.0;
	halv=two=0;
	if (!presch) {
		if (lq < 1.0) {
			/* reject */
			if (start) {
				hnew=pow(lq,1.0/5.0)*h/2.0;
				goto Init;
			} else {
				for (k=1; k<=4; k++) elmrow(1,n,k,k+1,r,r,-1.0);
				for (k=1; k<=3; k++) elmrow(1,n,k,k+1,r,r,-1.0);
				for (k=1; k<=4; k++) elmrow(1,n,k,k+1,rf,rf,-1.0);
				t -= h2;
				halv=1;
				hnew=h;
				goto Mstp;
			}
		} else {
			/* step size */
			if (lq < 2.0) {
				halv=1;
				hnew=h;
			} else {
				if (lq > 80.0)
					hnew=((lq > 5120.0) ? pow(lq/5.0,1.0/5.0) : 2.0)*h2;
				if (hnew > hmax) hnew=hmax;
				if (tend > t && tend-t < hnew) hnew=tend-t;
				two=(hnew == 2.0*h2);
			}
		}
	}
	Tryck:
	if (tp <= t) (*control)(&tp,t,h,hnew,r,err,n,tend);
	if (start) start=0;
	if (hnew == h2) t += h2;
	eci++;
	if (t < tend+h2)
		goto Mstp;
	else
		goto End;
	Miss:
	*fail = presch;
	if (!(*fail)) {
		if (eci > 1) t -= h2;
		halv=two=0;
		hnew=h2/2.0;
		if (start)
			goto Init;
		else
			goto Tryck;
	}
	End:
	free_integer_vector(ps1,1);
	free_integer_vector(ps2,1);
	free_real_vector(y,1);
	free_real_vector(z,1);
	free_real_vector(s1,1);
	free_real_vector(s2,1);
	free_real_vector(s3,1);
	free_real_vector(u1,1);
	free_real_vector(u3,1);
	free_real_vector(w1,1);
	free_real_vector(w2,1);
	free_real_vector(w3,1);
	free_real_vector(ehr,1);
	free_real_matrix(r,1,5,1);
	free_real_matrix(rf,1,5,1);
	free_real_matrix(a1,1,n,1);
	free_real_matrix(a2,1,n,1);
	free_real_matrix(kof,2,4,2);
}

int impexrecomp(real_t **a, real_t h, real_t t, real_t y[], int ps[],
		int n, int (*available)(real_t, real_t [], real_t **, int),
		void (*deriv)(real_t, real_t [], real_t [], int))
{
	/* this function is internally used by IMPEX */

	int i,j;
	real_t sl,aux[4],ss,*f1,*f2;

	sl=h/2.0;
	if (!(*available)(t,y,a,n)) {
		f1=allocate_real_vector(1,n);
		f2=allocate_real_vector(1,n);
		(*deriv)(t,y,f1,n);
		for (i=1; i<=n; i++) {
			ss=1.0e-6*y[i];
			if (fabs(ss) < 1.0e-6) ss=1.0e-6;
			y[i] += ss;
			(*deriv)(t,y,f2,n);
			for (j=1; j<=n; j++) a[j][i]=(f2[j]-f1[j])/ss;
			y[i] -= ss;
		}
		free_real_vector(f1,1);
		free_real_vector(f2,1);
	}
	for (i=1; i<=n; i++) {
		mulrow(1,n,i,i,a,a,-sl);
		a[i][i] += 1.0;
	}
	aux[2]=FLT_EPSILON;
	dec(a,n,aux,ps);
	if (aux[3] < n)
		return 1;
	else
		return 0;
}

int impexlargestep(int n, real_t y[], real_t t, real_t *t1,
		real_t *t2, real_t *t3, real_t s1[], real_t s2[], real_t s3[],
		real_t h, real_t h2, real_t z[], real_t u1[], real_t u3[],
		real_t w1[], real_t w2[], real_t w3[], int ps1[], int ps2[],
		real_t weights[], real_t **a1, real_t **a2, real_t eps,
		void (*deriv)(real_t, real_t [], real_t [], int),
		int (*available)(real_t, real_t [], real_t **, int))
{
	/* this function is internally used by IMPEX */

	int impexiterate(real_t [], real_t [], real_t **, real_t, real_t,
						real_t [], int [], int, real_t,
						void (*)(real_t, real_t [], real_t [], int),
						int (*)(real_t, real_t [], real_t **, int));
	real_t a,b,c;

	a=(t+h-(*t1))/((*t1)-(*t2));
	b=(t+h-(*t2))/((*t1)-(*t3));
	c=(t+h-(*t1))/((*t2)-(*t3))*b;
	b *= a;
	a += 1.0+b;
	b=a+c-1.0;
	mulvec(1,n,0,z,s1,a);
	elmvec(1,n,0,z,s2,-b);
	elmvec(1,n,0,z,s3,c);
	if (impexiterate(z,s1,a1,h,t+h/2.0,weights,ps1,
							n,eps,deriv,available)) return 1;
	dupvec(1,n,0,y,z);
	a=(t+h2-(*t1))/((*t1)-(*t2));
	b=(t+h2-(*t2))/((*t1)-(*t3));
	c=(t+h2-(*t1))/((*t2)-(*t3))*b;
	b *= a;
	a += 1.0+b;
	b=a+c-1.0;
	mulvec(1,n,0,z,s1,a);
	elmvec(1,n,0,z,s2,-b);
	elmvec(1,n,0,z,s3,c);
	if (impexiterate(z,y,a1,h,t+3.0*h/2.0,weights,ps1,
							n,eps,deriv,available)) return 1;
	dupvec(1,n,0,u3,u1);
	dupvec(1,n,0,u1,y);
	dupvec(1,n,0,s3,s2);
	dupvec(1,n,0,s2,s1);
	dupvec(1,n,0,s1,z);
	elmvec(1,n,0,z,w1,1.0);
	elmvec(1,n,0,z,s2,-1.0);
	if (impexiterate(z,w1,a2,h2,t+h,weights,ps2,
							n,eps,deriv,available)) return 1;
	(*t3)=(*t2);
	(*t2)=(*t1);
	(*t1)=t+h2;
	dupvec(1,n,0,w3,w2);
	dupvec(1,n,0,w2,w1);
	dupvec(1,n,0,w1,z);
	return 0;
}

int impexiterate(real_t z[], real_t y[], real_t **a, real_t h, real_t t,
					real_t weights[], int ps[], int n, real_t eps,
					void (*deriv)(real_t, real_t [], real_t [], int),
					int (*available)(real_t, real_t [], real_t **, int))
{
	/* this function is internally used by IMPEXLARGESTEP (IMPEX) */

	int i,it,lit,fail;
	real_t max,max1,conv,*dz,*f1,temp;

	f1=allocate_real_vector(1,n);
	dz=allocate_real_vector(1,n);
	for (i=1; i<=n; i++) z[i]=(z[i]+y[i])/2.0;
	it=lit=1;
	conv=1.0;
	fail=0;
	while (1) {
		deriv(t,z,f1,n);
		for (i=1; i<=n; i++) f1[i]=dz[i]=z[i]-h*f1[i]/2.0-y[i];
		sol(a,n,ps,dz);
		elmvec(1,n,0,z,dz,-1.0);
		max=0.0;
		for (i=1; i<=n; i++) {
			temp=weights[i]*dz[i];
			max += temp*temp;
		}
		max=sqrt(max);
		if (max*conv < eps/10.0) break;
		it++;
		if (it != 2) {
			conv=max/max1;
			if (conv > 0.2) {
				if (lit == 0) {
					fail=1;
					break;
				}
				lit=0;
				conv=1.0;
				it=1;
				if (impexrecomp(a,h,t,z,ps,n,available,deriv)) {
					fail=1;
					break;
				}
			}
		}
		max1=max;
	}
	if (!fail)
		for (i=1; i<=n; i++) z[i]=2.0*z[i]-y[i];
	free_real_vector(f1,1);
	free_real_vector(dz,1);
	return fail;
}

void impexbackdiff(int n, real_t u1[], real_t u3[], real_t w1[],
					real_t w2[], real_t w3[], real_t s1[], real_t s2[],
					real_t s3[], real_t **r, real_t **rf)
{
	/* this function is internally used by IMPEX */

	int i,k;
	real_t b0,b1,b2,b3;

	for (i=1; i<=n; i++) {
		b1=(u1[i]+2.0*s2[i]+u3[i])/4.0;
		b2=(w1[i]+2.0*w2[i]+w3[i])/4.0;
		b3=(s3[i]+2.0*u3[i]+s2[i])/4.0;
		b2=(b2-b1)/3.0;
		b0=b1-b2;
		b2 -= (s1[i]-2.0*s2[i]+s3[i])/16.0;
		b1=2.0*b3-(b2+rf[1][i])-(b0+r[1][i])/2.0;
		b3=0.0;
		for (k=1; k<=4; k++) {
			b1 -= b3;
			b3=r[k][i];
			r[k][i]=b0;
			b0 -= b1;
		}
		r[5][i]=b0;
		for (k=1; k<=4; k++) {
			b3=rf[k][i];
			rf[k][i]=b2;
			b2 -= b3;
		}
		rf[5][i]=b2;
	}
}
