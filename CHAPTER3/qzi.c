#include "../real.h"


void qzi(int n, real_t **a, real_t **b, real_t **x, real_t alfr[],
			real_t alfi[], real_t beta[], int iter[], real_t em[])
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void hshdecmul(int, real_t **, real_t **, real_t);
	void hestgl3(int, real_t **, real_t **, real_t **);
	void hsh2row3(int, int, int, int, int, real_t, real_t,
						real_t **, real_t **, real_t **);
	void hsh3row3(int, int, int, int, real_t, real_t, real_t,
						real_t **, real_t **, real_t **);
	void hsh2col(int, int, int, int, real_t, real_t, real_t **, real_t **);
	void hsh3col(int, int, int, int, real_t, real_t, real_t,
					real_t **, real_t **);
	void chsh2(real_t, real_t, real_t, real_t, real_t *, real_t *, real_t *);
	void comdiv(real_t, real_t, real_t, real_t, real_t *, real_t *);
	int i,q,m,m1,q1,j,k,k1,k2,k3,km1,stationary,goon,l,mr,mi,l1,out;
	real_t dwarf,eps,epsa,epsb,
			anorm,bnorm,ani,bni,constt,a10,a20,a30,b11,b22,b33,b44,a11,
			a12,a21,a22,a33,a34,a43,a44,b12,b34,old1,old2,
			an,bn,e,c,d,er,ei,a11r,a11i,a12r,a12i,a21r,a21i,a22r,a22i,
			cz,szr,szi,cq,sqr,sqi,ssr,ssi,tr,ti,bdr,bdi,r,
			betm,alfm,sl,sk,tkk,tkl,tlk,tll,almi,almr,slr,sli,skr,ski,
			dr,di,tkkr,tkki,tklr,tkli,tlkr,tlki,tllr,tlli,s;

	dwarf=em[0];
	eps=em[1];
	hshdecmul(n,a,b,dwarf);
	hestgl3(n,a,b,x);
	anorm=bnorm=0.0;
	for (i=1; i<=n; i++) {
		bni=0.0;
		iter[i]=0;
		ani = (i > 1) ? fabs(a[i][i-1]) : 0.0;
		for (j=i; j<=n; j++) {
			ani += fabs(a[i][j]);
			bni += fabs(b[i][j]);
		}
		if (ani > anorm) anorm=ani;
		if (bni > bnorm) bnorm=bni;
	}
	if (anorm == 0.0) anorm=eps;
	if (bnorm == 0.0) bnorm=eps;
	epsa=eps*anorm;
	epsb=eps*bnorm;
	m=n;
	out=0;
	do {
		i=q=m;
		while ((i > 1) ? fabs(a[i][i-1]) > epsa : 0) {
			q=i-1;
			i--;
		}
		if (q > 1) a[q][q-1]=0.0;
		goon=1;
		while (goon) {
			if (q >= m-1) {
				m=q-1;
				goon=0;
			} else {
				if (fabs(b[q][q]) <= epsb) {
					b[q][q]=0.0;
					q1=q+1;
					hsh2col(q,q,n,q,a[q][q],a[q1][q],a,b);
					a[q1][q]=0.0;
					q=q1;
				} else {
					goon=0;
					m1=m-1;
					q1=q+1;
					constt=0.75;
					(iter[m])++;
					stationary = (iter[m] == 1) ? 1 :
						(fabs(a[m][m-1]) >= constt*old1 &&
						fabs(a[m-1][m-2]) >= constt*old2);
					if (iter[m] > 30 && stationary) {
						for (i=1; i<=m; i++) iter[i] = -1;
						out=1;
						break;
					}
					if (iter[m] == 10 && stationary) {
						a10=0.0;
						a20=1.0;
						a30=1.1605;
					} else {
						b11=b[q][q];
						b22 = (fabs(b[q1][q1]) < epsb) ? epsb : b[q1][q1];
						b33 = (fabs(b[m1][m1]) < epsb) ? epsb : b[m1][m1];
						b44 = (fabs(b[m][m]) < epsb) ? epsb : b[m][m];
						a11=a[q][q]/b11;
						a12=a[q][q1]/b22;
						a21=a[q1][q]/b11;
						a22=a[q1][q1]/b22;
						a33=a[m1][m1]/b33;
						a34=a[m1][m]/b44;
						a43=a[m][m1]/b33;
						a44=a[m][m]/b44;
						b12=b[q][q1]/b22;
						b34=b[m1][m]/b44;
						a10=((a33-a11)*(a44-a11)-a34*a43+a43*b34*a11)/a21+
								a12-a11*b12;
						a20=(a22-a11-a21*b12)-(a33-a11)-(a44-a11)+a43*b34;
						a30=a[q+2][q1]/b22;
					}
					old1=fabs(a[m][m-1]);
					old2=fabs(a[m-1][m-2]);
					for (k=q; k<=m1; k++) {
						k1=k+1;
						k2=k+2;
						k3 = (k+3 > m) ? m : k+3;
						km1 = (k-1 < q) ? q : k-1;
						if (k != m1) {
							if (k == q)
								hsh3col(km1,km1,n,k,a10,a20,a30,a,b);
							else {
								hsh3col(km1,km1,n,k,a[k][km1],
											a[k1][km1],a[k2][km1],a,b);
								a[k1][km1]=a[k2][km1]=0.0;
							}
							hsh3row3(1,k3,n,k,b[k2][k2],b[k2][k1],
										b[k2][k],a,b,x);
							b[k2][k]=b[k2][k1]=0.0;
						} else {
							hsh2col(km1,km1,n,k,a[k][km1],a[k1][km1],a,b);
							a[k1][km1]=0.0;
						}
						hsh2row3(1,k3,k3,n,k,b[k1][k1],b[k1][k],a,b,x);
						b[k1][k]=0.0;
					}
				}
			}
		}	/* goon loop */
		if (out) break;
	} while (m >= 3);

	m=n;
	do {
		if ((m > 1) ? (a[m][m-1] == 0) : 1) {
			alfr[m]=a[m][m];
			beta[m]=b[m][m];
			alfi[m]=0.0;
			m--;
		} else {
			l=m-1;
			if (fabs(b[l][l]) <= epsb) {
				b[l][l]=0.0;
				hsh2col(l,l,n,l,a[l][l],a[m][l],a,b);
				a[m][l]=b[m][l]=0.0;
				alfr[l]=a[l][l];
				alfr[m]=a[m][m];
				beta[l]=b[l][l];
				beta[m]=b[m][m];
				alfi[m]=alfi[l]=0.0;
			} else
				if (fabs(b[m][m]) <= epsb) {
					b[m][m]=0.0;
					hsh2row3(1,m,m,n,l,a[m][m],a[m][l],a,b,x);
					a[m][l]=b[m][l]=0.0;
					alfr[l]=a[l][l];
					alfr[m]=a[m][m];
					beta[l]=b[l][l];
					beta[m]=b[m][m];
					alfi[m]=alfi[l]=0.0;
				} else {
					an=fabs(a[l][l])+fabs(a[l][m])+fabs(a[m][l])+
						fabs(a[m][m]);
					bn=fabs(b[l][l])+fabs(b[l][m])+fabs(b[m][m]);
					a11=a[l][l]/an;
					a12=a[l][m]/an;
					a21=a[m][l]/an;
					a22=a[m][m]/an;
					b11=b[l][l]/bn;
					b12=b[l][m]/bn;
					b22=b[m][m]/bn;
					e=a11/b11;
					c=((a22-e*b22)/b22-(a21*b12)/(b11*b22))/2.0;
					d=c*c+(a21*(a12-e*b12))/(b11*b22);
					if (d >= 0.0) {
						e += ((c < 0.0) ? c-sqrt(d) : c+sqrt(d));
						a11 -= e*b11;
						a12 -= e*b12;
						a22 -= e*b22;
						if (fabs(a11)+fabs(a12) >= fabs(a21)+fabs(a22))
							hsh2row3(1,m,m,n,l,a12,a11,a,b,x);
						else
							hsh2row3(1,m,m,n,l,a22,a21,a,b,x);
						if (an >= fabs(e)*bn)
							hsh2col(l,l,n,l,b[l][l],b[m][l],a,b);
						else
							hsh2col(l,l,n,l,a[l][l],a[m][l],a,b);
						a[m][l]=b[m][l]=0.0;
						alfr[l]=a[l][l];
						alfr[m]=a[m][m];
						beta[l]=b[l][l];
						beta[m]=b[m][m];
						alfi[m]=alfi[l]=0.0;
					} else {
						er=e+c;
						ei=sqrt(-d);
						a11r=a11-er*b11;
						a11i=ei*b11;
						a12r=a12-er*b12;
						a12i=ei*b12;
						a21r=a21;
						a21i=0.0;
						a22r=a22-er*b22;
						a22i=ei*b22;
						if (fabs(a11r)+fabs(a11i)+fabs(a12r)+fabs(a12i) >=
								fabs(a21r)+fabs(a22r)+fabs(a22i))
							chsh2(a12r,a12i,-a11r,-a11i,&cz,&szr,&szi);
						else
							chsh2(a22r,a22i,-a21r,-a21i,&cz,&szr,&szi);
						if (an >= (fabs(er)+fabs(ei))*bn)
							chsh2(cz*b11+szr*b12,szi*b12,szr*b22,szi*b22,
									&cq,&sqr,&sqi);
						else
							chsh2(cz*a11+szr*a12,szi*a12,cz*a21+szr*a22,
									szi*a22,&cq,&sqr,&sqi);
						ssr=sqr*szr+sqi*szi;
						ssi=sqr*szi-sqi*szr;
						tr=cq*cz*a11+cq*szr*a12+sqr*cz*a21+ssr*a22;
						ti=cq*szi*a12-sqi*cz*a21+ssi*a22;
						bdr=cq*cz*b11+cq*szr*b12+ssr*b22;
						bdi=cq*szi*b12+ssi*b22;
						r=sqrt(bdr*bdr+bdi*bdi);
						beta[l]=bn*r;
						alfr[l]=an*(tr*bdr+ti*bdi)/r;
						alfi[l]=an*(tr*bdi-ti*bdr)/r;
						tr=ssr*a11-sqr*cz*a12-cq*szr*a21+cq*cz*a22;
						ti = -ssi*a11-sqi*cz*a12+cq*szi*a21;
						bdr=ssr*b11-sqr*cz*b12+cq*cz*b22;
						bdi = -ssi*b11-sqi*cz*b12;
						r=sqrt(bdr*bdr+bdi*bdi);
						beta[m]=bn*r;
						alfr[m]=an*(tr*bdr+ti*bdi)/r;
						alfi[m]=an*(tr*bdi-ti*bdr)/r;
					}
				}
			m -= 2;
		}
	} while (m > 0);

	for (m=n; m>=1; m--)
		if (alfi[m] == 0.0) {
			alfm=alfr[m];
			betm=beta[m];
			b[m][m]=1.0;
			l1=m;
			for (l=m-1; l>=1; l--) {
				sl=0.0;
				for (j=l1; j<=m; j++)
					sl += (betm*a[l][j]-alfm*b[l][j])*b[j][m];
				if ((l != 1) ? (betm*a[l][l-1] == 0.0) : 1) {
					d=betm*a[l][l]-alfm*b[l][l];
					if (d == 0.0) d=(epsa+epsb)/2.0;
					b[l][m] = -sl/d;
				} else {
					k=l-1;
					sk=0.0;
					for (j=l1; j<=m; j++)
						sk += (betm*a[k][j]-alfm*b[k][j])*b[j][m];
					tkk=betm*a[k][k]-alfm*b[k][k];
					tkl=betm*a[k][l]-alfm*b[k][l];
					tlk=betm*a[l][k];
					tll=betm*a[l][l]-alfm*b[l][l];
					d=tkk*tll-tkl*tlk;
					if (d == 0.0) d=(epsa+epsb)/2.0;
					b[l][m]=(tlk*sk-tkk*sl)/d;
					b[k][m] = (fabs(tkk) >= fabs(tlk)) ?
								-(sk+tkl*b[l][m])/tkk :
								-(sl+tll*b[l][m])/tlk;
					l--;
				}
				l1=l;
			}
		} else {
			almr=alfr[m-1];
			almi=alfi[m-1];
			betm=beta[m-1];
			mr=m-1;
			mi=m;
			b[m-1][mr]=almi*b[m][m]/(betm*a[m][m-1]);
			b[m-1][mi]=(betm*a[m][m]-almr*b[m][m])/(betm*a[m][m-1]);
			b[m][mr]=0.0;
			b[m][mi] = -1.0;
			l1=m-1;
			for (l=m-2; l>=1; l--) {
				slr=sli=0.0;
				for (j=l1; j<=m; j++) {
					tr=betm*a[l][j]-almr*b[l][j];
					ti = -almi*b[l][j];
					slr += tr*b[j][mr]-ti*b[j][mi];
					sli += tr*b[j][mi]+ti*b[j][mr];
				}
				if ((l != 1) ? (betm*a[l][l-1] == 0.0) : 1) {
					dr=betm*a[l][l]-almr*b[l][l];
					di = -almi*b[l][l];
					comdiv(-slr,-sli,dr,di,&b[l][mr],&b[l][mi]);
				} else {
					k=l-1;
					skr=ski=0.0;
					for (j=l1; j<=m; j++) {
						tr=betm*a[k][j]-almr*b[k][j];
						ti = -almi*b[k][j];
						skr += tr*b[j][mr]-ti*b[j][mi];
						ski += tr*b[j][mi]+ti*b[j][mr];
					}
					tkkr=betm*a[k][k]-almr*b[k][k];
					tkki = -almi*b[k][k];
					tklr=betm*a[k][l]-almr*b[k][l];
					tkli = -almi*b[k][l];
					tlkr=betm*a[l][k];
					tlki=0.0;
					tllr=betm*a[l][l]-almr*b[l][l];
					tlli = -almi*b[l][l];
					dr=tkkr*tllr-tkki*tlli-tklr*tlkr;
					di=tkkr*tlli+tkki*tllr-tkli*tlkr;
					if (dr == 0.0 && di == 0.0) dr=(epsa+epsb)/2.0;
					comdiv(tlkr*skr-tkkr*slr+tkki*sli,
							tlkr*ski-tkkr*sli-tkki*slr,
							dr,di,&b[l][mr],&b[l][mi]);
					if (fabs(tkkr)+fabs(tkki) >= fabs(tlkr))
						comdiv(-skr-tklr*b[l][mr]+tkli*b[l][mi],
								-ski-tklr*b[l][mi]-tkli*b[l][mr],
								tkkr,tkki,&b[k][mr],&b[k][mi]);
					else
						comdiv(-slr-tllr*b[l][mr]+tlli*b[l][mi],
								-sli-tllr*b[l][mi]-tlli*b[l][mr],
								tlkr,tlki,&b[k][mr],&b[k][mi]);
					l--;
				}
				l1=l;
			}
			m--;
		}
	for (m=n; m>=1; m--)
		for (k=1; k<=n; k++) x[k][m]=matmat(1,m,k,m,x,b);
	for (m=n; m>=1; m--) {
		s=0.0;
		if (alfi[m] == 0.0) {
			for (k=1; k<=n; k++) {
				r=fabs(x[k][m]);
				if (r >= s) {
					s=r;
					d=x[k][m];
				}
			}
			for (k=1; k<=n; k++) x[k][m] /= d;
		} else {
			for (k=1; k<=n; k++) {
				r=fabs(x[k][m-1])+fabs(x[k][m]);
				an=x[k][m-1]/r;
				bn=x[k][m]/r;
				r *= sqrt(an*an+bn*bn);
				if (r >= s) {
					s=r;
					dr=x[k][m-1];
					di=x[k][m];
				}
			}
			for (k=1; k<=n; k++)
				comdiv(x[k][m-1],x[k][m],dr,di,&x[k][m-1],&x[k][m]);
			m--;
		}
	}
}
