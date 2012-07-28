/* ARK.C */
void ark(real_t *t, real_t *te, int *m0, int *m, real_t u[], void (*derivative)(int *, int *, real_t *, real_t []), real_t data[], void (*out)(int *, int *, real_t *, real_t *, real_t [], real_t []));
real_t arkmui(int i, int n, int p, real_t lambda[]);
real_t arklabda(int i, int j, int n, int p, real_t lambda[]);
/* ARKMAT.C */
void arkmat(real_t *t, real_t te, int m, int n, real_t **u, void (*der)(int, int, real_t, real_t **, real_t **), int type, int *order, real_t *spr, void (*out)(real_t, real_t, int, int, real_t **, int, int, real_t *));
/* DAVUPD.C */
void davupd(real_t h[], int n, real_t v[], real_t w[], real_t c1, real_t c2);
/* DIFFSYS.C */
void diffsys(real_t *x, real_t xe, int n, real_t y[], void (*derivative)(int, real_t, real_t [], real_t []), real_t aeta, real_t reta, real_t s[], real_t h0, void (*output)(int, real_t, real_t, real_t [], real_t []));
/* EFERK.C */
void eferk(real_t *x, real_t xe, int m, real_t y[], real_t *sigma, real_t phi, void (*derivative)(int, real_t []), real_t **j, void (*jacobian)(int, real_t **, real_t [], real_t *), int *k, int l, int aut, real_t aeta, real_t reta, real_t hmin, real_t hmax, int linear, void (*output)(real_t, real_t, int, real_t [], real_t **, int));
/* EFRK.C */
void efrk(real_t *t, real_t te, int m0, int m, real_t u[], real_t *sigma, real_t *phi, real_t *diameter, void (*derivative)(int, int, real_t, real_t []), int *k, real_t *step, real_t r, real_t l, real_t beta[], int thirdorder, real_t tol, void (*output)(int, int, real_t, real_t, real_t [], real_t *, real_t *, real_t *, int, real_t *, int, int));
/* EFSIRK.C */
void efsirk(real_t *x, real_t xe, int m, real_t y[], real_t *delta, void (*derivative)(int, real_t [], real_t *), void (*jacobian)(int, real_t **, real_t [], real_t *), real_t **j, int *n, real_t aeta, real_t reta, real_t hmin, real_t hmax, int linear, void (*output)(real_t, real_t, int, real_t [], real_t, real_t **, int));
/* EFT.C */
void eft(real_t *t, real_t te, int m0, int m, real_t u[], real_t (*sigma)(real_t, int, int), real_t phi, real_t (*diameter)(real_t, int, int), void (*derivative)(real_t, int, int, int, real_t []), int *k, real_t alfa, int norm, real_t (*aeta)(real_t, int, int), real_t (*reta)(real_t, int, int), real_t *eta, real_t *rho, real_t hmin, real_t *hstart, void (*out)(real_t, real_t, int, int, real_t [], int, real_t, real_t));
/* ELIMINAT.C */
void elimination(real_t **u, int lj, int uj, int ll, int ul, void (*residual)(int, int, int, int, real_t **), real_t a, real_t b, int *n, real_t discr[], int *k, real_t *rateconv, real_t *domeigval, void (*out)(real_t **, int, int, int, int, int *, real_t [], int, real_t, real_t));
real_t optpol(real_t x, real_t a, real_t b, real_t domeigval);
/* FEMHERMS.C */
void femhermsym(real_t x[], real_t y[], int n, real_t (*p)(real_t), real_t (*q)(real_t), real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
void femhermsymeval(int order, int l, real_t (*p)(real_t), real_t (*q)(real_t), real_t (*r)(real_t), real_t (*f)(real_t), real_t *a11, real_t *a12, real_t *a13, real_t *a14, real_t *a22, real_t *a23, real_t *a24, real_t *a33, real_t *a34, real_t *a44, real_t *b1, real_t *b2, real_t *b3, real_t *b4, real_t *xl, real_t xl1);
/* FEMLAG.C */
void femlag(real_t x[], real_t y[], int n, real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FEMLAGSK.C */
void femlagskew(real_t x[], real_t y[], int n, real_t (*q)(real_t), real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FEMLAGSP.C */
void femlagspher(real_t x[], real_t y[], int n, int nc, real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FEMLAGSY.C */
void femlagsym(real_t x[], real_t y[], int n, real_t (*p)(real_t), real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FLEMIN.C */
real_t flemin(int n, real_t x[], real_t g[], real_t h[], real_t (*funct)(int, real_t [], real_t []), real_t in[], real_t out[]);
/* FLEUPD.C */
void fleupd(real_t h[], int n, real_t v[], real_t w[], real_t c1, real_t c2);
/* GMS.C */
void gms(real_t *x, real_t xe, int r, real_t y[], real_t h, real_t hmin, real_t hmax, real_t *delta, void (*derivative)(int, real_t [], real_t *), void (*jacobian)(int, real_t **, real_t [], real_t *), real_t aeta, real_t reta, int *n, int *jev, int *lu, int nsjev, int linear, void (*out)(real_t, real_t, int, real_t [], real_t, int, int, int));
void gmsopconstruct(int *reeval, int *update, int r, real_t **hjac, real_t **h2jac2, real_t **rqz, real_t y[], real_t aux[], int ri[], int ci[], int *lu, int *jev, int *nsjev1, real_t *delta, real_t *alfa, real_t h1, real_t h0, real_t *s1, real_t *s2, void (*jacobian)(int, real_t **, real_t [], real_t *));
void gmscoefficient(real_t *xl1, real_t *xl0, real_t x0, int change, int n, real_t *q1, real_t *q2, real_t h1, real_t alfa, real_t **bd1, real_t **bd2, int strategy);
void gmsdiffscheme(int k, int count, int r, real_t fl[], real_t yl[], int *n, int *nsjev1, real_t y0[], real_t alfa, real_t **bd1, real_t **bd2, real_t h1, real_t y[], real_t **hjac, real_t **h2jac2, real_t **rqz, int ri[], int ci[], real_t *delta, void (*derivative)(int, real_t [], real_t *));
/* GSSNEWTO.C */
void gssnewton(int m, int n, real_t par[], real_t rv[], real_t **jjinv, int (*funct)(int, int, real_t [], real_t []), void (*jacobian)(int, int, real_t [], real_t [], real_t **), real_t in[], real_t out[]);
/* IMPEX.C */
void impex(int n, real_t t0, real_t tend, real_t y0[], void (*deriv)(real_t, real_t [], real_t [], int), int (*available)(real_t, real_t [], real_t **, int), real_t h0, real_t hmax, int presch, real_t eps, real_t weights[], void (*update)(real_t [], real_t [], int), int *fail, void (*control)(real_t *, real_t, real_t, real_t, real_t **, real_t [], int, real_t));
int impexrecomp(real_t **a, real_t h, real_t t, real_t y[], int ps[], int n, int (*available)(real_t, real_t [], real_t **, int), void (*deriv)(real_t, real_t [], real_t [], int));
int impexlargestep(int n, real_t y[], real_t t, real_t *t1, real_t *t2, real_t *t3, real_t s1[], real_t s2[], real_t s3[], real_t h, real_t h2, real_t z[], real_t u1[], real_t u3[], real_t w1[], real_t w2[], real_t w3[], int ps1[], int ps2[], real_t weights[], real_t **a1, real_t **a2, real_t eps, void (*deriv)(real_t, real_t [], real_t [], int), int (*available)(real_t, real_t [], real_t **, int));
int impexiterate(real_t z[], real_t y[], real_t **a, real_t h, real_t t, real_t weights[], int ps[], int n, real_t eps, void (*deriv)(real_t, real_t [], real_t [], int), int (*available)(real_t, real_t [], real_t **, int));
void impexbackdiff(int n, real_t u1[], real_t u3[], real_t w1[], real_t w2[], real_t w3[], real_t s1[], real_t s2[], real_t s3[], real_t **r, real_t **rf);
/* LINEMIN.C */
void linemin(int n, real_t x[], real_t d[], real_t nd, real_t *alfa, real_t g[], real_t (*funct)(int, real_t [], real_t []), real_t f0, real_t *f1, real_t df0, real_t *df1, int *evlmax, int strongsearch, real_t in[]);
/* LINIGER1.C */
void liniger1vs(real_t *x, real_t xe, int m, real_t y[], real_t *sigma, void (*derivative)(int, real_t [], real_t *), real_t **j, void (*jacobian)(int, real_t **, real_t [], real_t *), int itmax, real_t hmin, real_t hmax, real_t aeta, real_t reta, real_t info[], void (*output)(real_t, real_t, int, real_t [], real_t, real_t **, real_t []));
void liniger1vscoef(int m, real_t **a, real_t **j, real_t aux[], int pi[], real_t h, real_t sigma, real_t *mu, real_t *mu1, real_t *beta, real_t *p);
/* LINIGER2.C */
void liniger2(real_t *x, real_t xe, int m, real_t y[], real_t *sigma1, real_t *sigma2, real_t (*f)(int, real_t [], int, real_t *, real_t *), int (*evaluate)(int), real_t **j, void (*jacobian)(int, real_t **, real_t [], real_t *, real_t *), int *k, int itmax, real_t step, real_t aeta, real_t reta, void (*output)(real_t, real_t, int, real_t [], real_t, real_t, real_t **, int));
void liniger2coef(int m, real_t **j, real_t **a, real_t aux[], int pi[], real_t h, real_t sigma1, real_t sigma2, real_t *c0, real_t *c1, real_t *c2, real_t *c3, real_t *c4);
/* MARQUARD.C */
void marquardt(int m, int n, real_t par[], real_t g[], real_t **v, int (*funct)(int, int, real_t [], real_t []), void (*jacobian)(int, int, real_t [], real_t [], real_t **), real_t in[], real_t out[]);
/* MININ.C */
real_t minin(real_t *x, real_t *a, real_t *b, real_t (*fx)(real_t), real_t (*tolx)(real_t));
/* MININDER.C */
real_t mininder(real_t *x, real_t *y, real_t (*fx)(real_t), real_t (*dfx)(real_t), real_t (*tolx)(real_t));
/* MTAYLOR.C */
void modifiedtaylor(real_t *t, real_t te, int m0, int m, real_t u[], real_t (*sigma)(real_t, int, int), real_t taumin, void (*derivative)(real_t, int, int, int, real_t []), int *k, real_t data[], real_t alfa, int norm, real_t (*aeta)(real_t, int, int), real_t (*reta)(real_t, int, int), real_t *eta, real_t *rho, void (*out)(real_t, real_t, int, int, real_t [], int, real_t, real_t));
/* MULTISTE.C */
int multistep(real_t *x, real_t xend, real_t y[], real_t hmin, real_t hmax, real_t ymax[], real_t eps, int *first, real_t save[], void (*deriv)(real_t [], int, real_t, real_t []), int (*available)(int, real_t, real_t [], real_t **), real_t **jacobian, int stiff, int n, void (*out)(real_t, int, int, real_t, real_t []));
void multistepreset(real_t y[], real_t save[], real_t *x, real_t *ch, real_t *c, real_t *h, int *decomposed, real_t hmin, real_t hmax, real_t hold, real_t xold, int m, int k, int n);
void multisteporder(real_t a[], real_t save[], real_t *tolup, real_t *tol, real_t *toldwn, real_t *tolconv, real_t *a0, int *decompose, real_t eps, int k, int n);
void multistepstep(int *knew, real_t *chnew, real_t tolup, real_t tol, real_t toldwn, real_t delta[], real_t error, real_t lastdelta[], real_t y[], real_t ymax[], int fails, int m, int k, int n);
void multistepjacobian(int n, real_t x, real_t y[], real_t eps, real_t fixy[], real_t fixdy[], real_t dy[], real_t **jacobian, void (*deriv)(real_t [], int, real_t, real_t []), int (*available)(int, real_t, real_t [], real_t **), int *evaluate, int *decompose, int *evaluated);
/* NONLINFE.C */
void nonlinfemlagskew(real_t x[], real_t y[], int n, real_t (*f)(real_t, real_t, real_t), real_t (*fy)(real_t, real_t, real_t), real_t (*fz)(real_t, real_t, real_t), int nc, real_t e[]);
/* PEIDE.C */
void peide(int n, int m, int nobs, int *nbp, real_t par[], real_t res[], int bp[], real_t **jtjinv, real_t in[], real_t out[], int (*deriv)(int, int, real_t [], real_t [], real_t, real_t []), int (*jacdfdy)(int, int, real_t [], real_t [], real_t, real_t **), int (*jacdfdp)(int, int, real_t [], real_t [], real_t, real_t **), void (*callystart)(int, int, real_t [], real_t [], real_t []), void (*data)(int, real_t [], real_t [], int []), void (*monitor)(int, int, int, real_t [], real_t [], int, int));
int peidefunct(int nrow, int ncol, real_t par[], real_t res[], int n, int m, int nobs, int *nbp, int first, int *sec, int *max, int *nis, real_t eps1, int weight, int bp[], real_t save[], real_t ymax[], real_t y[], real_t **yp, real_t **fy, real_t **fp, int cobs[], real_t tobs[], real_t obs[], real_t in[], real_t aux[], int clean, int (*deriv)(int, int, real_t [], real_t [], real_t, real_t []), int (*jacdfdy)(int, int, real_t [], real_t [], real_t, real_t **), int (*jacdfdp)(int, int, real_t [], real_t [], real_t, real_t **), void (*callystart)(int, int, real_t [], real_t [], real_t []), void (*monitor)(int, int, int, real_t [], real_t [], int, int));
void peidereset(int n, int k, real_t hmin, real_t hmax, real_t hold, real_t xold, real_t y[], real_t save[], real_t *ch, real_t *x, real_t *h, int *decompose);
void peideorder(int n, int k, real_t eps, real_t a[], real_t save[], real_t *tol, real_t *tolup, real_t *toldwn, real_t *tolconv, real_t *a0, int *decompose);
void peidestep(int n, int k, int fails, real_t tolup, real_t toldwn, real_t tol, real_t error, real_t delta[], real_t lastdelta[], real_t y[], real_t ymax[], int *knew, real_t *chnew);
real_t peideinterpol(int startindex, int jump, int k, real_t tobsdif, real_t y[]);
/* PRAXIS.C */
void praxis(int n, real_t x[], real_t (*funct)(int, real_t []), real_t in[], real_t out[]);
void praxismin(int j, int nits, real_t *d2, real_t *x1, real_t *f1, int fk, int n, real_t x[], real_t **v, real_t *qa, real_t *qb, real_t *qc, real_t qd0, real_t qd1, real_t q0[], real_t q1[], int *nf, int *nl, real_t *fx, real_t m2, real_t m4, real_t dmin, real_t ldt, real_t reltol, real_t abstol, real_t small, real_t h, real_t (*funct)(int, real_t []));
real_t praxisflin(real_t l, int j, int n, real_t x[], real_t **v, real_t *qa, real_t *qb, real_t *qc, real_t qd0, real_t qd1, real_t q0[], real_t q1[], int *nf, real_t (*funct)(int, real_t []));
/* QUANEWB1.C */
void quanewbnd1(int n, int lw, int rw, real_t x[], real_t f[], int (*funct)(int, int, int, real_t [], real_t []), real_t in[], real_t out[]);
real_t quanewbnd1s(int i);
real_t quanewbnd1t(real_t x, int i);
/* QUANEWBN.C */
void quanewbnd(int n, int lw, int rw, real_t x[], real_t f[], real_t jac[], int (*funct)(int, int, int, real_t [], real_t []), real_t in[], real_t out[]);
/* RICHARDS.C */
void richardson(real_t **u, int lj, int uj, int ll, int ul, int inap, void (*residual)(int, int, int, int, real_t **), real_t a, real_t b, int *n, real_t discr[], int *k, real_t *rateconv, real_t *domeigval, void (*out)(real_t **, int, int, int, int, int *, real_t [], int, real_t, real_t));
/* RK1.C */
void rk1(real_t *x, real_t a, real_t b, real_t *y, real_t ya, real_t (*fxy)(real_t, real_t), real_t e[], real_t d[], int fi);
/* RK2.C */
void rk2(real_t *x, real_t a, real_t b, real_t *y, real_t ya, real_t *z, real_t za, real_t (*fxyz)(real_t, real_t, real_t), real_t e[], real_t d[], int fi);
/* RK2N.C */
void rk2n(real_t *x, real_t a, real_t b, real_t y[], real_t ya[], real_t z[], real_t za[], real_t (*fxyzj)(int, int, real_t, real_t [], real_t []), real_t e[], real_t d[], int fi, int n);
/* RK3.C */
void rk3(real_t *x, real_t a, real_t b, real_t *y, real_t ya, real_t *z, real_t za, real_t (*fxy)(real_t, real_t), real_t e[], real_t d[], int fi);
/* RK3N.C */
void rk3n(real_t *x, real_t a, real_t b, real_t y[], real_t ya[], real_t z[], real_t za[], real_t (*fxyj)(int, int, real_t, real_t []), real_t e[], real_t d[], int fi, int n);
/* RK4A.C */
void rk4a(real_t *x, real_t xa, real_t (*b)(real_t, real_t), real_t *y, real_t ya, real_t (*fxy)(real_t, real_t), real_t e[], real_t d[], int fi, int xdir, int pos);
void rk4arkstep(real_t *x, real_t xl, real_t h, real_t *y, real_t yl, real_t zl, real_t (*fxy)(real_t, real_t), int d, int invf, real_t *k0, real_t *k1, real_t *k2, real_t *k3, real_t *k4, real_t *k5, real_t *discr, real_t mu);
/* RK4NA.C */
void rk4na(real_t x[], real_t xa[], real_t (*b)(int, real_t []), real_t (*fxj)(int, int, real_t []), real_t e[], real_t d[], int fi, int n, int l, int pos);
void rk4narkstep(real_t h, int d, int n, int iv, real_t mu, real_t (*fxj)(int, int, real_t []), real_t x[], real_t xl[], real_t y[], real_t discr[], real_t **k);
/* RK5NA.C */
void rk5na(real_t x[], real_t xa[], real_t (*b)(int, real_t []), real_t (*fxj)(int, int, real_t []), real_t e[], real_t d[], int fi, int n, int l, int pos);
void rk5narkstep(real_t h, int d, int n, real_t mu, real_t (*fxj)(int, int, real_t []), real_t x[], real_t xl[], real_t y[], real_t discr[], real_t **k);
/* RKE.C */
void rke(real_t *x, real_t *xe, int n, real_t y[], void (*der)(int, real_t, real_t []), real_t data[], int fi, void (*out)(int, real_t, real_t, real_t [], real_t []));
/* RNK1MIN.C */
real_t rnk1min(int n, real_t x[], real_t g[], real_t h[], real_t (*funct)(int, real_t [], real_t []), real_t in[], real_t out[]);
/* RNK1UPD.C */
void rnk1upd(real_t h[], int n, real_t v[], real_t c);
/* ZEROIN.C */
int zeroin(real_t *x, real_t *y, real_t (*fx)(real_t), real_t (*tolx)(real_t));
/* ZEROINDE.C */
int zeroinder(real_t *x, real_t *y, real_t (*fx)(real_t), real_t (*dfx)(real_t), real_t (*tolx)(real_t));
/* ZEROINRA.C */
int zeroinrat(real_t *x, real_t *y, real_t (*fx)(real_t), real_t (*tolx)(real_t));
