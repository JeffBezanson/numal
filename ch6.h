/* AIRY.C */
void airy(real_t z, real_t *ai, real_t *aid, real_t *bi, real_t *bid, real_t *expon, int first);
/* AIRYZERO.C */
real_t airyzeros(int n, int d, real_t zai[], real_t vai[]);
/* ARCCOSH.C */
real_t arccosh(real_t x);
/* ARCSINH.C */
real_t arcsinh(real_t x);
/* ARCTANH.C */
real_t arctanh(real_t x);
/* BACKWARD.C */
void backward(real_t x, real_t p, real_t q, real_t i0, int nmax, real_t eps, real_t i[]);
/* BESPQA01.C */
void besspqa01(real_t a, real_t x, real_t *pa, real_t *qa, real_t *pa1, real_t *qa1);
/* BESSI0.C */
real_t bessi0(real_t x);
/* BESSI1.C */
real_t bessi1(real_t x);
/* BESSIAPL.C */
void bessiaplusn(real_t a, real_t x, int n, real_t ia[]);
/* BESSI.C */
void bessi(real_t x, int n, real_t i[]);
/* BESSJ0.C */
real_t bessj0(real_t x);
/* BESSJ1.C */
real_t bessj1(real_t x);
/* BESSJAPL.C */
void bessjaplusn(real_t a, real_t x, int n, real_t ja[]);
/* BESSJ.C */
void bessj(real_t x, int n, real_t j[]);
/* BESSK01.C */
void bessk01(real_t x, real_t *k0, real_t *k1);
/* BESSKA01.C */
void besska01(real_t a, real_t x, real_t *ka, real_t *ka1);
/* BESSKAPL.C */
void besskaplusn(real_t a, real_t x, int nmax, real_t kan[]);
/* BESSK.C */
void bessk(real_t x, int n, real_t k[]);
/* BESSPQ0.C */
void besspq0(real_t x, real_t *p0, real_t *q0);
/* BESSPQ1.C */
void besspq1(real_t x, real_t *p1, real_t *q1);
/* BESSY01.C */
void bessy01(real_t x, real_t *y0, real_t *y1);
/* BESSYA01.C */
void bessya01(real_t a, real_t x, real_t *ya, real_t *ya1);
/* BESSYAPL.C */
void bessyaplusn(real_t a, real_t x, int nmax, real_t yan[]);
/* BESSY.C */
void bessy(real_t x, int n, real_t y[]);
/* BESSZERO.C */
void besszeros(real_t a, int n, real_t z[], int d);
/* EIALPHA.C */
void eialpha(real_t x, int n, real_t alpha[]);
/* EI.C */
real_t ei(real_t x);
/* ENX.C */
void enx(real_t x, int n1, int n2, real_t a[]);
/* ERRORFUN.C */
void errorfunction(real_t x, real_t *erf, real_t *erfc);
/* FG.C */
void fg(real_t x, real_t *f, real_t *g);
/* FORWARD.C */
void forward(real_t x, real_t p, real_t q, real_t i0, real_t i1, int nmax, real_t i[]);
/* FRESNEL.C */
void fresnel(real_t x, real_t *c, real_t *s);
/* GAMMA.C */
real_t numal_gamma(real_t x);
/* IBPPLUSN.C */
void ibpplusn(real_t x, real_t p, real_t q, int nmax, real_t eps, real_t i[]);
/* IBQPLUSN.C */
void ibqplusn(real_t x, real_t p, real_t q, int nmax, real_t eps, real_t i[]);
/* INCBETA.C */
real_t incbeta(real_t x, real_t p, real_t q, real_t eps);
/* INCOMGAM.C */
void incomgam(real_t x, real_t a, real_t *klgam, real_t *grgam, real_t gam, real_t eps);
/* INVERRFN.C */
void inverseerrorfunction(real_t x, real_t oneminx, real_t *inverf);
/* IXPFIX.C */
void ixpfix(real_t x, real_t p, real_t q, int nmax, real_t eps, real_t i[]);
/* IXQFIX.C */
void ixqfix(real_t x, real_t p, real_t q, int nmax, real_t eps, real_t i[]);
/* LOGGAMMA.C */
real_t loggamma(real_t x);
/* LOGONEPL.C */
real_t logoneplusx(real_t x);
/* NEBESIAP.C */
void nonexpbessiaplusn(real_t a, real_t x, int n, real_t ia[]);
/* NEBESK01.C */
void nonexpbessk01(real_t x, real_t *k0, real_t *k1);
/* NEBESSI0.C */
real_t nonexpbessi0(real_t x);
/* NEBESSI1.C */
real_t nonexpbessi1(real_t x);
/* NEBESSI.C */
void nonexpbessi(real_t x, int n, real_t i[]);
/* NEBESSK.C */
void nonexpbessk(real_t x, int n, real_t k[]);
/* NEBSKA01.C */
void nonexpbesska01(real_t a, real_t x, real_t *ka, real_t *ka1);
/* NEBSKAPL.C */
void nonexpbesskaplusn(real_t a, real_t x, int nmax, real_t kan[]);
/* NESPBESI.C */
void nonexpspherbessi(real_t x, int n, real_t i[]);
/* NESPBESK.C */
void nonexpspherbessk(real_t x, int n, real_t k[]);
/* NONEXPEN.C */
void nonexpenx(real_t x, int n1, int n2, real_t a[]);
/* NONEXPER.C */
real_t nonexperfc(real_t x);
/* RECIPGAM.C */
real_t recipgamma(real_t x, real_t *odd, real_t *even);
/* SINCOSFG.C */
void sincosfg(real_t x, real_t *f, real_t *g);
/* SINCOSIN.C */
void sincosint(real_t x, real_t *si, real_t *ci);
/* SPBESSI.C */
void spherbessi(real_t x, int n, real_t i[]);
/* SPBESSJ.C */
void spherbessj(real_t x, int n, real_t j[]);
/* SPBESSK.C */
void spherbessk(real_t x, int n, real_t k[]);
/* SPBESSY.C */
void spherbessy(real_t x, int n, real_t y[]);
/* START.C */
int start(real_t x, int n, int t);
