/* EULER.C */
real_t euler(real_t (*ai)(int), real_t eps, int tim);
/* GSSJACWG.C */
void gssjacwghts(int n, real_t alfa, real_t beta, real_t x[], real_t w[]);
/* GSSLAGWG.C */
void gsslagwghts(int n, real_t alfa, real_t x[], real_t w[]);
/* GSSWTS.C */
void gsswts(int n, real_t zer[], real_t b[], real_t c[], real_t w[]);
/* GSSWTSSY.C */
void gsswtssym(int n, real_t zer[], real_t c[], real_t w[]);
/* INTEGRAL.C */
real_t integral(real_t a, real_t b, real_t (*fx)(real_t), real_t e[], int ua, int ub);
real_t integralqad(int transf, real_t (*fx)(real_t), real_t e[], real_t *x0, real_t *x1, real_t *x2, real_t *f0, real_t *f1, real_t *f2, real_t re, real_t ae, real_t b1);
void integralint(int transf, real_t (*fx)(real_t), real_t e[], real_t *x0, real_t *x1, real_t *x2, real_t *f0, real_t *f1, real_t *f2, real_t *sum, real_t re, real_t ae, real_t b1, real_t hmin);
/* JACOBNBN.C */
void jacobnbndf(int n, int lw, int rw, real_t x[], real_t f[], real_t jac[], real_t (*di)(int), int (*funct)(int, int, int, real_t [], real_t []));
/* JACOBNMF.C */
void jacobnmf(int n, int m, real_t x[], real_t f[], real_t **jac, real_t (*di)(int), void (*funct)(int, int, real_t [], real_t []));
/* JACOBNNF.C */
void jacobnnf(int n, real_t x[], real_t f[], real_t **jac, real_t (*di)(int), void (*funct)(int, real_t [], real_t []));
/* QADRAT.C */
real_t qadrat(real_t *x, real_t a, real_t b, real_t (*fx)(real_t), real_t e[]);
real_t lint(real_t *x, real_t (*fx)(real_t), real_t e[], real_t x0, real_t xn, real_t f0, real_t f2, real_t f3, real_t f5, real_t f6, real_t f7, real_t f9, real_t f14, real_t hmin, real_t hmax, real_t re, real_t ae);
/* RECCOF.C */
void reccof(int n, int m, real_t *x, real_t (*wx)(real_t), real_t b[], real_t c[], real_t l[], int sym);
/* SUMPOSSE.C */
real_t sumposseries(real_t (*ai)(real_t), int maxaddup, real_t maxzero, int maxrecurs, int machexp, int tim);
real_t sumposseriessumup(int bjk, real_t (*ai)(real_t), int maxaddup, real_t maxzero, int maxrecurs, int machexp, int tim, int recurs, int vl, int vl2, int vl4, int jj);
real_t sumposseriesbjk(int j, real_t i, real_t machexp, real_t (*ai)(real_t));
/* TRICUB.C */
real_t tricub(real_t xi, real_t yi, real_t xj, real_t yj, real_t xk, real_t yk, real_t (*g)(real_t, real_t), real_t re, real_t ae);
real_t tricubint(real_t ax1, real_t ay1, real_t af1, real_t ax2, real_t ay2, real_t af2, real_t ax3, real_t ay3, real_t af3, real_t bx1, real_t by1, real_t bf1, real_t bx2, real_t by2, real_t bf2, real_t bx3, real_t by3, real_t bf3, real_t px, real_t py, real_t pf, real_t (*g)(real_t, real_t), real_t re, real_t ae, real_t *surf, real_t surfmin);
