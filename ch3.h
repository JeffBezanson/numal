/* ALLJACZE.C */
void alljaczer(int n, real_t alfa, real_t beta, real_t zer[]);
/* ALLLAGZE.C */
void alllagzer(int n, real_t alfa, real_t zer[]);
/* ALLZEROP.C */
void allzerortpol(int n, real_t b[], real_t c[], real_t zer[], real_t em[]);
/* BAKCOMHE.C */
void bakcomhes(real_t **ar, real_t **ai, real_t tr[], real_t ti[], real_t del[], real_t **vr, real_t **vi, int n, int n1, int n2);
/* BAKHRMTR.C */
void bakhrmtri(real_t **a, int n, int n1, int n2, real_t **vecr, real_t **veci, real_t tr[], real_t ti[]);
/* BAKLBR.C */
void baklbr(int n, int n1, int n2, real_t d[], int inter[], real_t **vec);
/* BAKLBRCO.C */
void baklbrcom(int n, int n1, int n2, real_t d[], int inter[], real_t **vr, real_t **vi);
/* BAKREAH1.C */
void bakreahes1(real_t **a, int n, int index[], real_t v[]);
/* BAKREAH2.C */
void bakreahes2(real_t **a, int n, int n1, int n2, int index[], real_t **vec);
/* BAKSYMT1.C */
void baksymtri1(real_t a[], int n, int n1, int n2, real_t **vec);
/* BAKSYMT2.C */
void baksymtri2(real_t **a, int n, int n1, int n2, real_t **vec);
/* BOUNDS.C */
void bounds(int n, real_t a[], real_t re[], real_t im[], real_t rele, real_t abse, real_t recentre[], real_t imcentre[], real_t bound[]);
void kcluster(int k, int m, int n, real_t rc[], real_t re[], real_t im[], real_t recentre[], real_t imcentre[], real_t bound[], real_t clust[]);
/* CHLDEC1.C */
void chldec1(real_t a[], int n, real_t aux[]);
/* CHLDEC2.C */
void chldec2(real_t **a, int n, real_t aux[]);
/* CHLDECBN.C */
void chldecbnd(real_t a[], int n, int w, real_t aux[]);
/* CHLDECI1.C */
void chldecinv1(real_t a[], int n, real_t aux[]);
/* CHLDECI2.C */
void chldecinv2(real_t **a, int n, real_t aux[]);
/* CHLDECS1.C */
void chldecsol1(real_t a[], int n, real_t aux[], real_t b[]);
/* CHLDECS2.C */
void chldecsol2(real_t **a, int n, real_t aux[], real_t b[]);
/* CHLDECSB.C */
void chldecsolbnd(real_t a[], int n, int w, real_t aux[], real_t b[]);
/* CHLDETM1.C */
real_t chldeterm1(real_t a[], int n);
/* CHLDETM2.C */
real_t chldeterm2(real_t **a, int n);
/* CHLDTMBN.C */
real_t chldetermbnd(real_t a[], int n, int w);
/* CHLINV1.C */
void chlinv1(real_t a[], int n);
/* CHLINV2.C */
void chlinv2(real_t **a, int n);
/* CHLSOL1.C */
void chlsol1(real_t a[], int n, real_t b[]);
/* CHLSOL2.C */
void chlsol2(real_t **a, int n, real_t b[]);
/* CHLSOLBN.C */
void chlsolbnd(real_t a[], int n, int w, real_t b[]);
/* COMEIG1.C */
int comeig1(real_t **a, int n, real_t em[], real_t re[], real_t im[], real_t **vec);
/* COMEIGVA.C */
int comeigval(real_t **a, int n, real_t em[], real_t re[], real_t im[]);
/* COMKWD.C */
void comkwd(real_t pr, real_t pi, real_t qr, real_t qi, real_t *gr, real_t *gi, real_t *kr, real_t *ki);
/* COMVALQR.C */
int comvalqri(real_t **a, int n, real_t em[], real_t re[], real_t im[]);
/* COMVECHE.C */
void comveches(real_t **a, int n, real_t lambda, real_t mu, real_t em[], real_t u[], real_t v[]);
/* CONJGRAD.C */
void conjgrad(void (*matvec)(real_t [], real_t []), real_t x[], real_t r[], int l, int n, int (*goon)(int, real_t), int *iterate, real_t *norm2);
/* DECBND.C */
void decbnd(real_t a[], int n, int lw, int rw, real_t aux[], real_t m[], int p[]);
/* DEC.C */
void dec(real_t **a, int n, real_t aux[], int p[]);
/* DECINV.C */
void decinv(real_t **a, int n, real_t aux[]);
/* DECSOLBN.C */
void decsolbnd(real_t a[], int n, int lw, int rw, real_t aux[], real_t b[]);
/* DECSOL.C */
void decsol(real_t **a, int n, real_t aux[], real_t b[]);
/* DECSOLS2.C */
void decsolsym2(real_t **a, int n, real_t b[], real_t tol, int aux[]);
/* DECSOLST.C */
void decsolsymtri(real_t diag[], real_t co[], int n, real_t aux[], real_t b[]);
/* DECSOLTP.C */
void decsoltripiv(real_t sub[], real_t diag[], real_t super[], int n, real_t aux[], real_t b[]);
/* DECSOLTR.C */
void decsoltri(real_t sub[], real_t diag[], real_t super[], int n, real_t aux[], real_t b[]);
/* DECSYM2.C */
void decsym2(real_t **a, int n, real_t tol, int aux[], int p[], real_t detaux[]);
/* DECSYMTR.C */
void decsymtri(real_t diag[], real_t co[], int n, real_t aux[]);
/* DECTRI.C */
void dectri(real_t sub[], real_t diag[], real_t super[], int n, real_t aux[]);
/* DECTRIPI.C */
void dectripiv(real_t sub[], real_t diag[], real_t super[], int n, real_t aid[], real_t aux[], int piv[]);
/* DETERMBN.C */
real_t determbnd(real_t a[], int n, int lw, int rw, int sgndet);
/* DETERM.C */
real_t determ(real_t **a, int n, int sign);
/* DETMSYM2.C */
real_t determsym2(real_t detaux[], int n, int aux[]);
/* EIGCOM.C */
int eigcom(real_t **ar, real_t **ai, int n, real_t em[], real_t valr[], real_t vali[], real_t **vr, real_t **vi);
/* EIGHRM.C */
void eighrm(real_t **a, int n, int numval, real_t val[], real_t **vecr, real_t **veci, real_t em[]);
/* EIGSYM1.C */
void eigsym1(real_t a[], int n, int numval, real_t val[], real_t **vec, real_t em[]);
/* EIGSYM2.C */
void eigsym2(real_t **a, int n, int numval, real_t val[], real_t **vec, real_t em[]);
/* EIGVALCO.C */
int eigvalcom(real_t **ar, real_t **ai, int n, real_t em[], real_t valr[], real_t vali[]);
/* EIGVALHR.C */
void eigvalhrm(real_t **a, int n, int numval, real_t val[], real_t em[]);
/* EIGVALS1.C */
void eigvalsym1(real_t a[], int n, int numval, real_t val[], real_t em[]);
/* EIGVALS2.C */
void eigvalsym2(real_t **a, int n, int numval, real_t val[], real_t em[]);
/* EQILBR.C */
void eqilbr(real_t **a, int n, real_t em[], real_t d[], int inter[]);
/* EQILBRCO.C */
void eqilbrcom(real_t **a1, real_t **a2, int n, real_t em[], real_t d[], int inter[]);
/* ERBELM.C */
void erbelm(int n, real_t aux[], real_t nrminv);
/* FEMLAG.C */
void femlag(real_t x[], real_t y[], int n, real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FEMLAGSP.C */
void femlagspher(real_t x[], real_t y[], int n, int nc, real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* FEMLAGSY.C */
void femlagsym(real_t x[], real_t y[], int n, real_t (*p)(real_t), real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[]);
/* GSITSOLE.C */
void gssitisolerb(real_t **a, int n, real_t aux[], real_t b[]);
/* GSSELM.C */
void gsselm(real_t **a, int n, real_t aux[], int ri[], int ci[]);
/* GSSERB.C */
void gsserb(real_t **a, int n, real_t aux[], int ri[], int ci[]);
/* GSSINV.C */
void gssinv(real_t **a, int n, real_t aux[]);
/* GSSINVER.C */
void gssinverb(real_t **a, int n, real_t aux[]);
/* GSSITISO.C */
void gssitisol(real_t **a, int n, real_t aux[], real_t b[]);
/* GSSNRI.C */
void gssnri(real_t **a, int n, real_t aux[], int ri[], int ci[]);
/* GSSSOL.C */
void gsssol(real_t **a, int n, real_t aux[], real_t b[]);
/* GSSSOLER.C */
void gsssolerb(real_t **a, int n, real_t aux[], real_t b[]);
/* HESTGL2.C */
void hestgl2(int n, real_t **a, real_t **b);
/* HESTGL3.C */
void hestgl3(int n, real_t **a, real_t **b, real_t **x);
/* HOMSOL.C */
int homsol(real_t **a, int m, int n, real_t **v, real_t em[]);
/* HOMSOLSV.C */
void homsolsvd(real_t **u, real_t val[], real_t **v, int m, int n);
/* HSH2COL.C */
void hsh2col(int la, int lb, int u, int i, real_t a1, real_t a2, real_t **a, real_t **b);
/* HSH2ROW2.C */
void hsh2row2(int l, int ua, int ub, int j, real_t a1, real_t a2, real_t **a, real_t **b);
/* HSH2ROW3.C */
void hsh2row3(int l, int ua, int ub, int ux, int j, real_t a1, real_t a2, real_t **a, real_t **b, real_t **x);
/* HSH3COL.C */
void hsh3col(int la, int lb, int u, int i, real_t a1, real_t a2, real_t a3, real_t **a, real_t **b);
/* HSH3ROW2.C */
void hsh3row2(int l, int u, int j, real_t a1, real_t a2, real_t a3, real_t **a, real_t **b);
/* HSH3ROW3.C */
void hsh3row3(int l, int u, int ux, int j, real_t a1, real_t a2, real_t a3, real_t **a, real_t **b, real_t **x);
/* HSHCOMHE.C */
void hshcomhes(real_t **ar, real_t **ai, int n, real_t em[], real_t b[], real_t tr[], real_t ti[], real_t del[]);
/* HSHDECMU.C */
void hshdecmul(int n, real_t **a, real_t **b, real_t dwarf);
/* HSHHRMTR.C */
void hshhrmtri(real_t **a, int n, real_t d[], real_t b[], real_t bb[], real_t em[], real_t tr[], real_t ti[]);
/* HSHHRMTV.C */
void hshhrmtrival(real_t **a, int n, real_t d[], real_t bb[], real_t em[]);
/* HSHREABI.C */
void hshreabid(real_t **a, int m, int n, real_t d[], real_t b[], real_t em[]);
/* INV1.C */
real_t inv1(real_t **a, int n, int ri[], int ci[], int withnorm);
/* INV.C */
void inv(real_t **a, int n, int p[]);
/* ITISOL.C */
void itisol(real_t **a, real_t **lu, int n, real_t aux[], int ri[], int ci[], real_t b[]);
/* ITISOLER.C */
void itisolerb(real_t **a, real_t **lu, int n, real_t aux[], int ri[], int ci[], real_t b[]);
/* LSQDECOM.C */
void lsqdecomp(real_t **a, int n, int m, int n1, real_t aux[], real_t aid[], int ci[]);
/* LSQDGLIN.C */
void lsqdglinv(real_t **a, int m, real_t aid[], int ci[], real_t diag[]);
/* LSQINV.C */
void lsqinv(real_t **a, int m, real_t aid[], int c[]);
/* LSQORTDE.C */
void lsqortdec(real_t **a, int n, int m, real_t aux[], real_t aid[], int ci[]);
/* LSQORTDS.C */
void lsqortdecsol(real_t **a, int n, int m, real_t aux[], real_t diag[], real_t b[]);
/* LSQREFSO.C */
void lsqrefsol(real_t **a, real_t **qr, int n, int m, int n1, real_t aux[], real_t aid[], int ci[], real_t b[], real_t *ldx, real_t x[], real_t res[]);
/* LSQSOL.C */
void lsqsol(real_t **a, int n, int m, real_t aid[], int ci[], real_t b[]);
/* LUPZEROP.C */
void lupzerortpol(int n, int m, real_t b[], real_t c[], real_t zer[], real_t em[]);
/* MERGESOR.C */
void mergesort(real_t a[], int p[], int low, int up);
void merge(int lo, int ls, int rs, int p[], real_t a[], int hp[]);
/* ONENRMIN.C */
real_t onenrminv(real_t **a, int n);
/* ORTHOG.C */
void orthog(int n, int lc, int uc, real_t **x);
/* PRETFMMA.C */
void pretfmmat(real_t **a, int m, int n, real_t d[]);
/* PSDINV.C */
int psdinv(real_t **a, int m, int n, real_t em[]);
/* PSDINVSV.C */
void psdinvsvd(real_t **u, real_t val[], real_t **v, int m, int n, real_t em[]);
/* PSTTFMMA.C */
void psttfmmat(real_t **a, int n, real_t **v, real_t b[]);
/* QRICOM.C */
int qricom(real_t **a1, real_t **a2, real_t b[], int n, real_t em[], real_t val1[], real_t val2[], real_t **vec1, real_t **vec2);
/* QRIHRM.C */
int qrihrm(real_t **a, int n, real_t val[], real_t **vr, real_t **vi, real_t em[]);
/* QRISNGVA.C */
int qrisngval(real_t **a, int m, int n, real_t val[], real_t em[]);
/* QRISNGVB.C */
int qrisngvalbid(real_t d[], real_t b[], int n, real_t em[]);
/* QRISNGVD.C */
int qrisngvaldec(real_t **a, int m, int n, real_t val[], real_t **v, real_t em[]);
/* QRISNVDB.C */
int qrisngvaldecbid(real_t d[], real_t b[], int m, int n, real_t **u, real_t **v, real_t em[]);
/* QRISYM.C */
int qrisym(real_t **a, int n, real_t val[], real_t em[]);
/* QRISYMTR.C */
int qrisymtri(real_t **a, int n, real_t d[], real_t b[], real_t bb[], real_t em[]);
/* QRIVALHR.C */
int qrivalhrm(real_t **a, int n, real_t val[], real_t em[]);
/* QRIVALS1.C */
int qrivalsym1(real_t a[], int n, real_t val[], real_t em[]);
/* QRIVALS2.C */
int qrivalsym2(real_t **a, int n, real_t val[], real_t em[]);
/* QRIVALST.C */
int qrivalsymtri(real_t d[], real_t bb[], int n, real_t em[]);
/* QZI.C */
void qzi(int n, real_t **a, real_t **b, real_t **x, real_t alfr[], real_t alfi[], real_t beta[], int iter[], real_t em[]);
/* QZIVAL.C */
void qzival(int n, real_t **a, real_t **b, real_t alfr[], real_t alfi[], real_t beta[], int iter[], real_t em[]);
/* REAEIG1.C */
int reaeig1(real_t **a, int n, real_t em[], real_t val[], real_t **vec);
/* REAEIG3.C */
int reaeig3(real_t **a, int n, real_t em[], real_t val[], real_t **vec);
/* REAEIGVA.C */
int reaeigval(real_t **a, int n, real_t em[], real_t val[]);
/* REAQRI.C */
int reaqri(real_t **a, int n, real_t em[], real_t val[], real_t **vec);
/* REAVALQR.C */
int reavalqri(real_t **a, int n, real_t em[], real_t val[]);
/* REAVECHE.C */
void reaveches(real_t **a, int n, real_t lambda, real_t em[], real_t v[]);
/* ROWPERM.C */
void rowperm(int perm[], int low, int upp, int i, real_t **mat);
/* SELZEROP.C */
void selzerortpol(int n, int n1, int n2, real_t b[], real_t c[], real_t zer[], real_t em[]);
/* SOLBND.C */
void solbnd(real_t a[], int n, int lw, int rw, real_t m[], int p[], real_t b[]);
/* SOL.C */
void sol(real_t **a, int n, int p[], real_t b[]);
/* SOLELM.C */
void solelm(real_t **a, int n, int ri[], int ci[], real_t b[]);
/* SOLOVR.C */
int solovr(real_t **a, int m, int n, real_t x[], real_t em[]);
/* SOLSVDOV.C */
void solsvdovr(real_t **u, real_t val[], real_t **v, int m, int n, real_t x[], real_t em[]);
/* SOLSVDUN.C */
void solsvdund(real_t **u, real_t val[], real_t **v, int m, int n, real_t x[], real_t em[]);
/* SOLSYM2.C */
void solsym2(real_t **a, int n, real_t b[], int p[], real_t detaux[]);
/* SOLSYMTR.C */
void solsymtri(real_t diag[], real_t co[], int n, real_t b[]);
/* SOLTRI.C */
void soltri(real_t sub[], real_t diag[], real_t super[], int n, real_t b[]);
/* SOLTRIPI.C */
void soltripiv(real_t sub[], real_t diag[], real_t super[], int n, real_t aid[], int piv[], real_t b[]);
/* SOLUND.C */
int solund(real_t **a, int m, int n, real_t x[], real_t em[]);
/* SYMEIGIM.C */
void symeigimp(int n, real_t **a, real_t **vec, real_t val[], real_t lbound[], real_t ubound[], real_t aux[]);
/* TFMPREVE.C */
void tfmprevec(real_t **a, int n);
/* TFMREAHE.C */
void tfmreahes(real_t **a, int n, real_t em[], int index[]);
/* TFMSYMT1.C */
void tfmsymtri1(real_t a[], int n, real_t d[], real_t b[], real_t bb[], real_t em[]);
/* TFMSYMT2.C */
void tfmsymtri2(real_t **a, int n, real_t d[], real_t b[], real_t bb[], real_t em[]);
/* VALQRICO.C */
int valqricom(real_t **a1, real_t **a2, real_t b[], int n, real_t em[], real_t val1[], real_t val2[]);
/* VALSYMTR.C */
void valsymtri(real_t d[], real_t bb[], int n, int n1, int n2, real_t val[], real_t em[]);
real_t sturm(real_t d[], real_t bb[], int n, real_t x, int k, real_t machtol, real_t max, int *count, real_t *lb, real_t *ub);
/* VECPERM.C */
void vecperm(int perm[], int low, int upp, real_t vector[]);
/* VECSYMTR.C */
void vecsymtri(real_t d[], real_t b[], int n, int n1, int n2, real_t val[], real_t **vec, real_t em[]);
/* ZERPOL.C */
int zerpol(int n, real_t a[], real_t em[], real_t re[], real_t im[], real_t d[]);
int zerpolfunction(int n, real_t d[], real_t f[], real_t x, real_t y, real_t tol, int *it, real_t *newf);
