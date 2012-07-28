/* ABSMAXMA.C */
real_t absmaxmat(int lr, int ur, int lc, int uc, int *i, int *j, real_t **a);
/* CARPOL.C */
void carpol(real_t ar, real_t ai, real_t *r, real_t *c, real_t *s);
/* CHSH2.C */
void chsh2(real_t a1r, real_t a1i, real_t a2r, real_t a2i, real_t *c, real_t *sr, real_t *si);
/* COLCST.C */
void colcst(int l, int u, int j, real_t **a, real_t x);
/* COMABS.C */
real_t comabs(real_t xr, real_t xi);
/* COMCOLCS.C */
void comcolcst(int l, int u, int j, real_t **ar, real_t **ai, real_t xr, real_t xi);
/* COMDIV.C */
void comdiv(real_t xr, real_t xi, real_t yr, real_t yi, real_t *zr, real_t *zi);
/* COMEUCNR.C */
real_t comeucnrm(real_t **ar, real_t **ai, int lw, int n);
/* COMMATVE.C */
void commatvec(int l, int u, int i, real_t **ar, real_t **ai, real_t br[], real_t bi[], real_t *rr, real_t *ri);
/* COMMUL.C */
void commul(real_t ar, real_t ai, real_t br, real_t bi, real_t *rr, real_t *ri);
/* COMROWCS.C */
void comrowcst(int l, int u, int i, real_t **ar, real_t **ai, real_t xr, real_t xi);
/* COMSCL.C */
void comscl(real_t **a, int n, int n1, int n2, real_t im[]);
/* COMSQRT.C */
void comsqrt(real_t ar, real_t ai, real_t *pr, real_t *pi);
/* DUPCOLVE.C */
void dupcolvec(int l, int u, int j, real_t **a, real_t b[]);
/* DUPMAT.C */
void dupmat(int l, int u, int i, int j, real_t **a, real_t **b);
/* DUPROWVE.C */
void duprowvec(int l, int u, int i, real_t **a, real_t b[]);
/* DUPVEC.C */
void dupvec(int l, int u, int shift, real_t a[], real_t b[]);
/* DUPVECCO.C */
void dupveccol(int l, int u, int j, real_t a[], real_t **b);
/* DUPVECRO.C */
void dupvecrow(int l, int u, int i, real_t a[], real_t **b);
/* ELMCOL.C */
void elmcol(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* ELMCOLRO.C */
void elmcolrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* ELMCOLVE.C */
void elmcolvec(int l, int u, int i, real_t **a, real_t b[], real_t x);
/* ELMCOMCO.C */
void elmcomcol(int l, int u, int i, int j, real_t **ar, real_t **ai, real_t **br, real_t **bi, real_t xr, real_t xi);
/* ELMCOMRO.C */
void elmcomrowvec(int l, int u, int i, real_t **ar, real_t **ai, real_t br[], real_t bi[], real_t xr, real_t xi);
/* ELMCOMVE.C */
void elmcomveccol(int l, int u, int j, real_t ar[], real_t ai[], real_t **br, real_t **bi, real_t xr, real_t xi);
/* ELMROW.C */
void elmrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* ELMROWCO.C */
void elmrowcol(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* ELMROWVE.C */
void elmrowvec(int l, int u, int i, real_t **a, real_t b[], real_t x);
/* ELMVEC.C */
void elmvec(int l, int u, int shift, real_t a[], real_t b[], real_t x);
/* ELMVECCO.C */
void elmveccol(int l, int u, int i, real_t a[], real_t **b, real_t x);
/* ELMVECRO.C */
void elmvecrow(int l, int u, int i, real_t a[], real_t **b, real_t x);
/* FULMATVE.C */
void fulmatvec(int lr, int ur, int lc, int uc, real_t **a, real_t b[], real_t c[]);
/* FULSYMMA.C */
void fulsymmatvec(int lr, int ur, int lc, int uc, real_t a[], real_t b[], real_t c[]);
/* FULTAMVE.C */
void fultamvec(int lr, int ur, int lc, int uc, real_t **a, real_t b[], real_t c[]);
/* HSHCOLMA.C */
void hshcolmat(int lr, int ur, int lc, int uc, int i, real_t x, real_t **u, real_t **a);
/* HSHCOLTA.C */
void hshcoltam(int lr, int ur, int lc, int uc, int i, real_t x, real_t **u, real_t **a);
/* HSHCOMCO.C */
int hshcomcol(int l, int u, int j, real_t **ar, real_t **ai, real_t tol, real_t *k, real_t *c, real_t *s, real_t *t);
/* HSHCOMPR.C */
void hshcomprd(int i, int ii, int l, int u, int j, real_t **ar, real_t **ai, real_t **br, real_t **bi, real_t t);
/* HSHROWMA.C */
void hshrowmat(int lr, int ur, int lc, int uc, int i, real_t x, real_t **u, real_t **a);
/* HSHROWTA.C */
void hshrowtam(int lr, int ur, int lc, int uc, int i, real_t x, real_t **u, real_t **a);
/* HSHVECMA.C */
void hshvecmat(int lr, int ur, int lc, int uc, real_t x, real_t u[], real_t **a);
/* HSHVECTA.C */
void hshvectam(int lr, int ur, int lc, int uc, real_t x, real_t u[], real_t **a);
/* ICHCOL.C */
void ichcol(int l, int u, int i, int j, real_t **a);
/* ICHROW.C */
void ichrow(int l, int u, int i, int j, real_t **a);
/* ICHROWCO.C */
void ichrowcol(int l, int u, int i, int j, real_t **a);
/* ICHSEQ.C */
void ichseq(int l, int u, int il, int shift, real_t a[]);
/* ICHSEQVE.C */
void ichseqvec(int l, int u, int il, int shift, real_t a[]);
/* ICHVEC.C */
void ichvec(int l, int u, int shift, real_t a[]);
/* INFNRMCO.C */
real_t infnrmcol(int l, int u, int j, int *k, real_t **a);
/* INFNRMMA.C */
real_t infnrmmat(int lr, int ur, int lc, int uc, int *kr, real_t **a);
/* INFNRMRO.C */
real_t infnrmrow(int l, int u, int i, int *k, real_t **a);
/* INFNRMVE.C */
real_t infnrmvec(int l, int u, int *k, real_t a[]);
/* INIMAT.C */
void inimat(int lr, int ur, int lc, int uc, real_t **a, real_t x);
/* INIMATD.C */
void inimatd(int lr, int ur, int shift, real_t **a, real_t x);
/* INISYMD.C */
void inisymd(int lr, int ur, int shift, real_t a[], real_t x);
/* INISYMRO.C */
void inisymrow(int l, int u, int i, real_t a[], real_t x);
/* INIVEC.C */
void inivec(int l, int u, real_t a[], real_t x);
/* LNGINTAD.C */
void lngintadd(int u[], int v[], int sum[]);
/* LNGINTDI.C */
void lngintdivide(int u[], int v[], int quotient[], int remainder[]);
/* LNGINTMU.C */
void lngintmult(int u[], int v[], int product[]);
/* LNGINTPO.C */
void lngintpower(int u[], int exponent, int result[]);
/* LNGINTSU.C */
void lngintsubtract(int u[], int v[], int difference[]);
/* MATMAT.C */
real_t matmat(int l, int u, int i, int j, real_t **a, real_t **b);
/* MATTAM.C */
real_t mattam(int l, int u, int i, int j, real_t **a, real_t **b);
/* MATVEC.C */
real_t matvec(int l, int u, int i, real_t **a, real_t b[]);
/* MAXELMRO.C */
int maxelmrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* MULCOL.C */
void mulcol(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* MULROW.C */
void mulrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x);
/* MULVEC.C */
void mulvec(int l, int u, int shift, real_t a[], real_t b[], real_t x);
/* ONENRMCO.C */
real_t onenrmcol(int l, int u, int j, real_t **a);
/* ONENRMMA.C */
real_t onenrmmat(int lr, int ur, int lc, int uc, int *kc, real_t **a);
/* ONENRMRO.C */
real_t onenrmrow(int l, int u, int i, real_t **a);
/* ONENRMVE.C */
real_t onenrmvec(int l, int u, real_t a[]);
/* REASCL.C */
void reascl(real_t **a, int n, int n1, int n2);
/* RESVEC.C */
void resvec(int lr, int ur, int lc, int uc, real_t **a, real_t b[], real_t c[], real_t x);
/* ROTCOL.C */
void rotcol(int l, int u, int i, int j, real_t **a, real_t c, real_t s);
/* ROTCOMCO.C */
void rotcomcol(int l, int u, int i, int j, real_t **ar, real_t **ai, real_t cr, real_t ci, real_t s);
/* ROTCOMRO.C */
void rotcomrow(int l, int u, int i, int j, real_t **ar, real_t **ai, real_t cr, real_t ci, real_t s);
/* ROTROW.C */
void rotrow(int l, int u, int i, int j, real_t **a, real_t c, real_t s);
/* ROWCST.C */
void rowcst(int l, int u, int i, real_t **a, real_t x);
/* SCAPRD1.C */
real_t scaprd1(int la, int sa, int lb, int sb, int n, real_t a[], real_t b[]);
/* SCLCOM.C */
void sclcom(real_t **ar, real_t **ai, int n, int n1, int n2);
/* SEQVEC.C */
real_t seqvec(int l, int u, int il, int shift, real_t a[], real_t b[]);
/* SYMMATVE.C */
real_t symmatvec(int l, int u, int i, real_t a[], real_t b[]);
/* SYMRESVE.C */
void symresvec(int lr, int ur, int lc, int uc, real_t a[], real_t b[], real_t c[], real_t x);
/* TAMMAT.C */
real_t tammat(int l, int u, int i, int j, real_t **a, real_t **b);
/* TAMVEC.C */
real_t tamvec(int l, int u, int i, real_t **a, real_t b[]);
/* VECVEC.C */
real_t vecvec(int l, int u, int shift, real_t a[], real_t b[]);
