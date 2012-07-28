/* ALLCHEPO.C */
void allchepol(int n, real_t x, real_t t[]);
/* ALLORTPO.C */
void allortpol(int n, real_t x, real_t b[], real_t c[], real_t p[]);
/* ALLORTPS.C */
void allortpolsym(int n, real_t x, real_t c[], real_t p[]);
/* CHEPOL.C */
real_t chepol(int n, real_t x);
/* CHEPOLSU.C */
real_t chepolsum(int n, real_t x, real_t a[]);
/* CHSPOL.C */
void chspol(int n, real_t a[]);
/* COMFOUS1.C */
void comfouser1(int n, real_t theta, real_t ar[], real_t ai[], real_t *rr, real_t *ri);
/* COMFOUS2.C */
void comfouser2(int n, real_t theta, real_t ar[], real_t ai[], real_t *rr, real_t *ri);
/* COMFOUSE.C */
void comfouser(int n, real_t theta, real_t a[], real_t *rr, real_t *ri);
/* COSSER.C */
real_t cosser(int n, real_t theta, real_t a[]);
/* DERPOL.C */
void derpol(int n, int k, real_t x, real_t a[]);
/* FOUSER1.C */
real_t fouser1(int n, real_t theta, real_t a[], real_t b[]);
/* FOUSER2.C */
real_t fouser2(int n, real_t theta, real_t a[], real_t b[]);
/* FOUSER.C */
real_t fouser(int n, real_t theta, real_t a[]);
/* GRNNEW.C */
void grnnew(int n, real_t x[], real_t a[]);
/* INTCHS.C */
void intchs(int n, real_t a[], real_t b[]);
/* JFRAC.C */
real_t jfrac(int n, real_t a[], real_t b[]);
/* LINTFMPO.C */
void lintfmpol(real_t p, real_t q, int n, real_t a[]);
/* NEWGRN.C */
void newgrn(int n, real_t x[], real_t a[]);
/* NORDERPO.C */
void norderpol(int n, int k, real_t x, real_t a[]);
/* ODDCHEPO.C */
real_t oddchepolsum(int n, real_t x, real_t a[]);
/* ORTPOL.C */
real_t ortpol(int n, real_t x, real_t b[], real_t c[]);
/* ORTPOLSY.C */
real_t ortpolsym(int n, real_t x, real_t c[]);
/* POL.C */
real_t pol(int n, real_t x, real_t a[]);
/* POLCHS.C */
void polchs(int n, real_t a[]);
/* POLSHTCH.C */
void polshtchs(int n, real_t a[]);
/* SHTCHSPO.C */
void shtchspol(int n, real_t a[]);
/* SINSER.C */
real_t sinser(int n, real_t theta, real_t b[]);
/* SUMORTPO.C */
real_t sumortpol(int n, real_t x, real_t b[], real_t c[], real_t a[]);
/* SUMORTPS.C */
real_t sumortpolsym(int n, real_t x, real_t c[], real_t a[]);
/* TAYPOL.C */
void taypol(int n, int k, real_t x, real_t a[]);
