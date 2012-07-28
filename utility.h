/* FREE_IM.C */
void free_integer_matrix(int **m, int lr, int ur, int lc);
/* FREE_IV.C */
void free_integer_vector(int *v, int l);
/* FREE_RM.C */
void free_real_matrix(real_t **m, int lr, int ur, int lc);
/* FREE_RV.C */
void free_real_vector(real_t *v, int l);
/* IMATRIX.C */
int **allocate_integer_matrix(int lr, int ur, int lc, int uc);
/* IVECTOR.C */
int *allocate_integer_vector(int l, int u);
/* RMATRIX.C */
real_t **allocate_real_matrix(int lr, int ur, int lc, int uc);
/* RVECTOR.C */
real_t *allocate_real_vector(int l, int u);
