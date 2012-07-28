#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void lsqortdec(real_t **, int, int, real_t [], real_t [], int []);
	void lsqinv(real_t **, int, real_t [], int []);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	int i,j,piv[3];
	real_t **a,**c,**t,aid[3],aux[6];

	a=allocate_real_matrix(1,5,1,2);
	c=allocate_real_matrix(1,5,1,2);
	t=allocate_real_matrix(1,2,1,2);
	aux[2]=1.0e-6;
	a[1][1]=c[1][1] = -2.0;
	a[2][1]=c[2][1] = -1.0;
	a[3][1]=c[3][1]=1.0;
	a[4][1]=c[4][1]=2.0;
	a[5][1]=c[5][1]=1.0;
	a[1][2]=a[2][2]=a[3][2]=a[4][2]=c[1][2]=c[2][2]=c[3][2]=c[4][2]=1.0;
	a[5][2]=c[5][2]=2.0;
	lsqortdec(a,5,2,aux,aid,piv);
	if (aux[3] == 2) {
		lsqinv(a,2,aid,piv);
		t[1][1]=a[1][1];
		t[2][2]=a[2][2];
		t[2][1]=t[1][2]=a[1][2];
		for (j=1; j<=2; j++)
			for (i=1; i<=5; i++) a[i][j]=matmat(1,2,i,j,c,t);
		printf("Aux[2, 3, 5] =  %e   %e   %e\n"
				"\nInverse:\n   %e   %e\n  %e    %e\n"
				"\nCheck:  S' * (S * T) :\n   %e   %e\n  %e   %e\n",
				aux[2],aux[3],aux[5],t[1][1],t[1][2],t[2][1],t[2][2],
				tammat(1,5,1,1,c,a),tammat(1,5,1,2,c,a),
				tammat(1,5,2,1,c,a),tammat(1,5,2,2,c,a));
	}
	free_real_matrix(a,1,5,1);
	free_real_matrix(c,1,5,1);
	free_real_matrix(t,1,2,1);
}

