#include <stdio.h>
void main ()
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void lsqortdec(real_t **, int, int, real_t [], real_t [], int []);
	void lsqsol(real_t **, int, int, real_t [], int [], real_t []);
	void lsqdglinv(real_t **, int, real_t [], int [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	int i,piv[3];
	real_t sum,temp,**a,**c,b[6],x[6],diag[3],aid[3],aux[6];

	a=allocate_real_matrix(1,5,1,2);
	c=allocate_real_matrix(1,5,1,2);
	aux[2]=1.0e-6;
	a[1][1]=c[1][1] = -2.0;
	a[2][1]=c[2][1] = -1.0;
	a[3][1]=c[3][1]=1.0;
	a[4][1]=c[4][1]=2.0;
	a[5][1]=c[5][1]=1.0;
	a[1][2]=a[2][2]=a[3][2]=a[4][2]=c[1][2]=c[2][2]=c[3][2]=c[4][2]=1.0;
	a[5][2]=c[5][2]=2.0;
	b[1]=x[1]=0.0;
	b[2]=x[2]=1.0;
	b[3]=x[3]=b[4]=x[4]=2.0;
	b[5]=x[5]=3.0;
	lsqortdec(a,5,2,aux,aid,piv);
	if (aux[3] == 2) {
		lsqsol(a,5,2,aid,piv,x);
		lsqdglinv(a,2,aid,piv,diag);
		sum=0.0;
		for (i=1; i<=5; i++) {
			temp=b[i]-c[i][1]*x[1]-c[i][2]*x[2];
			sum += temp*temp;
		}
		printf("Aux[2, 3, 5] =  %e   %e   %e\n"
				"LSQ solution :  %e   %e\n"
				"Residue (delivered) :  %e\n"
				"Residue (checked)   :  %e\n"
				"Diagonal of inverse M'M :  %e   %e\n",
				aux[2],aux[3],aux[5],x[1],x[2],sqrt(vecvec(3,5,0,x,x)),
				sqrt(sum),diag[1],diag[2]);
	}
	free_real_matrix(a,1,5,1);
	free_real_matrix(c,1,5,1);
}

