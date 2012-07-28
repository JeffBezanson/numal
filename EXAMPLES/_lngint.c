#include <stdio.h>
void main ()
{
	void lngintadd(int [], int [], int []);
	void lngintsubtract(int [], int [], int []);
	void lngintmult(int [], int [], int []);
	void lngintdivide(int [], int [], int [], int []);
	void lngintpower(int [], int, int []);
	int i,u[100],v[100],r1[100],r2[100];

	u[0]=5;
	u[1]=33;
	u[2]=u[3]=u[4]=u[5]=70;
	v[0]=2;
	v[1]=4;
	v[2]=44;
	printf("\nInput numbers:\n");
	for (i=1; i<=u[0]; i++) printf("%4d",u[i]);
	printf("\n");
	for (i=1; i<=v[0]; i++) printf("%4d",v[i]);
	printf("\n\nadd: ");

	lngintadd(u,v,r1);
	for (i=1; i<=r1[0]; i++) printf("%4d",r1[i]);
	printf("\nsubtract: ");
   lngintsubtract(u,v,r1);
	for (i=1; i<=r1[0]; i++) printf("%4d",r1[i]);
	printf("\nmultiple: ");
	lngintmult(u,v,r1);
	for (i=1; i<=r1[0]; i++) printf("%4d",r1[i]);
	printf("\ndivide:  Quotient: ");
	lngintdivide(u,v,r1,r2);
	for (i=1; i<=r1[0]; i++) printf("%4d",r1[i]);
	printf("  Remainder: ");
	for (i=1; i<=r2[0]; i++) printf("%4d",r2[i]);
	printf("\npower: ");
	lngintpower(v,4,r1);
	for (i=1; i<=r1[0]; i++) printf("%4d",r1[i]);

	/*  solution : sum:          33 70 70 75 14
						subtract:     33 70 70 66 26
						mult:    1 49 65 93 93 90 80
						div:   quotient:  7 59 16 82  remainder: 2 62
						power:      3 88 62 60 24 96
	*/
}

