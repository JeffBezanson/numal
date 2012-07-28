#include <stdio.h>
void main ()
{
	void bessj(real_t, int, real_t []);
	real_t bessj0(real_t);
	real_t bessj1(real_t);
	int i;
	real_t j[2];
	real_t x[4]={1.0, 5.0, 10.0, 25.0};

	printf("BESSJ0 delivers:\n");
	for (i=0; i<=3; i++) {
		bessj(x[i],1,j);
		printf(" %6.1f   %+7.2e   %+7.2e\n",
				x[i],j[0]-bessj0(x[i]),j[1]-bessj1(x[i]));
	}
}

