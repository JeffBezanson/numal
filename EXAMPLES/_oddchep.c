#include <stdio.h>
void main ()
{
	real_t oddchepolsum(int, real_t, real_t []);
	real_t a[2] = {0.5, 0.2};

	printf("ODDCHEPOLSUM delivers:  %-7.2f%-7.2f%-7.2f",
			oddchepolsum(1,-1.0,a),oddchepolsum(1,0.0,a),
			oddchepolsum(1,1.0,a));
}

