#include <stdio.h>
void main ()
{
	real_t chepolsum(int, real_t, real_t []);
	real_t a[3] = {1.0, 0.5, 0.25};

	printf("CHEPOLSUM delivers:  %-6.2f%-6.2f%-6.2f",
			chepolsum(2,-1.0,a),chepolsum(2,0.0,a),chepolsum(2,1.0,a));
}

