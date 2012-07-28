#include <stdio.h>
void main ()
{
	void allchepol(int, real_t, real_t []);
	real_t chepol(int, real_t);
	real_t t[3];

	allchepol(2,-1.0,t);
	printf("Delivers:\n %4.0f%4.0f%4.0f",t[0],t[1],t[2]);
	allchepol(2,0.0,t);
	printf("\n %4.0f%4.0f%4.0f",t[0],t[1],t[2]);
	allchepol(2,1.0,t);
	printf("\n %4.0f%4.0f%4.0f",t[0],t[1],t[2]);
	printf("\n\n %-6.0f%-6.0f%-6.0f",
			chepol(2,-1.0),chepol(2,0.0),chepol(2,1.0));
}

