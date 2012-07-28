#include <stdio.h>
void main ()
{
	void polchs(int, real_t []);
	void chspol(int, real_t []);
	real_t a[3] = {1.0, 2.0, 3.0};

	printf("         a[0]  a[1]  a[2]\ninput:"
			"   %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
	polchs(2,a);
	printf("\npolchs:  %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
	chspol(2,a);
	printf("\nchspol:  %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
}

