#include <stdio.h>
void main ()
{
	void polshtchs(int, real_t []);
	void shtchspol(int, real_t []);
	real_t a[3];

	a[0]=1.0;	a[1]=2.0;	a[2]=3.0;
	printf("            a[0]  a[1]  a[2]\ninput:"
			"      %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
	polshtchs(2,a);
	printf("\npolshtchs:  %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
	shtchspol(2,a);
	printf("\nshtchspol:  %-6.2f%-6.2f%-6.2f",a[0],a[1],a[2]);
}

