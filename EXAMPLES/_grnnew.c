#include <stdio.h>
void main ()
{
	void grnnew(int, real_t [], real_t []);
	void newgrn(int, real_t [], real_t []);
	real_t x[2] = {1.0, 2.0},  a[3] = {1.0, 2.0, 3.0};

	printf("           a[0]   a[1]   a[2]\ninput:"
			"  %7.2f%7.2f%7.2f",a[0],a[1],a[2]);
	grnnew(2,x,a);
	printf("\ngrnnew: %7.2f%7.2f%7.2f",a[0],a[1],a[2]);
	newgrn(2,x,a);
	printf("\nnewgrn: %7.2f%7.2f%7.2f",a[0],a[1],a[2]);
}

