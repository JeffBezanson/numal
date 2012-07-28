#include <stdio.h>

void main ()
{
	void ini(int, int, int []);
	int s[3];

	ini(2,20,s);
	printf("INI selects out of 0,1,...,20 the numbers:\n"
			"  %4d %4d %4d\n",s[0],s[1],s[2]);
}

