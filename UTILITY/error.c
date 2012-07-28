#include "../real.h"
#include <stdlib.h>
#include <stdio.h>

void system_error(char error_message[])
{
	void exit(int);

	printf("%s",error_message);
	exit(1);
}
