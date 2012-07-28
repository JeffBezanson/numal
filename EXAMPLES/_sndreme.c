#include <stdio.h>

void main ()
{
	int s[3];
	real_t em[2],g[8];

	printf("SNDREMEZ delivers:\n\n");
	g[0]=10.0;  g[1]=12.0;  g[2] = -15.0;  g[3] = -10.0;
	g[4] = -14.0;  g[5]=15.0;  g[6]=10.0;  g[7]=11.0;
	em[0]=10.0;  s[0]=0.0;  s[1]=3.0;  s[2]=6.0;
	printf("The numbers:\n S[J]:  %d  %d  %d\n G[S[j]]:"
			"  %4.0f %4.0f %4.0f\n\n",
			s[0],s[1],s[2],g[s[0]],g[s[1]],g[s[2]]);
	sndremez(2,7,s,g,em);
	printf("are exchanged with:\n S[J]:"
			"  %d  %d  %d\n G[S[j]]:  %4.0f %4.0f %4.0f\n\n"
			"The reference set of function values is:\n"
			" %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n",
			s[0],s[1],s[2],g[s[0]],g[s[1]],g[s[2]],
			g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
}

