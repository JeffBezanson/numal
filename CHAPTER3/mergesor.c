#include "../real.h"
void mergesort(real_t a[], int p[], int low, int up)
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void merge(int, int, int, int [], real_t [], int []);
	int i,lo,step,stap,umlp1,umsp1,rest,restv,*hp;

	hp=allocate_integer_vector(low,up);
	for (i=low; i<=up; i++) p[i]=i;
	restv=0;
	umlp1=up-low+1;
	step=1;
	do {
		stap=2*step;
		umsp1=up-stap+1;
		for (lo=low; lo<=umsp1; lo += stap)
			merge(lo,step,step,p,a,hp);
		rest=up-lo+1;
		if (rest > restv && restv > 0)
			merge(lo,rest-restv,restv,p,a,hp);
		restv=rest;
		step *= 2;
	} while (step < umlp1);
	free_integer_vector(hp,low);
}

void merge(int lo, int ls, int rs, int p[], real_t a[], int hp[])
{
	/* this procedure is used internally by MERGESORT */

	int l,r,lout,rout,i,pl,pr;

	l=lo;
	r=lo+ls;
	lout=rout=0;
	i=lo;
	do {
		pl=p[l];
		pr=p[r];
		if (a[pl] > a[pr]) {
			hp[i]=pr;
			r++;
			rout = (r == lo+ls+rs);
		} else {
			hp[i]=pl;
			l++;
			lout = (l == lo+ls);
		}
		i++;
	} while (!(lout || rout));
	if (rout) {
		for (i=lo+ls-1; i>=l; i--) p[i+rs]=p[i];
		r=l+rs;
	}
	for (i=r-1; i>=lo; i--) p[i]=hp[i];
}
