#include "../real.h"


real_t sumposseries(real_t (*ai)(real_t), int maxaddup, real_t maxzero,
						int maxrecurs, int machexp, int tim)
{
	real_t sumposseriessumup(int, real_t (*)(real_t), int, real_t, int,
				int, int, int, int, int, int, int);
	int recurs,vl,vl2,vl4;

	recurs=0;
	vl=1000;
	vl2=2*vl;
	vl4=2*vl2;
	return sumposseriessumup(0,ai,maxaddup,maxzero,maxrecurs,
					machexp,tim,recurs,vl,vl2,vl4,0);
}

real_t sumposseriessumup(int bjk, real_t (*ai)(real_t),
				int maxaddup, real_t maxzero, int maxrecurs,
				int machexp, int tim, int recurs, int vl,
				int vl2, int vl4, int jj)
{
	/* this function is internally used by SUMPOSSERIES */

	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t sumposseriesbjk(int, real_t, real_t, real_t (*)(real_t));
	int j,transform,j2,jodd,k,n,t;
	real_t i,sum,nextterm,*v,mn,mp,ds,esum,m[16],temp,vj;

	i=maxaddup+1;
	j=1;
	transform=0;
	while (1) {
		temp = (bjk) ? sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
		if (temp <= maxzero) {
			if (j >= tim) break;
			j++;
			i++;
		} else {
			if (recurs != maxrecurs) transform=1;
			break;
		}
	}
	if (!transform) {
		sum=i=0.0;
		j=0;
		do {
			i++;
			nextterm = (bjk) ?
				sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
			j = (nextterm <= maxzero) ? j+1 : 0;
			sum += nextterm;
		} while (j < tim);
		return sum;
	}
	/* transform series */
	v=allocate_real_vector(1,vl);
	j2=0;
	jodd=1;
	/* euler */
	n=t=j=0;
	jj=j+1;
	if (jodd) {
		jodd=0;
		recurs++;
		temp=vj=sumposseriessumup(1,ai,maxaddup,maxzero,
					maxrecurs,machexp,tim,recurs,vl,vl2,vl4,jj);
		recurs--;
		if (jj <= vl)
			v[jj]=temp;
		else
			if (jj <= vl2) v[jj-vl]=temp;
	} else {
		jodd=1;
		if (jj > vl4) {
			recurs++;
			vj = -sumposseriessumup(1,ai,maxaddup,maxzero,
						maxrecurs,machexp,tim,recurs,vl,vl2,vl4,jj);
			recurs--;
		} else {
			j2++;
			i=j2;
			if (jj > vl2) {
				temp = (bjk) ?
							sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
				vj = -(v[j2-vl]-temp)/2.0;
			}
			else {
				temp = (bjk) ?
							sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
				temp=v[(jj <= vl) ? jj : jj-vl]=(v[j2]-temp)/2.0;
				vj = -temp;
			}
		}
	}
	m[0]=vj;
	esum=m[0]/2.0;
	do {
		j++;
		jj=j+1;
		if (jodd) {
			jodd=0;
			recurs++;
			temp=vj=sumposseriessumup(1,ai,maxaddup,maxzero,
							maxrecurs,machexp,tim,recurs,vl,vl2,vl4,jj);
			recurs--;
			if (jj <= vl)
				v[jj]=temp;
			else
				if (jj <= vl2) v[jj-vl]=temp;
		} else {
			jodd=1;
			if (jj > vl4) {
				recurs++;
				vj = -sumposseriessumup(1,ai,maxaddup,maxzero,
							maxrecurs,machexp,tim,recurs,vl,vl2,vl4,jj);
				recurs--;
			} else {
				j2++;
				i=j2;
				if (jj > vl2) {
					temp = (bjk) ?
								sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
					vj = -(v[j2-vl]-temp)/2.0;
				} else {
					temp = (bjk) ?
								sumposseriesbjk(jj,i,machexp,ai) : (*ai)(i);
					temp=v[(jj <= vl) ? jj : jj-vl]=(v[j2]-temp)/2.0;
					vj = -temp;
				}
			}
		}
		mn=vj;
		for (k=0; k<=n; k++) {
			mp=(mn+m[k])/2.0;
			m[k]=mn;
			mn=mp;
		}
		if (fabs(mn) < fabs(m[n]) && n < 15) {
			ds=mn/2.0;
			n++;
			m[n]=mn;
		} else
			ds=mn;
		esum += ds;
		t = (fabs(ds) < maxzero) ? t+1 : 0;
	} while (t < tim);
	free_real_vector(v,1);
	return esum;
}

real_t sumposseriesbjk(int j, real_t i, real_t machexp,
					real_t (*ai)(real_t))
{
	/* this function is internally used by SUMPOSSERIES */

	real_t coeff;

	if (i > machexp) return 0.0;
	coeff=pow(2.0,i-1.0);
	return coeff*(*ai)(j*coeff);
}

