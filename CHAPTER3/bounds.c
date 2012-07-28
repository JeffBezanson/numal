#include "../real.h"

#include "../utility.h"

void bounds(int n, real_t a[], real_t re[], real_t im[], real_t rele,
				real_t abse, real_t recentre[], real_t imcentre[],
				real_t bound[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void kcluster(int, int, int, real_t [], real_t [], real_t [],
						real_t [], real_t [], real_t [], real_t []);
	int i,j,k,index1,index2,goon,place,clustin;
	real_t h,min,recent,imcent,xk,yk,zk,corr,*rc,*c,*rce,*clust,
			boundin,*wa1,*wa2,temp1,temp2;

	rc=allocate_real_vector(0,n);
	c=allocate_real_vector(0,n);
	rce=allocate_real_vector(0,n);
	clust=allocate_real_vector(1,n);
	rc[0]=c[0]=a[n];
	rce[0]=fabs(c[0]);
	k=0;
	for (i=1; i<=n; i++) {
		rc[i]=rce[i]=0.0;
		c[i]=a[n-i];
	}
	while (k < n) {
		k++;
		xk=re[k];
		yk=im[k];
		zk=xk*xk+yk*yk;
		for (j=k; j>=1; j--) rce[j] += rce[j-1]*sqrt(zk);
		if (yk == 0.0)
			for (j=k; j>=1; j--) rc[j] -= xk*rc[j-1];
		else {
			k++;
			if (k <= n && xk == re[k] && yk == -im[k]) {
				xk = -2.0*xk;
				for (j=k; j>=1; j--) rce[j] += rce[j-1]*sqrt(zk);
				for (j=k; j>=2; j--) rc[j] += xk*rc[j-1]+zk*rc[j-2];
				rc[1] += xk*rc[0];
			}
		}
	}
	rc[0]=rce[0];
	corr=1.06*FLT_MIN;
	for (i=1; i<=n-1; i++)
		rc[i]=fabs(rc[i]-c[i])+rce[i]*corr*(n+i-2)+
				rele*fabs(c[i])+abse;
	rc[n]=fabs(rc[n]-c[n])+rce[n]*corr*(n-1)+rele*fabs(c[n])+abse;
	for (i=1; i<=n; i++)
		kcluster(1,i,n,rc,re,im,recentre,imcentre,bound,clust);
	goon=1;
	while (goon) {
		index1=index2=0;
		min=FLT_MAX;
		i=n-clust[n]+1;
		while (i >= 2) {
			j=i;
			recent=recentre[i];
			imcent=imcentre[i];
			while (j >= 2) {
				j -= clust[j-1];
				temp1=recent-recentre[j];
				temp2=imcent-imcentre[j];
				h=sqrt(temp1*temp1+temp2*temp2);
				if (h < bound[i]+bound[j] && h <= min) {
					index1=j;
					index2=i;
					min=h;
				}
			}
			i -= clust[i-1];
		}
		if (index1 == 0)
			goon=0;
		else {
			if (imcentre[index1] == 0.0) {
				if (imcentre[index2] != 0.0) clust[index2] *= 2.0;
			}
			else
				if (imcentre[index2] == 0.0) clust[index1] *= 2.0;
			k=index1+clust[index1];
			if (k != index2) {
				/*	shift */
				wa1=allocate_real_vector(1,clust[index2]);
				wa2=allocate_real_vector(1,clust[index2]);
				clustin=clust[index2];
				boundin=bound[index2];
				imcent=imcentre[index2];
				recent=recentre[index2];
				for (j=1; j<=clustin; j++) {
					place=index2+j-1;
					wa1[j]=re[place];
					wa2[j]=im[place];
				}
				for (j=index2-1; j>=k; j--) {
					place=j+clustin;
					re[place]=re[j];
					im[place]=im[j];
					clust[place]=clust[j];
					bound[place]=bound[j];
					recentre[place]=recentre[j];
					imcentre[place]=imcentre[j];
				}
				for (j=k+clustin-1; j>=k; j--) {
					place=j+1-k;
					re[j]=wa1[place];
					im[j]=wa2[place];
					bound[j]=boundin;
					clust[j]=clustin;
					recentre[j]=recent;
					imcentre[j]=imcent;
				}
				free_real_vector(wa1,1);
				free_real_vector(wa2,1);
			}	/* end of shift */
			k=clust[index1]+clust[k];
			kcluster(k,index1,n,rc,re,im,recentre,imcentre,
						bound,clust);
		}
	}
	free_real_vector(rc,0);
	free_real_vector(c,0);
	free_real_vector(rce,0);
	free_real_vector(clust,1);
}

void kcluster(int k, int m, int n, real_t rc[], real_t re[],
					real_t im[], real_t recentre[], real_t imcentre[],
					real_t bound[], real_t clust[])
{
	/* this function is used internally by BOUNDS */

	int i,stop,l,nonzero;
	real_t recent,imcent,d,prod,rad,gr,r,*dist,s,h1,h2,temp1,temp2;

	dist=allocate_real_vector(m,m+k-1);
	recent=re[m];
	imcent=im[m];
	stop=m+k-1;
	l = (imcent == 0.0) ? 0 : ((imcent > 0.0) ? 1 : -1);
	nonzero = (l != 0);
	for (i=m+1; i<=stop; i++) {
		recent += re[i];
		if (nonzero) {
			nonzero=(l == ((im[i] == 0.0) ? 0 :
						((im[i]>0.0) ? 1 : -1)));
			imcent += im[i];
		}
	}
	recent /= k;
	imcent = (nonzero ? imcent/k : 0.0);
	d=0.0;
	rad=0.0;
	for (i=m; i<=stop; i++) {
		recentre[i]=recent;
		imcentre[i]=imcent;
		temp1=re[i]-recent;
		temp2=im[i]-imcent;
		dist[i]=sqrt(temp1*temp1+temp2*temp2);
		if (d < dist[i]) d=dist[i];
	}
	s=sqrt(recent*recent+imcent*imcent);
	h1=rc[1];
	h2=rc[0];
	for (i=2; i<=n; i++) h1=h1*s+rc[i];
	for (i=1; i<=m-1; i++) {
		temp1=re[i]-recent;
		temp2=im[i]-imcent;
		h2 *= fabs(sqrt(temp1*temp1+temp2*temp2));
	}
	for (i=m+k; i<=n; i++) {
		temp1=re[i]-recent;
		temp2=im[i]-imcent;
		h2 *= fabs(sqrt(temp1*temp1+temp2*temp2));
	}
	gr=fabs((h1 == 0.0) ? 0.0 : ((h2 == 0.0) ? 10.0 : h1/h2));
	if (gr > 0.0)
		do {
			r=rad;
			rad=d+exp(log(1.1*gr)/k);
			if (rad == r) rad *= exp(log(1.1)/k);
			s=sqrt(recent*recent+imcent*imcent)+rad;
			h1=rc[1];
			h2=rc[0];
			for (i=2; i<=n; i++) h1=h1*s+rc[i];
			for (i=1; i<=m-1; i++) {
				temp1=re[i]-recent;
				temp2=im[i]-imcent;
				h2 *= fabs(sqrt(temp1*temp1+temp2*temp2)-rad);
			}
			for (i=m+k; i<=n; i++) {
				temp1=re[i]-recent;
				temp2=im[i]-imcent;
				h2 *= fabs(sqrt(temp1*temp1+temp2*temp2)-rad);
			}
			gr=(h1 == 0.0) ? 0.0 : ((h2 == 0.0) ? -10.0 : h1/h2);
			prod=1.0;
			for (i=m; i<=stop; i++) prod *= (rad-dist[i]);
		} while (prod <= gr);
	for (i=m; i<=stop; i++) {
		bound[i]=rad;
		clust[i]=k;
	}
	free_real_vector(dist,m);
}
