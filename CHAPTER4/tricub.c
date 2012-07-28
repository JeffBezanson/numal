#include "../real.h"


real_t tricub(real_t xi, real_t yi, real_t xj, real_t yj, real_t xk,
				real_t yk, real_t (*g)(real_t,real_t), real_t re, real_t ae)
{
	real_t tricubint(real_t, real_t, real_t, real_t,
			real_t, real_t, real_t, real_t, real_t,
			real_t, real_t, real_t, real_t, real_t,
			real_t, real_t, real_t, real_t, real_t,
			real_t, real_t, real_t (*)(real_t, real_t),
			real_t, real_t, real_t *, real_t);
	real_t surf,surfmin,xz,yz,gi,gj,gk;

	surf=0.5*fabs(xj*yk-xk*yj+xi*yj-xj*yi+xk*yi-xi*yk);
	surfmin=surf*re;
	re *= 30.0;
	ae=30.0*ae/surf;
	xz=(xi+xj+xk)/3.0;
	yz=(yi+yj+yk)/3.0;
	gi=(*g)(xi,yi);
	gj=(*g)(xj,yj);
	gk=(*g)(xk,yk);
	xi *= 0.5;
	yi *= 0.5;
	xj *= 0.5;
	yj *= 0.5;
	xk *= 0.5;
	yk *= 0.5;
	return tricubint(xi,yi,gi,xj,yj,gj,xk,yk,gk,
						xj+xk,yj+yk,(*g)(xj+xk,yj+yk),
						xk+xi,yk+yi,(*g)(xk+xi,yk+yi),
						xi+xj,yi+yj,(*g)(xi+xj,yi+yj),
						0.5*xz,0.5*yz,(*g)(xz,yz),
						g,re,ae,&surf,surfmin)/60.0;
}

real_t tricubint(real_t ax1, real_t ay1, real_t af1, real_t ax2,
			real_t ay2, real_t af2, real_t ax3, real_t ay3, real_t af3,
			real_t bx1, real_t by1, real_t bf1, real_t bx2, real_t by2,
			real_t bf2, real_t bx3, real_t by3, real_t bf3, real_t px,
			real_t py, real_t pf, real_t (*g)(real_t, real_t),
			real_t re, real_t ae, real_t *surf, real_t surfmin)
{
	/* this function is internally used by TRICUB */

	real_t e,i3,i4,i5,a,b,c,sx1,sy1,sx2,sy2,sx3,sy3,cx1,cy1,cf1,
			cx2,cy2,cf2,cx3,cy3,cf3,dx1,dy1,df1,dx2,dy2,df2,dx3,dy3,
			df3,result;

	a=af1+af2+af3;
	b=bf1+bf2+bf3;
	i3=3.0*a+27.0*pf+8.0*b;
	e=fabs(i3)*re+ae;
	if (*surf < surfmin || fabs(5.0*a+45.0*pf-i3) < e)
		return i3*(*surf);
	else {
		cx1=ax1+px;
		cy1=ay1+py;
		cf1=(*g)(cx1,cy1);
		cx2=ax2+px;
		cy2=ay2+py;
		cf2=(*g)(cx2,cy2);
		cx3=ax3+px;
		cy3=ay3+py;
		cf3=(*g)(cx3,cy3);
		c=cf1+cf2+cf3;
		i4=a+9.0*pf+4.0*b+12.0*c;
		if (fabs(i3-i4) < e)
			return i4*(*surf);
		else {
			sx1=0.5*bx1;
			sy1=0.5*by1;
			dx1=ax1+sx1;
			dy1=ay1+sy1;
			df1=(*g)(dx1,dy1);
			sx2=0.5*bx2;
			sy2=0.5*by2;
			dx2=ax2+sx2;
			dy2=ay2+sy2;
			df2=(*g)(dx2,dy2);
			sx3=0.5*bx3;
			sy3=0.5*by3;
			dx3=ax3+sx3;
			dy3=ay3+sy3;
			df3=(*g)(dx3,dy3);
			i5=(51.0*a+2187.0*pf+276.0*b+972.0*c-
					768.0*(df1+df2+df3))/63.0;
			if (fabs(i4-i5) < e)
				return i5*(*surf);
			else {
				(*surf) *= 0.25;
				result=tricubint(sx1,sy1,bf1,sx2,sy2,bf2,sx3,sy3,bf3,
							dx1,dy1,df1,dx2,dy2,df2,dx3,dy3,df3,px,py,pf,
							g,re,ae,surf,surfmin)+
						tricubint(ax1,ay1,af1,sx3,sy3,bf3,sx2,sy2,bf2,
							dx1,dy1,df1,ax1+sx2,ay1+sy2,
							(*g)(ax1+sx2,ay1+sy2),ax1+sx3,ay1+sy3,
							(*g)(ax1+sx3,ay1+sy3),0.5*cx1,0.5*cy1,cf1,
							g,re,ae,surf,surfmin)+
						tricubint(ax2,ay2,af2,sx3,sy3,bf3,sx1,sy1,bf1,
							dx2,dy2,df2,ax2+sx1,ay2+sy1,
							(*g)(ax2+sx1,ay2+sy1),ax2+sx3,ay2+sy3,
							(*g)(ax2+sx3,ay2+sy3),0.5*cx2,0.5*cy2,cf2,
							g,re,ae,surf,surfmin)+
						tricubint(ax3,ay3,af3,sx1,sy1,bf1,sx2,sy2,bf2,
							dx3,dy3,df3,ax3+sx2,ay3+sy2,
							(*g)(ax3+sx2,ay3+sy2),ax3+sx1,ay3+sy1,
							(*g)(ax3+sx1,ay3+sy1),0.5*cx3,0.5*cy3,cf3,
							g,re,ae,surf,surfmin);
				(*surf) *= 4.0;
				return result;
			}
		}
	}
}

