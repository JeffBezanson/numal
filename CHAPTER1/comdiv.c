#include "../real.h"


void comdiv(real_t xr, real_t xi, real_t yr, real_t yi, real_t *zr, real_t *zi)
{
	real_t h,d;

	if (fabs(yi) < fabs(yr)) {
		if (yi == 0.0) {
			*zr=xr/yr;
			*zi=xi/yr;
		} else {
			h=yi/yr;
			d=h*yi+yr;
			*zr=(xr+h*xi)/d;
			*zi=(xi-h*xr)/d;
		}
	} else {
		h=yr/yi;
		d=h*yr+yi;
		*zr=(xr*h+xi)/d;
		*zi=(xi*h-xr)/d;
	}
}
