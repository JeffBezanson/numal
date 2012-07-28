#ifndef __REAL_H_
#define __REAL_H_

#ifdef __INTEL_COMPILER_
#include <mathimf.h>
#define M_PI   3.14159265358979323846
#define M_E    2.71828182845904523536
#define M_LN2  0.69314718055994530942
#define M_LN10 2.30258509299404568402
#else
#include <math.h>
#endif

#ifdef USE_DOUBLE

#ifndef HAVE_REAL_T
typedef double real_t;
#define HAVE_REAL_T
#endif

#define FLT_EPSILON 2.2204460492503131e-16
#define FLT_MAX 1.7976931348623157e+308
#define FLT_MIN 2.2250738585072014e-308

#else

#ifndef HAVE_REAL_T
typedef float real_t;
#define HAVE_REAL_T
#endif

#define FLT_EPSILON 1.1920928e-7
#define FLT_MAX 3.402823466e+38
#define FLT_MIN 1.175494351e-38

#endif

#endif
