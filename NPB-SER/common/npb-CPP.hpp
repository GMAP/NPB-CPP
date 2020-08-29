/**
 * NASA Advanced Supercomputing Parallel Benchmarks C++
 *
 * based on NPB 3.3.1
 *
 * original version and technical report:
 * http://www.nas.nasa.gov/Software/NPB/
 *
 * C++ version:
 *      Dalvan Griebler <dalvangriebler@gmail.com>
 *      Gabriell Alves de Araujo <hexenoften@gmail.com>
 *      Júnior Löff <loffjh@gmail.com>
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>

typedef int boolean;
typedef struct { double real; double imag; } dcomplex;

#define TRUE	1
#define FALSE	0

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define	pow2(a) ((a)*(a))

/* old version of the complex number operations */
#define get_real(c) c.real
#define get_imag(c) c.imag
#define cadd(c,a,b) (c.real = a.real + b.real, c.imag = a.imag + b.imag)
#define csub(c,a,b) (c.real = a.real - b.real, c.imag = a.imag - b.imag)
#define cmul(c,a,b) (c.real = a.real * b.real - a.imag * b.imag, \
                     c.imag = a.real * b.imag + a.imag * b.real)
#define crmul(c,a,b) (c.real = a.real * b, c.imag = a.imag * b)

/* latest version of the complex number operations */
#define dcomplex_create(r,i) (dcomplex){r, i}
#define dcomplex_add(a,b) (dcomplex){(a).real+(b).real, (a).imag+(b).imag}
#define dcomplex_sub(a,b) (dcomplex){(a).real-(b).real, (a).imag-(b).imag}
#define dcomplex_mul(a,b) (dcomplex){((a).real*(b).real)-((a).imag*(b).imag),\
	((a).real*(b).imag)+((a).imag*(b).real)}
#define dcomplex_mul2(a,b) (dcomplex){(a).real*(b), (a).imag*(b)}
static inline dcomplex dcomplex_div(dcomplex z1, dcomplex z2){
	double a = z1.real;
	double b = z1.imag;
	double c = z2.real;
	double d = z2.imag;
	double divisor = c*c + d*d;
	double real = (a*c + b*d) / divisor;
	double imag = (b*c - a*d) / divisor;
	dcomplex result = (dcomplex){real, imag};
	return result;
}
#define dcomplex_div2(a,b) (dcomplex){(a).real/(b), (a).imag/(b)}
#define dcomplex_abs(x)    sqrt(((x).real*(x).real) + ((x).imag*(x).imag))
#define dconjg(x)          (dcomplex){(x).real, -1.0*(x).imag}

extern double randlc(double *, double);
extern void vranlc(int, double *, double, double *);
extern void timer_clear(int);
extern void timer_start(int);
extern void timer_stop(int);
extern double timer_read(int);

extern void c_print_results(char* name,
		char class_npb,
		int n1,
		int n2,
		int n3,
		int niter,
		double t,
		double mops,
		char* optype,
		int passed_verification,
		char* npbversion,
		char* compiletime,
		char* compilerversion,
		char* cc,
		char* clink,
		char* c_lib,
		char* c_inc,
		char* cflags,
		char* clinkflags,
		char* rand);
