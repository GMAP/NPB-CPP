/*
MIT License

Copyright (c) 2021 Parallel Applications Modelling Group - GMAP 
	GMAP website: https://gmap.pucrs.br
	
	Pontifical Catholic University of Rio Grande do Sul (PUCRS)
	Av. Ipiranga, 6681, Porto Alegre - Brazil, 90619-900

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

------------------------------------------------------------------------------

The original NPB 3.4.1 version was written in Fortran and belongs to: 
	http://www.nas.nasa.gov/Software/NPB/

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>
*/ 

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <tuple>
#include <complex>
#include <algorithm>
#include <execution>
#include <span>
#include <ranges>
#include <numeric>

typedef int boolean;
typedef struct { double real; double imag; } dcomplex;

class CountIterator {
  public:
    CountIterator(const int size){
      iter_vec.resize(size);
      std::iota(iter_vec.begin(), iter_vec.end(), 0);
    };

    std::vector<int>::iterator front(){
      return iter_vec.begin();
    }

    std::vector<int>::iterator tail(){
      return iter_vec.end();
    }

    std::vector<int>::reverse_iterator rfront(){
      return iter_vec.rbegin();
    }

    std::vector<int>::reverse_iterator rtail(){
      return iter_vec.rend();
    }

  private:
    std::vector<int> iter_vec;
};

#define TRUE	1
#define FALSE	0

#define max_npb(a,b) (((a) > (b)) ? (a) : (b))
#define min_npb(a,b) (((a) < (b)) ? (a) : (b))
#define	pow2(a) ((a)*(a))

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

extern double randlc(double &, double);
extern void vranlc(int, double &, double, double *);
extern void timer_clear(int);
extern void timer_start(int);
extern void timer_stop(int);
extern double timer_read(int);

extern void c_print_results(std::string name,
		char class_npb,
		int n1,
		int n2,
		int n3,
		int niter,
		double t,
		double mops,
		std::string optype,
		int passed_verification,
		std::string npbversion,
		std::string compiletime,
		std::string compilerversion,
		std::string cc,
		std::string clink,
		std::string c_lib,
		std::string c_inc,
		std::string cflags,
		std::string clinkflags,
		std::string rand);
