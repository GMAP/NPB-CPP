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

Authors of the Fortran code:
	E. Barszcz
	P. Frederickson
	A. Woo
	M. Yarrow

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>
	Arthur S. Bianchessi <arthur.bianchessi@edu.pucrs.br>
	Leonardo Mallmann <leonardo.mallmann@edu.pucrs.br>
*/

#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"

#define NM (2+(1<<LM)) /* actual dimension including ghost cells for communications */
#define NV (ONE*(2+(1<<NDIM1))*(2+(1<<NDIM2))*(2+(1<<NDIM3))) /* size of rhs array */
#define NR (((NV+NM*NM+5*NM+7*LM+6)/7)*8) /* size of residual array */
#define MAXLEVEL (LT_DEFAULT+1) /* maximum number of levels */
#define M (NM+1) /* set at m=1024, can handle cases up to 1024^3 case */
#define MM (10)
#define	A (pow(5.0,13.0))
#define	X (314159265.0)
#define T_INIT 0
#define T_BENCH 1
#define T_MG3P 2
#define T_PSINV 3
#define T_RESID 4
#define T_RESID2 5
#define T_RPRJ3 6
#define T_INTERP 7
#define T_NORM2 8
#define T_COMM3 9
#define T_LAST 10

/* global variables */
std::vector<int> nx(MAXLEVEL+1, 0);
std::vector<int> ny(MAXLEVEL+1, 0);
std::vector<int> nz(MAXLEVEL+1, 0);
std::vector<int> m1(MAXLEVEL+1, 0);
std::vector<int> m2(MAXLEVEL+1, 0);
std::vector<int> m3(MAXLEVEL+1, 0);
std::vector<int> ir(MAXLEVEL+1, 0);
std::vector<int> debug_vec(8, 0);
std::vector<double> u(NR);
std::vector<double> v(NV);
std::vector<double> r(NR);
static int is1, is2, is3, ie1, ie2, ie3, lt, lb;
static boolean timeron;

CountIterator iter(X);
auto policy = std::execution::par;

struct Returnable {
	double s, rnmu;
	Returnable(double s, double rnmu) : s(s), rnmu(rnmu) {}

	Returnable operator+(const Returnable& other) const {
		return Returnable(s + other.s, std::max(rnmu, other.rnmu));
	}
};

/* function prototypes */
static void bubble(std::vector<double> &ten, std::vector<int> &j1, std::vector<int> &j2, std::vector<int> &j3, int m, int ind);
static void comm3(std::span<double> u, int n1, int n2, int n3);
static void interp(std::span<double> z, int mm1, int mm2, int mm3, std::span<double> u, int n1, int n2, int n3, int k);
static void mg3P(std::vector<double> &u, std::vector<double> &v, std::vector<double> &r, std::vector<double> &a, std::vector<double> &c, int n1, int n2, int n3, int k);
static void norm2u3(std::span<double> r, int n1, int n2, int n3, double* rnm2, double* rnmu, int nx, int ny, int nz);
static double power(double a, int n);
static void psinv(std::span<double> r, std::span<double> u, int n1, int n2, int n3, std::vector<double> &c, int k);
static void rep_nrm(std::span<double> u, int n1, int n2, int n3, std::string title, int kk);
static void resid(std::span<double> u, std::span<double> v, std::span<double> r, int n1, int n2, int n3, std::vector<double> &a, int k);
static void rprj3(std::span<double> r, int m1k, int m2k, int m3k, std::span<double> s, int m1j, int m2j, int m3j, int k);
static void setup(int* n1, int* n2, int* n3, int k);
static void showall(std::span<double> z, int n1, int n2, int n3);
static void zero3(std::span<double> z, int n1, int n2, int n3);
static void zran3(std::vector<double> &z, int n1, int n2, int n3, int nx, int ny, int k);

/* mg */
int main(int argc, char *argv[]){
	/*
	 * -------------------------------------------------------------------------
	 * k is the current level. it is passed down through subroutine args
	 * and is not global. it is the current iteration
	 * -------------------------------------------------------------------------
	 */
	int k;
	double t, tinit, mflops;
	std::vector<double> a(4), c(4);

	double rnm2, rnmu, epsilon;
	int n1, n2, n3, nit;
	double nn, verify_value, err;
	boolean verified;
	char class_npb;

	std::vector<std::string> t_names(10);
	double tmax;

	for(int i=T_INIT; i<T_LAST; i++){
		timer_clear(i);
	}

	timer_start(T_INIT);	

	/*
	 * ----------------------------------------------------------------------
	 * read in and broadcast input data
	 * ----------------------------------------------------------------------
	 */
	std::fstream fp;
	fp.open("timer.flag", std::ios::in);
	if(fp.is_open()){
		timeron = TRUE;
		t_names[T_INIT] = "init";
		t_names[T_BENCH] = "benchmk";
		t_names[T_MG3P] = "mg3P";
		t_names[T_PSINV] = "psinv";
		t_names[T_RESID] = "resid";
		t_names[T_RPRJ3] = "rprj3";
		t_names[T_INTERP] = "interp";
		t_names[T_NORM2] = "norm2";
		t_names[T_COMM3] = "comm3";
		fp.close();
	}else{
		timeron = FALSE;
	}
	fp.open("mg.input", std::ios::in);
	if(fp.is_open()){
		std::cout
      		<< " Reading from input file mg.input"
      		<< std::endl;
		if(!(fp >> lt)){
			std::cout
				<< " Error in reading elements"
				<< std::endl;
			exit(1);
		}
		while(fp.get() != '\n');
		if(!(fp >> nx[lt] >> ny[lt] >> nz[lt])){
			std::cout
				<< " Error in reading elements"
				<< std::endl;
			exit(1);
		}
		while(fp.get() != '\n');
		if(!(fp >> nit)){
			std::cout
				<< " Error in reading elements"
				<< std::endl;
			exit(1);
		}
		while(fp.get() != '\n');
		std::for_each(debug_vec.begin(), debug_vec.end(), [&fp](int &debug){
			if(!(fp >> debug)){
				std::cout
					<< " Error in reading elements"
					<< std::endl;
				exit(1);
			}
		});
		fp.close();
	}else{
		std::cout
      << " No input file. Using compiled defaults"
      << std::endl;
		lt = LT_DEFAULT;
		nit = NIT_DEFAULT;
		nx[lt] = NX_DEFAULT;
		ny[lt] = NY_DEFAULT;
		nz[lt] = NZ_DEFAULT;
		std::for_each(debug_vec.begin(), debug_vec.end(), [&fp](int &debug){
			debug = DEBUG_DEFAULT;
		});
	}

	if((nx[lt] != ny[lt]) || (nx[lt] != nz[lt])){
		class_npb = 'U';
	}else if(nx[lt] == 32 && nit == 4){
		class_npb = 'S';
	}else if(nx[lt] == 128 && nit == 4){
		class_npb = 'W';
	}else if(nx[lt] == 256 && nit == 4){
		class_npb = 'A';
	}else if(nx[lt] == 256 && nit == 20){
		class_npb = 'B';
	}else if(nx[lt] == 512 && nit == 20){
		class_npb = 'C';
	}else if(nx[lt] == 1024 && nit == 50){  
		class_npb = 'D';
	}else if(nx[lt] == 2048 && nit == 50){  
		class_npb = 'E';	
	}else{
		class_npb = 'U';
	}

	/*
	 * ---------------------------------------------------------------------
	 * use these for debug info:
	 * ---------------------------------------------------------------------
	 * debug_vec(0) = 1 !=> report all norms
	 * debug_vec(1) = 1 !=> some setup information
	 * debug_vec(1) = 2 !=> more setup information
	 * debug_vec(2) = k => at level k or below, show result of resid
	 * debug_vec(3) = k => at level k or below, show result of psinv
	 * debug_vec(4) = k => at level k or below, show result of rprj
	 * debug_vec(5) = k => at level k or below, show result of interp
	 * debug_vec(6) = 1 => (unused)
	 * debug_vec(7) = 1 => (unused)
	 * ---------------------------------------------------------------------
	 */
	a[0] = -8.0/3.0;
	a[1] =  0.0;
	a[2] =  1.0/6.0;
	a[3] =  1.0/12.0;

	if(class_npb == 'A' || class_npb == 'S' || class_npb =='W'){
		/* coefficients for the s(a) smoother */
		c[0] =  -3.0/8.0;
		c[1] =  +1.0/32.0;
		c[2] =  -1.0/64.0;
		c[3] =   0.0;
	}else{
		/* coefficients for the s(b) smoother */
		c[0] =  -3.0/17.0;
		c[1] =  +1.0/33.0;
		c[2] =  -1.0/61.0;
		c[3] =   0.0;
	}

	lb = 1;
	k = lt;

	setup(&n1,&n2,&n3,k);

	// zero3(u,n1,n2,n3);
	std::fill(u.begin(), u.end(), 0.0);
	zran3(v,n1,n2,n3,nx[lt],ny[lt],k);

	norm2u3(v,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	std::cout
    << std::endl
    << std::endl
    << "NAS Parallel Benchmarks 4.1 Serial C++ version - MG Benchmark"
    << std::endl
    << std::endl;
	std::cout
    << " Size: "
    << std::setw(3) << nx[lt] << "x"
    << std::setw(3) << ny[lt] << "x"
    << std::setw(3) << nz[lt]
    << " (class_npb " << std::setw(1) << class_npb << ")"
    << std::endl;
	std::cout
    << " Iterations: "
    << std::setw(3) << nit
    << std::endl;

	resid(u,v,r,n1,n2,n3,a,k);
	norm2u3(r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	/*
	 * ---------------------------------------------------------------------
	 * one iteration for startup
	 * ---------------------------------------------------------------------
	 */
	mg3P(u,v,r,a,c,n1,n2,n3,k);
	resid(u,v,r,n1,n2,n3,a,k);

	setup(&n1,&n2,&n3,k);

	std::fill(u.begin(), u.end(), 0.0);
	zran3(v,n1,n2,n3,nx[lt],ny[lt],k);

	timer_stop(T_INIT);
	tinit = timer_read(T_INIT);
	std::cout
		<< " Initialization time: "
		<< std::setw(15) << std::setprecision(3) << std::fixed << tinit
		<< " seconds"
		<< std::endl;

	for(int i=T_BENCH; i<T_LAST; i++){
		timer_clear(i);
	}
	timer_start(T_BENCH);

	if(timeron){timer_start(T_RESID2);}
	resid(u,v,r,n1,n2,n3,a,k);
	if(timeron){timer_stop(T_RESID2);}
	norm2u3(r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	for(int it = 1; it <= nit; it++){
		if((it==1)||(it==nit)||((it%5)==0)){
	  		std::cout
				<< "  iter "
				<< std::setw(3) << it
				<< std::endl;
		}
		if(timeron){timer_start(T_MG3P);}
		mg3P(u,v,r,a,c,n1,n2,n3,k);
		if(timeron){timer_stop(T_MG3P);}
		if(timeron){timer_start(T_RESID2);}
		resid(u,v,r,n1,n2,n3,a,k);
		if(timeron){timer_stop(T_RESID2);}
	}
	norm2u3(r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	timer_stop(T_BENCH);
	t = timer_read(T_BENCH);    	

	verified = FALSE;
	verify_value = 0.0;	

	std::cout
		<< " Benchmark completed"
		<< std::endl;

	epsilon = 1.0e-8;
	if(class_npb != 'U'){
		if(class_npb == 'S'){
			verify_value = 0.5307707005734e-04;
		}else if(class_npb == 'W'){
			verify_value = 0.6467329375339e-05;
		}else if(class_npb == 'A'){
			verify_value = 0.2433365309069e-05;
		}else if(class_npb == 'B'){
			verify_value = 0.1800564401355e-05;
		}else if(class_npb == 'C'){
			verify_value = 0.5706732285740e-06;
		}else if(class_npb == 'D'){
			verify_value = 0.1583275060440e-09;
		}else if(class_npb == 'E'){
			verify_value = 0.8157592357404e-10; 
		}

		err = fabs(rnm2-verify_value) / verify_value;
		if(err <= epsilon){
			verified = TRUE;
			std::cout
				<< " VERIFICATION SUCCESSFUL" 
				<< std::endl;
			std::cout
				<< " L2 Norm is "
				<< std::setw(20) << std::setprecision(13) << std::scientific << rnm2
				<< std::endl;
			std::cout
				<< " Error is   "
				<< std::setw(20) << std::setprecision(13) << std::scientific << err
				<< std::endl;
		}else{
			verified = FALSE;
			std::cout
        << " VERIFICATION FAILED"
        << std::endl;
			std::cout
        << " L2 Norm is             "
        << std::setw(20) << std::setprecision(13) << std::scientific << rnm2
        << std::endl;
			std::cout
        << " The correct L2 Norm is "
        << std::setw(20) << std::setprecision(13) << std::scientific << verify_value
        << std::endl;
		}
	}else{
		verified = FALSE;
		std::cout
      << " Problem size unknown"
      << std::endl;
		std::cout
      << " NO VERIFICATION PERFORMED"
      << std::endl;
	}

	nn = 1.0*nx[lt]*ny[lt]*nz[lt];

	if(t!=0.0){
		mflops = 58.0*nit*nn*1.0e-6/t;
	}else{
		mflops = 0.0;
	}

	c_print_results("MG",
			class_npb,
			nx[lt],
			ny[lt],
			nz[lt],
			nit,
			t,
			mflops,
			"          floating point",
			verified,
			NPBVERSION,
			COMPILETIME,
			COMPILERVERSION,
			CS1,
			CS2,
			CS3,
			CS4,
			CS5,
			CS6,
			CS7);

	/*
	 * ---------------------------------------------------------------------
	 * more timers
	 * ---------------------------------------------------------------------
	 */
	if(timeron){
		tmax = timer_read(T_BENCH);
		if(tmax==0.0){tmax=1.0;}
		std::cout
			<< "  SECTION   Time (secs)"
			<< std::endl;
		for(int i=T_BENCH; i<T_LAST; i++){
			t = timer_read(i);
			if(i==T_RESID2){
				t = timer_read(T_RESID) - t;
				std::cout
					<< "    --> "
					<< std::setw(8) << "mg-resid" << ":"
					<< std::setw(9) << std::setprecision(3) << t
					<< "  (" << std::setw(6) << std::setprecision(2) << t*100.0/tmax << "%)"
					<< std::endl;
			}else{
				std::cout
					<< "  " << std::left << std::setw(8) << t_names[i] << ":"
					<< std::right << std::setw(9) << std::setprecision(3) << t
					<< "  (" << std::setw(6) << std::setprecision(2) << t*100.0/tmax << "%)"
					<< std::endl;
			}
		}
	}

	return 0;
}

/*
 * ---------------------------------------------------------------------
 * bubble does a bubble sort in direction dir
 * ---------------------------------------------------------------------
 */
static void bubble(std::vector<double> &ten, std::vector<int> &j1, std::vector<int> &j2, std::vector<int> &j3, int m, int ind){
	double temp, MMI, MMI1;
	int j_temp;

	if (ind == 1) {
		for(int i = 0; i < m-1; i++){
			MMI = MM + i;
			MMI1 = MM + (i + 1);
			if (ten[MMI] > ten[MMI1]) {
				std::swap(ten[MMI], ten[MMI1]);
				std::swap(j1[MMI], j1[MMI1]);
				std::swap(j2[MMI], j2[MMI1]);
				std::swap(j3[MMI], j3[MMI1]);
			} else { return; }
		}
	} else {
		for(int i = 0; i < m-1; i++){
			if (ten[i] < ten[i+1]) {
				std::swap(ten[i], ten[i+1]);
				std::swap(j1[i], j1[i+1]);
				std::swap(j2[i], j2[i+1]);
				std::swap(j3[i], j3[i+1]);
			}
			else { return; }
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * comm3 organizes the communication on all borders 
 * ---------------------------------------------------------------------
 */
static void comm3(std::span<double> u, int n1, int n2, int n3){
	if (timeron) {timer_start(T_COMM3);}
	std::for_each_n(policy, iter.front() + 1, n3 - 2, [&] (int i3){
		int k3 = i3 * n1 * n2;

		/* axis = 1 */
		for(int i2 = 1; i2 < n2-1; i2++){
			int k2 = i2 * n1;
			u[k3 + k2] = u[k3 + k2 + (n1 - 2)];
			u[k3 + k2 + (n1 - 1)] = u[k3 + k2 + 1];
		}

		int k2 = n1 * (n2-1);
		/* axis = 2 */
		for(int i1 = 0; i1 < n1; i1++){
			u[k3 + i1] = u[k3 + (k2 - n1) + i1];
			u[k3 + k2 + i1] = u[k3 + n1 + i1];
		}
	});

	/* axis = 3 */
	std::for_each_n(policy, iter.front(), n2, [&] (int i2){
		int k3 = n1 * n2;
		int k3m = k3 * (n3 - 1);
		int k3mm = k3 * (n3 - 2);
		int k2 = i2 * n1;
		for(int i1 = 0; i1 < n1; i1++){
			u[k2 + i1] = u[k3mm + k2 + i1];
			u[k3m + k2 + i1] = u[k3 + k2 + i1];
		}
	});
	if (timeron) {timer_stop(T_COMM3);}
}

/*
 * --------------------------------------------------------------------
 * interp adds the trilinear interpolation of the correction
 * from the coarser grid to the current approximation: u = u + Qu'
 *     
 * observe that this  implementation costs  16A + 4M, where
 * A and M denote the costs of addition and multiplication.  
 * note that this vectorizes, and is also fine for cache 
 * based machines. vector machines may get slightly better 
 * performance however, with 8 separate "do i1" loops, rather than 4.
 * --------------------------------------------------------------------
 */
static void interp(std::span<double> z, int mm1, int mm2, int mm3, std::span<double> u, int n1, int n2, int n3, int k){
	int d1, d2, d3, t1, t2, t3;

	/* 
	 * --------------------------------------------------------------------
	 * note that m = 1037 in globals.h but for this only need to be
	 * 535 to handle up to 1024^3
	 * integer m
	 * parameter( m=535 )
	 * --------------------------------------------------------------------
	 */

	if(timeron){timer_start(T_INTERP);}
	if(n1 != 3 && n2 != 3 && n3 != 3){
		std::for_each_n(policy, iter.front(), mm3 - 1, [&] (int i3){
			std::vector<double> z1(M), z2(M), z3(M);
  			int kz2, kz3, ku2, ku3;

			kz3 = i3*mm2*mm1;
			ku3 = i3*n2*n1*2;
			for(int i2 = 0; i2 < mm2-1; i2++){
				kz2 = i2*mm1;
				ku2 = i2*n1*2;
				for(int i1 = 0; i1 < mm1; i1++){
					z1[i1] = z[kz3+(kz2+mm1)+i1] + z[kz3+kz2+i1];
					z2[i1] = z[(kz3+mm1*mm2)+kz2+i1] + z[kz3+kz2+i1];
					z3[i1] = z[(kz3+mm1*mm2)+(kz2+mm1)+i1] + z[(kz3+mm1*mm2)+kz2+i1] + z1[i1];
				}
				for(int i1 = 0; i1 < mm1-1; i1++){
					u[(ku3+ku2+(i1*2))] += z[kz3+kz2+i1];
					u[(ku3+ku2+(i1*2+1))] += 0.5*(z[kz3+kz2+i1+1]+z[kz3+kz2+i1]);
				}
				for(int i1 = 0; i1 < mm1-1; i1++){
					u[(ku3)+(ku2+n1)+(i1*2)] += 0.5 * z1[i1];
					u[(ku3)+(ku2+n1)+(i1*2+1)] += 0.25*( z1[i1] + z1[i1+1] );
				}
				for(int i1 = 0; i1 < mm1-1; i1++){
					u[(ku3+n1*n2)+(ku2)+(i1*2)] += 0.5 * z2[i1];
					u[(ku3+n1*n2)+(ku2)+(i1*2+1)] += 0.25*( z2[i1] + z2[i1+1] );
				}
				for(int i1 = 0; i1 < mm1-1; i1++){
					u[(ku3+n1*n2)+(ku2+n1)+(i1*2)] += 0.25 * z3[i1];
					u[(ku3+n1*n2)+(ku2+n1)+(i1*2+1)] += 0.125*( z3[i1] + z3[i1+1] );
				}
			}
		});
	}else{
		if(n1 == 3){
			d1 = 2;
			t1 = 1;
		}else{
			d1 = 1;
			t1 = 0;
		}      
		if(n2 == 3){
			d2 = 2;
			t2 = 1;
		}else{
			d2 = 1;
			t2 = 0;
		}          
		if(n3 == 3){
			d3 = 2;
			t3 = 1;
		}else{
			d3 = 1;
			t3 = 0;
		}
		std::for_each_n(policy, iter.front()+d3, mm3-d3, [&] (int i3) {
  			int kz2, kz3, ku2, ku3;
			kz3 = (i3-1)*mm1*mm2;
			ku3 = (i3-d3-1)*n1*n2*2;
			for(int i2 = d2; i2 <= mm2-1; i2++){
				kz2 = (i2-1)*n1;
				ku2 = (i2-d2-1)*n1;
				for(int i1 = d1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-d1-1)] += z[kz3+kz2+i1-1];
				}
				for(int i1 = 1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-t1-1)] += 0.5*(z[kz3+kz2+i1]+z[kz3+kz2+i1-1]);
				}
			}
			for(int i2 = 1; i2 <= mm2-1; i2++){
				kz2 = i2*n1;
				ku2 = (i2-t2-1)*n1;
				for (int i1 = d1; i1 <= mm1-1; i1++) {
					u[ku3+ku2+(2*i1-d1-1)] +=
						0.5*(z[kz3+kz2+i1-1]+z[kz3+kz2-n1+i1-1]);
				}
				for(int i1 = 1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-t1-1)] +=
						0.25*(z[kz3+kz2+i1]+z[kz3+kz2-n1+i1]
								+z[kz3+kz2+i1-1]+z[kz3+kz2-n1+i1-1]);
				}
			}
		});
		std::for_each_n(policy, iter.front()+1, mm3-1, [&] (int i3) {
  			int kz2, kz3, ku2, ku3;
			kz3 = (i3-1)*mm1*mm2;
			ku3 = (i3-t3-1)*n1*n2*2;
			for(int i2 = d2; i2 <= mm2-1; i2++){
				kz2 = (i2-1)*n1;
				ku2 = (i2-d2-1)*n1;
				for(int i1 = d1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-d1-1)] +=
						+0.5*(z[(kz3+mm1*mm2)+kz2+i1-1]+z[kz3+kz2+i1-1]);
				}
				for(int i1 = 1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-t1-1)] +=
						+0.25*(z[(kz3+mm1*mm2)+kz2+i1]+z[(kz3+mm1*mm2)+kz2+i1-1]
								+z[(kz3+mm1*mm2)+kz2+i1]+z[kz3+kz2+i1-1]);
				}
			}
			for(int i2 = 1; i2 <= mm2-1; i2++){
				kz2 = (i2-1)*n1;
				ku2 = (i2-t2-1)*n1;
				for (int i1 = d1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-d1-1)] +=
						0.25*(z[kz3+mm1*mm2+kz2+n1+i1-1]+z[kz3+mm1*mm2+kz2+i1-1]
								+z[kz3+kz2+n1+i1-1]+z[kz3+kz2+i1-1]);
				}
				for(int i1 = 1; i1 <= mm1-1; i1++){
					u[ku3+ku2+(2*i1-t1-1)] +=
						0.125*(z[kz3+mm1*mm2+kz2+n1+i1]+z[kz3+mm1*mm2+kz2+i1]
								+z[kz3+mm1*mm2+kz2+n1+i1-1]+z[kz3+mm1*mm2+kz2+i1-1]
								+z[kz3+kz2+n1+i1]+z[kz3+kz2+i1]
								+z[kz3+kz2+n1+i1-1]+z[kz3+kz2+i1-1]);
				}
			}
		});
	}
	if(timeron){timer_stop(T_INTERP);}

	if(debug_vec[0] >= 1){
		rep_nrm(z,mm1,mm2,mm3,"z: inter",k-1);
		rep_nrm(u,n1,n2,n3,"u: inter",k);
	}
	if(debug_vec[5] >= k){
		showall(z,mm1,mm2,mm3);
		showall(u,n1,n2,n3);
	}
}

/* 
 * --------------------------------------------------------------------
 * multigrid v-cycle routine
 * --------------------------------------------------------------------
 */
static void mg3P(std::vector<double> &u, std::vector<double> &v, std::vector<double> &r, std::vector<double> &a, std::vector<double> &c, int n1, int n2, int n3, int k){
	int j;

	/*
	 * --------------------------------------------------------------------
	 * down cycle.
	 * restrict the residual from the find grid to the coarse
	 * -------------------------------------------------------------------
	 */
	for(int k = lt; k >= lb+1; k--){
		j = k-1;
		rprj3({r.begin()+ir[k], r.end()}, m1[k], m2[k], m3[k],{r.begin()+ir[j], r.end()}, m1[j], m2[j], m3[j], k);
	}

	k = lb;
	/*
	 * --------------------------------------------------------------------
	 * compute an approximate solution on the coarsest grid
	 * --------------------------------------------------------------------
	 */
	std::fill_n(policy, u.begin()+ir[k], m1[k]*m2[k]*m3[k], 0.0);
	psinv({r.begin()+ir[k], r.end()}, {u.begin()+ir[k], u.end()}, m1[k], m2[k], m3[k], c, k);

	for(int k = lb+1; k <= lt-1; k++){
		j = k-1;
		/*
		 * --------------------------------------------------------------------
		 * prolongate from level k-1  to k
		 * -------------------------------------------------------------------
		 */
		std::fill_n(policy, u.begin()+ir[k], m1[k]*m2[k]*m3[k], 0.0);
		interp({u.begin()+ir[j], u.end()}, m1[j], m2[j], m3[j], {u.begin()+ir[k], u.end()}, m1[k], m2[k], m3[k], k);
		/*
		 * --------------------------------------------------------------------
		 * compute residual for level k
		 * --------------------------------------------------------------------
		 */
		resid({u.begin()+ir[k], u.end()}, {r.begin()+ir[k], r.end()}, {r.begin()+ir[k], r.end()}, m1[k], m2[k], m3[k], a, k);	
		/*
		 * --------------------------------------------------------------------
		 * apply smoother
		 * --------------------------------------------------------------------
		 */
		psinv({r.begin()+ir[k], r.end()}, {u.begin()+ir[k], u.end()}, m1[k], m2[k], m3[k], c, k);
	}

	j = lt - 1;
	k = lt;
	interp({u.begin()+ir[j], u.end()}, m1[j], m2[j], m3[j], u, n1, n2, n3, k);	
	resid(u, v, r, n1, n2, n3, a, k);	
	psinv(r, u, n1, n2, n3, c, k);
}

/*
 * ---------------------------------------------------------------------
 * norm2u3 evaluates approximations to the l2 norm and the
 * uniform (or l-infinity or chebyshev) norm, under the
 * assumption that the boundaries are periodic or zero. add the
 * boundaries in with half weight (quarter weight on the edges
 * and eighth weight at the corners) for inhomogeneous boundaries.
 * ---------------------------------------------------------------------
 */
static void norm2u3(std::span<double> r, int n1, int n2, int n3, double* rnm2, double* rnmu, int nx, int ny, int nz){
	double s;
	double dn;

	if(timeron){timer_start(T_NORM2);}
	dn = 1.0*nx*ny*nz;

	Returnable result = std::transform_reduce(policy, iter.front()+1, iter.front()+n3-1, Returnable(0.0, 0.0), 
	std::plus<Returnable>(),
	[&] (int i3) {
		int ku3, ku2;
		double s = 0, rnmu = 0;
		ku3 = i3*n1*n2;
		for(int i2 = 1; i2 < n2-1; i2++){
			ku2 = i2*n1;
			for(int i1 = 1; i1 < n1-1; i1++){
				s += r[ku3+ku2+i1] * r[ku3+ku2+i1];
				rnmu = std::max(fabs(r[ku3+ku2+i1]), rnmu);
			}
		}
		return Returnable(s, rnmu);
	});
	s 	  = result.s;
	*rnmu = result.rnmu;



	*rnm2 = sqrt(s/dn);
	if(timeron){timer_stop(T_NORM2);}
}

/*
 * ---------------------------------------------------------------------
 * power raises an integer, disguised as a double
 * precision real, to an integer power
 * ---------------------------------------------------------------------
 */
static double power(double a, int n){
	double aj;
	int nj;
	double rdummy;
	double power;

	power = 1.0;
	nj = n;
	aj = a;

	while(nj != 0){
		if((nj%2)==1){rdummy = randlc(power, aj);}
		rdummy = randlc(aj, aj);
		nj = nj/2;
	}

	return power;
}

/*
 * --------------------------------------------------------------------
 * psinv applies an approximate inverse as smoother: u = u + Cr
 * 
 * this  implementation costs  15A + 4M per result, where
 * A and M denote the costs of Addition and Multiplication.  
 * presuming coefficient c(3) is zero (the NPB assumes this,
 * but it is thus not a general case), 2A + 1M may be eliminated,
 * resulting in 13A + 3M.
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void psinv(std::span<double> r, std::span<double> u, int n1, int n2, int n3, std::vector<double> &c, int k){
	if(timeron){timer_start(T_PSINV);}
	std::for_each_n(policy, iter.front()+1, n3-2, [&] (int i3){
  		int kr2, kr3, kr3P, kr3M;
		std::vector<double> r1(M), r2(M);
		kr3 = i3*n2*n1;
		kr3P = (i3+1)*n2*n1;
		kr3M = (i3-1)*n2*n1;
		for(int i2 = 1; i2 < n2-1; i2++) {
      		kr2 = i2*n1;
			for(int i1 = 0; i1 < n1; i1++){
				r1[i1] = r[kr3+kr2-n1+i1] + r[kr3+kr2+n1+i1]
					+ r[kr3M+kr2+i1] + r[kr3P+kr2+i1];
				r2[i1] = r[kr3M+kr2-n1+i1] + r[kr3M+kr2+n1+i1]
					+ r[kr3P+kr2-n1+i1] + r[kr3P+kr2+n1+i1];
			}
			for(int i1 = 1; i1 < n1-1; i1++){
				u[kr3+kr2+i1] +=
					+ c[0] * r[kr3+kr2+i1]
					+ c[1] * ( r[kr3+kr2+i1-1] + r[kr3+kr2+i1+1]
							+ r1[i1] )
					+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
				/*
				 * --------------------------------------------------------------------
				 * assume c(3) = 0    (enable line below if c(3) not= 0)
				 * --------------------------------------------------------------------
				 * > + c(3) * ( r2(i1-1) + r2(i1+1) )
				 * --------------------------------------------------------------------
				 */
			}
		}
	});
	if(timeron){timer_stop(T_PSINV);}

	/*
	 * --------------------------------------------------------------------
	 * exchange boundary points
	 * --------------------------------------------------------------------
	 */
	comm3(u,n1,n2,n3);

	if(debug_vec[0] >= 1){
		rep_nrm(u,n1,n2,n3,"   psinv",k);
	}

	if(debug_vec[3] >= k){
		showall(u,n1,n2,n3);
	}
}

/*
 * ---------------------------------------------------------------------
 * report on norm
 * ---------------------------------------------------------------------
 */
static void rep_nrm(std::span<double> u, int n1, int n2, int n3, std::string title, int kk){
	double rnm2, rnmu;

	norm2u3(u,n1,n2,n3,&rnm2,&rnmu,nx[kk],ny[kk],nz[kk]);
	std::cout
    << (" Level%2d in %8s: norms =%21.14e%21.14e\n", kk, title, rnm2, rnmu)
    << std::setw(2) << kk
    << " in "
    << std::setw(8) << title << ":"
    << " norms ="
    << std::setw(21) << std::setprecision(14) << std::scientific << rnm2
    << std::setw(21) << std::setprecision(14) << std::scientific << rnmu
    << std::endl;
}

/*
 * --------------------------------------------------------------------
 * resid computes the residual: r = v - Au
 *
 * this  implementation costs  15A + 4M per result, where
 * A and M denote the costs of addition (or subtraction) and 
 * multiplication, respectively. 
 * presuming coefficient a(1) is zero (the NPB assumes this,
 * but it is thus not a general case), 3A + 1M may be eliminated,
 * resulting in 12A + 3M.
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void resid(std::span<double> u, std::span<double> v, std::span<double> r, int n1, int n2, int n3, std::vector<double> &a, int k){
	if(timeron){timer_start(T_RESID);}
	std::for_each_n(policy, iter.front()+1, n3-2, [&] (int i3) {
		int ku2, ku3, ku3P, ku3M;
		std::vector<double> u1(M), u2(M);

		ku3 = i3*n2*n1;
		ku3P = (i3+1)*n2*n1;
		ku3M = (i3-1)*n2*n1;
		for(int i2 = 1; i2 < n2-1; i2++){
      		ku2 = i2*n1;
			for(int i1 = 0; i1 < n1; i1++){
				u1[i1] = u[ku3+ku2-n1+i1] + u[ku3+ku2+n1+i1]
					+ u[ku3M+ku2+i1] + u[ku3P+ku2+i1];
				u2[i1] = u[ku3M+ku2-n1+i1] + u[ku3M+ku2+n1+i1]
					+ u[ku3P+ku2-n1+i1] + u[ku3P+ku2+n1+i1];
			}
			for(int i1 = 1; i1 < n1-1; i1++){
				r[ku3+ku2+i1] = v[ku3+ku2+i1]
					- a[0] * u[ku3+ku2+i1]
					/*
					 * ---------------------------------------------------------------------
					 * assume a(1) = 0 (enable 2 lines below if a(1) not= 0)
					 * ---------------------------------------------------------------------
					 * > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
					 * > + u1(i1) )
					 * ---------------------------------------------------------------------
					 */
					- a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
					- a[3] * ( u2[i1-1] + u2[i1+1] );
			}
		}
	});
	if(timeron){timer_stop(T_RESID);}

	/*
	 * --------------------------------------------------------------------
	 * exchange boundary data
	 * --------------------------------------------------------------------
	 */
	comm3(r,n1,n2,n3);

	if(debug_vec[0] >= 1){
		rep_nrm(r,n1,n2,n3,"   resid",k);
	}

	if(debug_vec[2] >= k){
		showall(r,n1,n2,n3);
	}
}

/*
 * --------------------------------------------------------------------
 * rprj3 projects onto the next coarser grid, 
 * using a trilinear finite element projection: s = r' = P r
 *     
 * this  implementation costs 20A + 4M per result, where
 * A and M denote the costs of addition and multiplication.  
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void rprj3(std::span<double> r, int m1k, int m2k, int m3k, std::span<double> s, int m1j, int m2j, int m3j, int k){
	int i3, i2, i1, d1, d2, d3, j;
  int ku2, ku3, ku3P, ku3P2;
  int ks2, ks3;

  double x2, y2;

	if(timeron){timer_start(T_RPRJ3);}
	if(m1k == 3){
		d1 = 2;
	}else{
		d1 = 1;
	}
	if(m2k == 3){
		d2 = 2;
	}else{
		d2 = 1;
	}
	if(m3k == 3){
		d3 = 2;
	}else{
		d3 = 1;
	}
	std::for_each_n(policy, iter.front()+1, m3j-2, [&] (int j3){
		int i3, i2, i1;
		int ku2, ku3, ku3P, ku3P2;
  		int ks2, ks3;
		std::vector<double> x1(M), y1(M);
		double x2, y2;

		i3 = 2*j3-d3;
		ku3 = i3*m1k*m2k;
		ku3P = (i3+1)*m1k*m2k;
		ku3P2 = (i3+2)*m1k*m2k;
		ks3 = j3*m1j*m2j;
		for(int j2 = 1; j2 < m2j-1; j2++){
			i2 = 2*j2-d2;
			ku2 = i2*m1k;
			ks2 = j2*m1j;
			for(int j1 = 1; j1 < m1j; j1++){
				i1 = 2*j1-d1;
				x1[i1] = r[ku3P+ku2+i1] + r[ku3P+(ku2+2*m1k)+i1]
					+ r[ku3+ku2+m1k+i1] + r[ku3P2+ku2+m1k+i1];
				y1[i1] = r[ku3+ku2+i1] + r[ku3P2+ku2+i1]
					+ r[ku3+(ku2+2*m1k)+i1] + r[ku3P2+(ku2+2*m1k)+i1];
			}
			for(int j1 = 1; j1 < m1j-1; j1++){
				i1 = 2*j1-d1;
				y2 = r[ku3+ku2+i1+1] + r[ku3P2+ku2+i1+1]
					+ r[ku3+(ku2+2*m1k)+i1+1] + r[ku3P2+(ku2+2*m1k)+i1+1];
				x2 = r[ku3P+ku2+i1+1] + r[ku3P+(ku2+2*m1k)+i1+1]
					+ r[ku3+ku2+m1k+i1+1] + r[ku3P2+ku2+m1k+i1+1];
				s[ks3+ks2+j1] =
					0.5 * r[ku3P+ku2+m1k+i1+1]
					+ 0.25 * ( r[ku3P+ku2+m1k+i1] + r[ku3P+ku2+m1k+i1+2] + x2)
					+ 0.125 * ( x1[i1] + x1[i1+2] + y2)
					+ 0.0625 * ( y1[i1] + y1[i1+2] );
			}
		}
	});
	if(timeron){timer_stop(T_RPRJ3);}

	j=k-1;
	comm3(s,m1j,m2j,m3j);

	if(debug_vec[0] >= 1){
		rep_nrm(s,m1j,m2j,m3j,"   rprj3",k-1);	
	}

	if(debug_vec[4] >= k){
		showall(s,m1j,m2j,m3j);
	}
}

static void setup(int* n1, int* n2, int* n3, int k){
	std::vector<int> mi((MAXLEVEL+1)*3), ng((MAXLEVEL+1)*3);

	ng[lt*3] = nx[lt];
	ng[lt*3+1] = ny[lt];
	ng[lt*3+2] = nz[lt];
	std::for_each_n(policy, iter.front(), 3, [&] (int ax){
		for(int k = lt-1; k >= 1; k--){
			ng[k*3+ax] = ng[(k+1)*3+ax]/2;
		}
	});
	std::for_each_n(policy, iter.rtail()-lt-1, lt, [&] (int k){
		nx[k] = ng[k*3];
		ny[k] = ng[k*3+1];
		nz[k] = ng[k*3+2];
	});
	
	std::for_each_n(policy, iter.rtail()-lt-1, lt, [&] (int k){
		for (int ax = 0; ax < 3; ax++){
			mi[k*3+ax] = 2 + ng[k*3+ax];
		}

		m1[k] = mi[k*3];
		m2[k] = mi[k*3+1];
		m3[k] = mi[k*3+2];
	});

	k = lt;
	is1 = 2 + ng[k*3] - ng[lt*3];
	ie1 = 1 + ng[k*3];
	*n1 = 3 + ie1 - is1;
	is2 = 2 + ng[k*3+1] - ng[lt*3+1];
	ie2 = 1 + ng[k*3+1];
	*n2 = 3 + ie2 - is2;
	is3 = 2 + ng[k*3+2] - ng[lt*3+2];
	ie3 = 1 + ng[k*3+2];
	*n3 = 3 + ie3 - is3;

	ir[lt] = 0;
	for(int j = lt-1; j >= 1; j--){
		ir[j] = ir[j+1]+ONE*m1[j+1]*m2[j+1]*m3[j+1];
	}

	if(debug_vec[1] >= 1){
		std::cout
			<< " in setup, \n"
			<< std::endl;
		std::cout
			<< "   k  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3"
			<< std::endl;
		std::cout
			<< std::setw(4)
			<< k << lt
			<< ng[k*3] << ng[k*3+1] << ng[k*3+2]
			<< *n1 << *n2 << *n3
			<< is1 << is2 << is3
			<< ie1 << ie2 << ie3
			<< std::endl;
	}
}

static void showall(std::span<double> z, int n1, int n2, int n3){
	int m1, m2, m3;
  	int ku2, ku3;

	m1 = std::min(n1,18);
	m2 = std::min(n2,14);
	m3 = std::min(n3,18);

	std::cout << std::endl;
	for(int i3 = 0; i3 < m3; i3++){
    	ku3 = i3*n1*n2;
		for(int i2 = 0; i2 < m2; i2++){
      		ku2 = i2*n1;
			for(int i1 = 0; i1 < m1; i1++){	
				std::cout
          			<< std::setw(6) << std::setprecision(3) << z[ku3+ku2+i1];
			}
			std::cout << std::endl;
		}
		std::cout
			<< " - - - - - - - "
			<< std::endl;
	}
	std::cout << std::endl;
}

/*
 * ---------------------------------------------------------------------
 * zran3 loads +1 at ten randomly chosen points,
 * loads -1 at a different ten random points,
 * and zero elsewhere.
 * ---------------------------------------------------------------------
 */
static void zran3(std::vector<double> &z, int n1, int n2, int n3, int nx, int ny, int k){
	int i0, m0, m1;
  	int ku2, ku3;

	int i, i1, d1, e1, e2, e3;
	double xx, x0, x1, a1, a2, ai, best;

	std::vector<double> ten(2*MM);
	std::vector<int> j1(2*MM), j2(2*MM), j3(2*MM), jg(2*MM*4);

	a1 = power(A, nx);
	a2 = power(A, nx*ny);

	std::fill(z.begin(), z.end(), 0.0);

	i = is1-2+nx*(is2-2+ny*(is3-2));

	ai = power(A, i);
	d1 = ie1 - is1 + 1;
	e1 = ie1 - is1 + 2;
	e2 = ie2 - is2 + 2;
	e3 = ie3 - is3 + 2;
	x0 = X;
	randlc(x0, ai);
	for(int i3 = 1; i3 < e3; i3++){
		x1 = x0;
		for(int i2 = 1; i2 < e2; i2++){
			xx = x1;
			vranlc(d1, xx, A, &(z[(i3*n2+i2)*n1+1]));
			randlc(x1,a1);
		}
		randlc(x0, a2);
	}

	/*
	 * ---------------------------------------------------------------------
	 * each processor looks for twenty candidates
	 * ---------------------------------------------------------------------
	 */	
	for(int i = 0; i < MM; i++){
		ten[MM+i] = 0.0;
		j1[MM+i] = 0;
		j2[MM+i] = 0;
		j3[MM+i] = 0;
		ten[i] = 1.0;
		j1[i] = 0;
		j2[i] = 0;
		j3[i] = 0;
	}

	for(int i3 = 1; i3 < n3-1; i3++){
    	ku3 = i3*n1*n2;
		for(int i2 = 1; i2 < n2-1; i2++){
      		ku2 = i2*n1;
			for(int i1 = 1; i1 < n1-1; i1++){
				if(z[ku3+ku2+i1] > ten[1*MM]){
					ten[1*MM] = z[ku3+ku2+i1];
					j1[MM] = i1;
					j2[MM] = i2;
					j3[MM] = i3;
					bubble(ten, j1, j2, j3, MM, 1);
				}
				if(z[ku3+ku2+i1] < ten[0]){
					ten[0] = z[ku3+ku2+i1];
					j1[0] = i1;
					j2[0] = i2;
					j3[0] = i3;
					bubble(ten, j1, j2, j3, MM, 0);
				}
			}
		}
	}

	/*
	 * ---------------------------------------------------------------------
	 * now which of these are globally best?
	 * ---------------------------------------------------------------------
	 */	
	i1 = MM - 1;
	i0 = MM - 1;
	for(int i = MM - 1; i >= 0; i--){
		best = 0.0;
		int jg_const=(MM+i)*4;
		jg[jg_const] = 0;
		if(best < ten[1*MM+i1]){
			jg[jg_const] = 0;
			jg[jg_const+1] = is1 - 2 + j1[MM+i1];
			jg[jg_const+2] = is2 - 2 + j2[MM+i1];
			jg[jg_const+3] = is3 - 2 + j3[MM+i1];
			i1 = i1-1;
		}else{
			jg[jg_const] = 0;
			jg[jg_const+1] = 0;
			jg[jg_const+2] = 0;
			jg[jg_const+3] = 0;
		}
		best = 1.0;
		jg_const=i*4;
		if(best > ten[i0]){
			jg[jg_const] = 0;
			jg[jg_const+1] = is1 - 2 + j1[i0];
			jg[jg_const+2] = is2 - 2 + j2[i0];
			jg[jg_const+3] = is3 - 2 + j3[i0];
			i0 = i0-1;
		}else{
			jg[jg_const] = 0;
			jg[jg_const+1] = 0;
			jg[jg_const+2] = 0;
			jg[jg_const+3] = 0;
		}
	}
	m1 = 0;
	m0 = 0;

	std::fill(z.begin(), z.end(), 0.0);
	for (int i = MM-1; i >= m0; i--){
		int jg_const=(MM+i)*4;
		z[   jg[jg_const+3]*n2*n1
			+jg[jg_const+2]*n1
			+jg[jg_const+1]] = -1.0;
	}
	for(int i = MM-1; i >= m1; i--){
		int jg_const=i*4;
		z[   jg[jg_const+3]*n2*n1
			+jg[jg_const+2]*n1
			+jg[jg_const+1]] = +1.0;
	}
	comm3(z, n1, n2, n3);
}
