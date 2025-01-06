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
	M. Yarrow
	C. Kuszmaul

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

/*
 * ---------------------------------------------------------------------
 * note: please observe that in the routine conj_grad three 
 * implementations of the sparse matrix-vector multiply have
 * been supplied. the default matrix-vector multiply is not
 * loop unrolled. the alternate implementations are unrolled
 * to a depth of 2 and unrolled to a depth of 8. please
 * experiment with these to find the fastest for your particular
 * architecture. if reporting timing results, any of these three may
 * be used without penalty.
 * ---------------------------------------------------------------------
 * class specific parameters: 
 * it appears here for reference only.
 * these are their values, however, this info is imported in the npbparams.h
 * include file, which is written by the sys/setparams.c program.
 * ---------------------------------------------------------------------
 */
#define NZ (NA*(NONZER+1)*(NONZER+1))
#define NAZ (NA*(NONZER+1))
#define T_INIT 0
#define T_BENCH 1
#define T_CONJ_GRAD 2
#define T_LAST 3

/* global variables */
std::vector<int> acol(NAZ, 0);
std::vector<double> aelt(NAZ, 0);
std::vector<int> colidx(NZ, 0);
std::vector<int> rowstr(NA+1, 0);
std::vector<int> iv(NA, 0);
std::vector<int> arow(NA, 0);
std::vector<double> a(NZ, 0);
std::vector<double> x(NA+2, 1);
std::vector<double> z(NA+2, 0);
std::vector<double> p(NA+2, 0);
std::vector<double> q(NA+2, 0);
std::vector<double> r(NA+2, 0);
/* Auxiliar vars for index calc */
int kx2_const = NONZER+1;
int kx2;
/* Auxiliar vars for index calc */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static double amult;
static double tran;
static boolean timeron;

auto policy = std::execution::par;
CountIterator iter(NZ); 

/* function prototypes */
static void conj_grad(std::vector<int> &colidx,
		std::vector<int> &rowstr,
		std::vector<double> &x,
		std::vector<double> &z,
		std::vector<double> &a,
		std::vector<double> &p,
		std::vector<double> &q,
		std::vector<double> &r,
		double &rnorm);
static int icnvrt(double x,
		int ipwr2);
static void makea(int n,
		int nz,
		std::vector<double> &a,
		std::vector<int> &colidx,
		std::vector<int> &rowstr,
		int firstrow,
		int lastrow,
		int firstcol,
		int lastcol,
		std::vector<int> &arow,
		std::vector<int> &acol,
		std::vector<double> &aelt,
		std::vector<int> &iv);
static void sparse(std::vector<double> &a,
		std::vector<int> &colidx,
		std::vector<int> &rowstr,
		int n,
		int nz,
		int nozer,
		std::vector<int> &arow,
		std::vector<int> &acol,
		std::vector<double> &aelt,
		int firstrow,
		int lastrow,
		std::vector<int> &nzloc,
		double rcond,
		double shift);
static void sprnvc(int n,
		int nz,
		int nn1,
		std::vector<double> &v,
		std::vector<int> &iv);
static void vecset(int n,
		std::vector<double> &v,
		std::vector<int> &iv,
		int &nzv,
		int i,
		double val);

/* cg */
int main(int argc, char **argv){
	double zeta;
	double rnorm;
	double norm_temp1, norm_temp2;
	double t, mflops, tmax;
	char class_npb;
	boolean verified;
	double zeta_verify_value, epsilon, err;

	std::vector<std::string> t_names(T_LAST);

	for(int i=0; i<T_LAST; i++){
		timer_clear(i);
	}

	timeron = std::filesystem::exists("timer.flag");
	if(timeron){
		t_names[T_INIT] = "init";
		t_names[T_BENCH] = "benchmk";
		t_names[T_CONJ_GRAD] = "conjgd";
	}

	timer_start(T_INIT);

	firstrow = 0;
	lastrow  = NA-1;
	firstcol = 0;
	lastcol  = NA-1;

	if(NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0){
		class_npb = 'S';
		zeta_verify_value = 8.5971775078648;
	}else if(NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0){
		class_npb = 'W';
		zeta_verify_value = 10.362595087124;
	}else if(NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0){
		class_npb = 'A';
		zeta_verify_value = 17.130235054029;
	}else if(NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0){
		class_npb = 'B';
		zeta_verify_value = 22.712745482631;
	}else if(NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0){
		class_npb = 'C';
		zeta_verify_value = 28.973605592845;
	}else if(NA == 1500000 && NONZER == 21 && NITER == 100 && SHIFT == 500.0){
		class_npb = 'D';
		zeta_verify_value = 52.514532105794;
	}else if(NA == 9000000 && NONZER == 26 && NITER == 100 && SHIFT == 1500.0){
		class_npb = 'E';
		zeta_verify_value = 77.522164599383;
	}else{
		class_npb = 'U';
	}

	std::cout
		<< std::endl
		<< std::endl;
	std::cout << " NAS Parallel Benchmarks 4.1 Serial C++ version - CG Benchmark";
	std::cout
		<< std::endl
		<< std::endl;
	std::cout
		<< " Size: "
		<< std::setw(11) << NA
		<< std::endl;
	std::cout
		<< " Iterations: "
		<< std::setw(5) << NITER
		<< std::endl;

	naa = NA;
	nzz = NZ;

	/* initialize random number generator */
	tran    = 314159265.0;
	amult   = 1220703125.0;
	zeta    = randlc( tran, amult );

	makea(naa, 
			nzz, 
			a, 
			colidx, 
			rowstr, 
			firstrow, 
			lastrow, 
			firstcol, 
			lastcol, 
			arow, 
			acol, 
			aelt,
			iv);

	/*
	 * ---------------------------------------------------------------------
	 * note: as a result of the above call to makea:
	 * values of j used in indexing rowstr go from 0 --> lastrow-firstrow
	 * values of colidx which are col indexes go from firstcol --> lastcol
	 * so:
	 * shift the col index vals from actual (firstcol --> lastcol) 
	 * to local, i.e., (0 --> lastcol-firstcol)
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front(), (lastrow - firstrow + 1), [&] (int j) {
		std::transform(colidx.begin()+rowstr[j], colidx.begin()+rowstr[j+1], colidx.begin()+rowstr[j], [&] (int colidx_el) {
			return colidx_el - firstcol;
		});
	});

	zeta = 0.0;

	/*
	 * -------------------------------------------------------------------
	 * ---->
	 * do one iteration untimed to init all code and data page tables
	 * ----> (then reinit, start timing, to niter its)
	 * -------------------------------------------------------------------*/
	for(int it = 1; it <= 1; it++){
		/* the call to the conjugate gradient routine */
		conj_grad(colidx, rowstr, x, z, a, p, q, r, rnorm);

		/*
		 * --------------------------------------------------------------------
		 * zeta = shift + 1/(x.z)
		 * so, first: (x.z)
		 * also, find norm of z
		 * so, first: (z.z)
		 * --------------------------------------------------------------------
		 */
		norm_temp1 = 0.0;
		norm_temp2 = 0.0;

		norm_temp1 = std::transform_reduce(policy,
				x.cbegin(),
				x.cbegin()+(lastcol - firstcol + 1),
				z.cbegin(),
				0.0,
				std::plus<double>(),
				[] (double x_el, double z_el) {
					return x_el*z_el;
				});

		norm_temp2 = std::transform_reduce(policy, 
				z.cbegin(), 
				z.cbegin()+(lastcol-firstcol+1), 
				z.cbegin(), 
				0.0, 
				std::plus<double>(), 
				std::multiplies<double>());

		norm_temp2 = 1.0 / sqrt(norm_temp2);

		/* normalize z to obtain x */
		std::transform(z.begin(), z.begin()+(lastcol - firstcol + 1), x.begin(), [norm_temp2] (double z_el) {
			return norm_temp2 * z_el;
		});
	} /* end of do one iteration untimed */

	/* set starting vector to (1, 1, .... 1) */	
	std::fill(x.begin(), x.begin()+(NA+1), 1.0);

	zeta = 0.0;

	timer_stop(T_INIT);

	std::cout
		<< " Initialization time = "
		<< std::setw(15) << std::setprecision(3) << timer_read(T_INIT)
		<< " seconds"
		<< std::endl;

	timer_start(T_BENCH);

	/*
	 * --------------------------------------------------------------------
	 * ---->
	 * main iteration for inverse power method
	 * ---->
	 * --------------------------------------------------------------------
	 */
	for(int it = 1; it <= NITER; it++){
		/* the call to the conjugate gradient routine */
		if(timeron) timer_start(T_CONJ_GRAD);
		conj_grad(colidx, rowstr, x, z, a, p, q, r, rnorm);
		if(timeron) timer_stop(T_CONJ_GRAD);

		/*
		 * --------------------------------------------------------------------
		 * zeta = shift + 1/(x.z)
		 * so, first: (x.z)
		 * also, find norm of z
		 * so, first: (z.z)
		 * --------------------------------------------------------------------
		 */
		norm_temp1 = 0.0;
		norm_temp2 = 0.0;

		norm_temp1 = std::transform_reduce(policy,
				x.cbegin(),
				x.cbegin()+(lastcol - firstcol + 1),
				z.cbegin(),
				0.0,
				std::plus<double>(),
				[] (double x_el, double z_el) {
					return x_el*z_el;
				});

		norm_temp2 = std::transform_reduce(policy, 
				z.cbegin(), 
				z.cbegin()+(lastcol-firstcol+1), 
				z.cbegin(), 
				0.0, 
				std::plus<double>(), 
				std::multiplies<double>());

		norm_temp2 = 1.0 / sqrt(norm_temp2);
		zeta = SHIFT + 1.0 / norm_temp1;
		if(it==1){
			std::cout
				<< std::endl
				<< "   iteration           ||r||                 zeta"
				<< std::endl;
		}
		std::cout
			<< "    "
			<< std::setw(5) << it
			<< "       "
			<< std::setw(20) << std::setprecision(14) << std::scientific << rnorm
			<< std::setw(20) << std::setprecision(13) << std::scientific << zeta
			<< std::endl;

		/* normalize z to obtain x */
		std::transform(policy, z.begin(), z.begin()+(lastcol - firstcol + 1), x.begin(), [norm_temp2] (double z_el) {
			return norm_temp2 * z_el;
		});

	} /* end of main iter inv pow meth */

	timer_stop(T_BENCH);

	/*
	 * --------------------------------------------------------------------
	 * end of timed section
	 * --------------------------------------------------------------------
	 */

	t = timer_read(T_BENCH);

	std::cout
		<< " Benchmark completed"
		<< std::endl;

	epsilon = 1.0e-10;
	if(class_npb != 'U'){
		err = fabs(zeta - zeta_verify_value) / zeta_verify_value;
		if(err <= epsilon){
			verified = TRUE;
			std::cout
				<< " VERIFICATION SUCCESSFUL"
				<< std::endl;
			std::cout
				<< " Zeta is    "
				<< std::setw(20) << std::setprecision(13) << std::scientific << zeta
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
				<< " Zeta                "
				<< std::setw(20) << std::setprecision(13) << std::scientific << zeta
				<< std::endl;
			std::cout
				<< " The correct zeta is "
				<< std::setw(20) << std::setprecision(13) << std::scientific << zeta_verify_value
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
	if(t != 0.0){
		mflops = (double)(2.0*NITER*NA)
			* (3.0+(double)(NONZER*(NONZER+1))
					+ 25.0
					* (5.0+(double)(NONZER*(NONZER+1)))+3.0)
			/ t / 1000000.0;
	}else{
		mflops = 0.0;
	}
	c_print_results("CG",
			class_npb,
			NA,
			0,
			0,
			NITER,
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
		if(tmax == 0.0){tmax = 1.0;}
		std::cout
			<< "  SECTION   Time (secs)"
			<< std::endl;
		for(int i = 0; i < T_LAST; i++){
			t = timer_read(i);
			if(i == T_INIT){
				std::cout
					<< "  " << std::setw(8) << t_names[i] << ":"
					<< std::setw(9) << std::setprecision(3) << t
					<< std::endl;
			}else{
				std::cout
					<< "  " << std::setw(8) << t_names[i] << ":"
					<< std::setw(9) << std::setprecision(3) << t
					<< "  (" << std::setw(6) << std::setprecision(2) << t*100.0/tmax << "%)"
					<< std::endl;
				if(i == T_CONJ_GRAD){
					t = tmax - t;
					std::cout
						<< "    --> "
						<< std::setw(8) << "rest" << ":"
						<< std::setw(9) << std::setprecision(3) << t
						<< "  (" << std::setw(6) << std::setprecision(2) << t*100.0/tmax << "%)"
						<< std::endl;
				}
			}
		}
	}

	return 0;
}

/*
 * ---------------------------------------------------------------------
 * floating point arrays here are named as in NPB1 spec discussion of 
 * CG algorithm
 * ---------------------------------------------------------------------
 */
static void conj_grad(std::vector<int> &colidx,
		std::vector<int> &rowstr,
		std::vector<double> &x,
		std::vector<double> &z,
		std::vector<double> &a,
		std::vector<double> &p,
		std::vector<double> &q,
		std::vector<double> &r,
		double &rnorm){
	int cgitmax;
	double d, sum, rho, rho0, alpha, beta;

	cgitmax = 25;

	rho = 0.0;

	/* initialize the CG algorithm */
	std::fill(policy, q.begin(), q.end(), 0);
	std::fill(policy, z.begin(), z.end(), 0);
	std::copy(policy, x.begin(), x.end(), r.begin());
	std::copy(policy, r.begin(), r.end(), p.begin());

	/*
	 * --------------------------------------------------------------------
	 * rho = r.r
	 * now, obtain the norm of r: First, sum squares of r elements locally...
	 * --------------------------------------------------------------------
	 */

	rho = std::transform_reduce(policy, r.begin(), r.begin()+(lastcol-firstcol+1), r.begin(), 0.0, std::plus<double>(), std::multiplies<double>());



	/* the conj grad iteration loop */
	for(int cgit = 1; cgit <= cgitmax; cgit++){
		/*
		 * ---------------------------------------------------------------------
		 * q = A.p
		 * the partition submatrix-vector multiply: use workspace w
		 * ---------------------------------------------------------------------
		 * 
		 * note: this version of the multiply is actually (slightly: maybe %5) 
		 * faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
		 * below. on the Cray t3d, the reverse is TRUE, i.e., the 
		 * unrolled-by-two version is some 10% faster.  
		 * the unrolled-by-8 version below is significantly faster
		 * on the Cray t3d - overall speed of code is 1.5 times faster.
		 */

		std::for_each(policy,
			iter.front(),
			iter.front()+(lastrow - firstrow + 1),
			[&] (int index) {
				q[index] = std::transform_reduce(
					a.cbegin()+rowstr[index],
					a.cbegin()+rowstr[index+1],
					colidx.cbegin()+rowstr[index],
					0.0,
					std::plus<double>(),
					[&p] (double a_el, int colidx_el) {
						return a_el*p[colidx_el];
					});
			});

		/*
		 * --------------------------------------------------------------------
		 * obtain p.q
		 * --------------------------------------------------------------------
		 */

		d = std::transform_reduce(policy,
				p.cbegin(),
				p.cbegin()+(lastcol - firstcol + 1),
				q.cbegin(),
				0.0,
				std::plus<double>(),
				[] (double p_el, double q_el) {
					return p_el*q_el;
				});

		/*
		 * --------------------------------------------------------------------
		 * obtain alpha = rho / (p.q)
		 * -------------------------------------------------------------------
		 */
		alpha = rho / d;

		/*
		 * --------------------------------------------------------------------
		 * save a temporary of rho
		 * --------------------------------------------------------------------
		 */
		rho0 = rho;

		/*
		 * ---------------------------------------------------------------------
		 * obtain z = z + alpha*p
		 * and    r = r - alpha*q
		 * ---------------------------------------------------------------------
		 */
		rho = 0.0;

		std::transform(policy,
			p.cbegin(),
			p.cbegin()+(lastcol - firstcol + 1),
			z.begin(),
			z.begin(),
			[alpha] (const double p_el, double z_el) {
				return z_el + alpha*p_el;
			});
		
		std::transform(policy,
			q.cbegin(),
			q.cbegin()+(lastcol - firstcol + 1),
			r.begin(),
			r.begin(),
			[alpha] (const double q_el, double r_el) {
				return r_el - alpha*q_el;
			});


		/*
		 * ---------------------------------------------------------------------
		 * rho = r.r
		 * now, obtain the norm of r: first, sum squares of r elements locally...
		 * ---------------------------------------------------------------------
		 */
		rho = std::transform_reduce(policy, r.begin(), r.begin()+(lastcol-firstcol+1), r.begin(), 0.0, std::plus<double>(), std::multiplies<double>());
		

		/*
		 * ---------------------------------------------------------------------
		 * obtain beta
		 * ---------------------------------------------------------------------
		 */
		beta = rho / rho0;

		/*
		 * ---------------------------------------------------------------------
		 * p = r + beta*p
		 * ---------------------------------------------------------------------
		 */

		std::transform(policy,
			r.cbegin(),
			r.cbegin()+(lastcol - firstcol + 1),
			p.cbegin(),
			p.begin(),
			[beta] (const double r_el, const double p_el) {
				return r_el + beta*p_el;
			});
	} /* end of do cgit=1, cgitmax */

	/*
	 * ---------------------------------------------------------------------
	 * compute residual norm explicitly: ||r|| = ||x - A.z||
	 * first, form A.z
	 * the partition submatrix-vector multiply
	 * ---------------------------------------------------------------------
	 */

	std::for_each(policy,
		iter.front(),
		iter.front()+(lastrow - firstrow + 1),
		[&] (int index) {
			r[index] = std::transform_reduce(
				a.cbegin()+rowstr[index],
				a.cbegin()+rowstr[index+1],
				colidx.cbegin()+rowstr[index],
				0.0,
				std::plus<double>(),
				[&z] (double a_el, int colidx_el) {
					return a_el*z[colidx_el];
				});
		});

	sum = 0.0;

	/*
	 * ---------------------------------------------------------------------
	 * at this point, r contains A.z
	 * ---------------------------------------------------------------------
	 */

	sum = std::transform_reduce(policy,
		x.cbegin(),
		x.cbegin()+(lastcol - firstcol + 1),
		r.cbegin(),
		0.0,
		std::plus<double>(),
		[] (const double x_el, const double r_el) {
			double d = x_el - r_el;
			return d*d;
		});

	rnorm = sqrt(sum);
}

/*
 * ---------------------------------------------------------------------
 * scale a double precision number x in (0,1) by a power of 2 and chop it
 * ---------------------------------------------------------------------
 */
static int icnvrt(double x, int ipwr2){
	return (int)(ipwr2 * x);
}

/*
 * ---------------------------------------------------------------------
 * generate the test problem for benchmark 6
 * makea generates a sparse matrix with a
 * prescribed sparsity distribution
 *
 * parameter    type        usage
 *
 * input
 *
 * n            i           number of cols/rows of matrix
 * nz           i           nonzeros as declared array size
 * rcond        r*8         condition number
 * shift        r*8         main diagonal shift
 *
 * output
 *
 * a            r*8         array for nonzeros
 * colidx       i           col indices
 * rowstr       i           row pointers
 *
 * workspace
 *
 * iv, arow, acol i
 * aelt           r*8
 * ---------------------------------------------------------------------
 */
static void makea(int n,
		int nz,
		std::vector<double> &a,
		std::vector<int> &colidx,
		std::vector<int> &rowstr,
		int firstrow,
		int lastrow,
		int firstcol,
		int lastcol,
		std::vector<int> &arow,
		std::vector<int> &acol,
		std::vector<double> &aelt,
		std::vector<int> &iv){
	int nzv, nn1;
	std::vector<int> ivc(NONZER+1);
	std::vector<double> vc(NONZER+1);

	/*
	 * --------------------------------------------------------------------
	 * nonzer is approximately  (int(sqrt(nnza /n)));
	 * --------------------------------------------------------------------
	 * nn1 is the smallest power of two not less than n
	 * --------------------------------------------------------------------
	 */
	nn1 = 1;
	do{
		nn1 *= 2;
	}while(nn1 < n);

	/*
	 * -------------------------------------------------------------------
	 * generate nonzero positions and save for the use in sparse
	 * -------------------------------------------------------------------
	 */
	for(int iouter = 0; iouter < n; iouter++){
		nzv = NONZER;
		sprnvc(n, nzv, nn1, vc, ivc);
		vecset(n, vc, ivc, nzv, iouter+1, 0.5);
		arow[iouter] = nzv;
		kx2 = iouter*kx2_const;
		for(int ivelt = 0; ivelt < nzv; ivelt++){
			acol[kx2+ivelt] = ivc[ivelt] - 1;
			aelt[kx2+ivelt] = vc[ivelt];
		}
	}

	/*
	 * ---------------------------------------------------------------------
	 * ... make the sparse matrix from list of elements with duplicates
	 * (iv is used as  workspace)
	 * ---------------------------------------------------------------------
	 */
	sparse(a,
			colidx,
			rowstr,
			n,
			nz,
			NONZER,
			arow,
			acol,
			aelt,
			firstrow,
			lastrow,
			iv,
			RCOND,
			SHIFT);
}

/*
 * ---------------------------------------------------------------------
 * rows range from firstrow to lastrow
 * the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
 * ---------------------------------------------------------------------
 */
static void sparse(std::vector<double> &a,
		std::vector<int> &colidx,
		std::vector<int> &rowstr,
		int n,
		int nz,
		int nozer,
		std::vector<int> &arow,
		std::vector<int> &acol,
		std::vector<double> &aelt,
		int firstrow,
		int lastrow,
		std::vector<int> &nzloc,
		double rcond,
		double shift){	
	int nrows;
	/*
	 * ---------------------------------------------------
	 * generate a sparse matrix from a list of
	 * [col, row, element] tri
	 * ---------------------------------------------------
	 */
	int i, j, j1, j2, nza, k, kk, nzrow, jcol;
	double size, scale, ratio, va;
	boolean goto_40;

	/*
	 * --------------------------------------------------------------------
	 * how many rows of result
	 * --------------------------------------------------------------------
	 */
	nrows = lastrow - firstrow + 1;

	/*
	 * --------------------------------------------------------------------
	 * ...count the number of triples in each row
	 * --------------------------------------------------------------------
	 */
	for(int j = 0; j < nrows+1; j++){
		rowstr[j] = 0;
	}
	for(int i = 0; i < n; i++){
		kx2 = i*kx2_const;
		for(int nza = 0; nza < arow[i]; nza++){
			j = acol[kx2+nza] + 1;
			rowstr[j] += arow[i];
		}
	}
	rowstr[0] = 0;
	for(int j = 1; j < nrows+1; j++){
		rowstr[j] += rowstr[j-1];
	} /* Can be made using transform */
	nza = rowstr[nrows] - 1;

	/*
	 * ---------------------------------------------------------------------
	 * ... rowstr(j) now is the location of the first nonzero
	 * of row j of a
	 * ---------------------------------------------------------------------
	 */
	if(nza > nz){
		std::cout
			<< "Space for matrix elements exceeded in sparse"
			<< std::endl;
		std::cout
			<< "nza, nzmax = " << nza << ", " << nz
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	/*
	 * ---------------------------------------------------------------------
	 * ... preload data pages
	 * ---------------------------------------------------------------------
	 */
	for(int j = 0; j < nrows; j++){
		std::fill(a.begin()+rowstr[j], a.begin()+rowstr[j+1], 0.0);
		std::fill(colidx.begin()+rowstr[j], colidx.begin()+rowstr[j+1], -1);
		nzloc[j] = 0;
	}

	/*
	 * ---------------------------------------------------------------------
	 * ... generate actual values by summing duplicates
	 * ---------------------------------------------------------------------
	 */
	size = 1.0;
	ratio = pow(rcond, (1.0 / (double)(n)));
	int pos2 = 0;
	for(int i = 0; i < n; i++){
		kx2 = i*kx2_const;
		for(int nza = 0; nza < arow[i]; nza++){
			j = acol[kx2+nza];

			scale = size * aelt[kx2+nza];
			for(int nzrow = 0; nzrow < arow[i]; nzrow++){
				jcol = acol[kx2+nzrow];
				va = aelt[kx2+nzrow] * scale;

				/*
				 * --------------------------------------------------------------------
				 * ... add the identity * rcond to the generated matrix to bound
				 * the smallest eigenvalue from below by rcond
				 * --------------------------------------------------------------------
				 */
				if(jcol == j && j == i){
					va += rcond - shift;
				}

				goto_40 = FALSE;

				for(k = rowstr[j]; k < rowstr[j+1]; k++){
					if(colidx[k] > jcol){
						/*
						 * ----------------------------------------------------------------
						 * ... insert colidx here orderly
						 * ----------------------------------------------------------------
						 */
						for(int kk = rowstr[j+1]-2; kk >= k; kk--){
							if(colidx[kk] > -1){
								a[kk+1] = a[kk];
								colidx[kk+1] = colidx[kk];
							}
						}

						colidx[k] = jcol;
						a[k]  = 0.0;
						goto_40 = TRUE;
						break;
					}else if(colidx[k] == -1){
						colidx[k] = jcol;
						goto_40 = TRUE;
						break;
					}else if(colidx[k] == jcol){
						/*
						 * --------------------------------------------------------------
						 * ... mark the duplicated entry
						 * -------------------------------------------------------------
						 */
						nzloc[j] = nzloc[j] + 1;
						goto_40 = TRUE;
						break;
					}
				}
				if(goto_40 == FALSE){
					std::cout
						<< "internal error in sparse: i=" << i
						<< std::endl;
					exit(EXIT_FAILURE);
				}
				a[k] += va;
			}
		}
		size *= ratio;
	}

	/*
	 * ---------------------------------------------------------------------
	 * ... remove empty entries and generate final results
	 * ---------------------------------------------------------------------
	 */
	for(int j = 1; j < nrows; j++){nzloc[j] += nzloc[j-1];}

	for(int j = 0; j < nrows; j++){
		if(j > 0){
			j1 = rowstr[j] - nzloc[j-1];
		}else{
			j1 = 0;
		}
		j2 = rowstr[j+1] - nzloc[j];
		nza = rowstr[j];
		for(int k = j1; k < j2; k++){
			a[k] = a[nza];
			colidx[k] = colidx[nza];
			nza += 1;
		}
	}
	for(int j = 1; j < nrows+1; j++){
		rowstr[j] -= nzloc[j-1];
	}
	nza = rowstr[nrows] - 1;
}

/*
 * ---------------------------------------------------------------------
 * generate a sparse n-vector (v, iv)
 * having nzv nonzeros
 *
 * mark(i) is set to 1 if position i is nonzero.
 * mark is all zero on entry and is reset to all zero before exit
 * this corrects a performance bug found by John G. Lewis, caused by
 * reinitialization of mark on every one of the n calls to sprnvc
 * ---------------------------------------------------------------------
 */
static void sprnvc(int n, int nz, int nn1, std::vector<double> &v, std::vector<int> &iv){
	int nzv, i;
	double vecelt, vecloc;

	nzv = 0;

	while(nzv < nz){
		vecelt = randlc(tran, amult);

		/*
		 * --------------------------------------------------------------------
		 * generate an integer between 1 and n in a portable manner
		 * --------------------------------------------------------------------
		 */
		vecloc = randlc(tran, amult);
		i = icnvrt(vecloc, nn1) + 1;
		if(i>n){continue;}

		/*
		 * --------------------------------------------------------------------
		 * was this integer generated already?
		 * --------------------------------------------------------------------
		 */
		boolean was_gen = FALSE;

		auto resFind = std::find_if(iv.begin(), iv.begin()+nzv, [i] (int el) { return el == i;});

		if (resFind != iv.begin()+nzv) {was_gen = TRUE;continue;}

		v[nzv] = vecelt;
		iv[nzv] = i;
		nzv += 1;
	}
}

/*
 * --------------------------------------------------------------------
 * set ith element of sparse vector (v, iv) with
 * nzv nonzeros to val
 * --------------------------------------------------------------------
 */
static void vecset(int n, std::vector<double> &v, std::vector<int> &iv, int &nzv, int i, double val){
	int k;
	boolean set;

	set = FALSE;
	for(int k = 0; k < nzv; k++){
		if(iv[k] == i){
			v[k] = val;
			set  = TRUE;
		}
	}

	if(set == FALSE){
		v[nzv]  = val;
		iv[nzv] = i;
		nzv     += 1;
	}
}
