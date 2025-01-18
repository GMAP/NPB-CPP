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
	S. Weeratunga
	V. Venkatakrishnan
	E. Barszcz
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

/*
 * ---------------------------------------------------------------------
 * driver for the performance evaluation of the solver for
 * five coupled parabolic/elliptic partial differential equations
 * ---------------------------------------------------------------------
 * parameters which can be overridden in runtime config file
 * isiz1,isiz2,isiz3 give the maximum size
 * ipr = 1 to print out verbose information
 * omega = 2.0 is correct for all classes
 * tolrsd is tolerance levels for steady state residuals
 * ---------------------------------------------------------------------
 * field variables and residuals
 * to improve cache performance, second two dimensions padded by 1 
 * for even number sizes only.
 * note: corresponding array (called "v") in routines blts, buts, 
 * and l2norm are similarly padded
 * ---------------------------------------------------------------------
 */
#define IPR_DEFAULT 1
#define OMEGA_DEFAULT 1.2
#define TOLRSD1_DEF 1.0e-08
#define TOLRSD2_DEF 1.0e-08
#define TOLRSD3_DEF 1.0e-08
#define TOLRSD4_DEF 1.0e-08
#define TOLRSD5_DEF 1.0e-08
#define C1 1.40e+00
#define C2 0.40e+00
#define C3 1.00e-01
#define C4 1.00e+00
#define C5 1.40e+00
#define T_TOTAL 1
#define T_RHSX 2
#define T_RHSY 3
#define T_RHSZ 4
#define T_RHS 5
#define T_JACLD 6
#define T_BLTS 7
#define T_JACU 8
#define T_BUTS 9
#define T_ADD 10
#define T_L2NORM 11
#define T_LAST 11

/* global variables */
std::vector<double> u((ISIZ3)*(ISIZ2/2*2+1)*(ISIZ1/2*2+1)*(5));
std::vector<double> rsd((ISIZ3)*(ISIZ2/2*2+1)*(ISIZ1/2*2+1)*(5));
std::vector<double> frct((ISIZ3)*(ISIZ2/2*2+1)*(ISIZ1/2*2+1)*(5));
std::vector<double> a((ISIZ2)*(ISIZ1/2*2+1)*(5)*(5));
std::vector<double> b((ISIZ2)*(ISIZ1/2*2+1)*(5)*(5));
std::vector<double> c((ISIZ2)*(ISIZ1/2*2+1)*(5)*(5));
std::vector<double> d((ISIZ2)*(ISIZ1/2*2+1)*(5)*(5));
std::vector<double> flux((ISIZ1)*(5));
std::vector<double> qs((ISIZ3)*(ISIZ2/2*2+1)*(ISIZ1/2*2+1));
std::vector<double> rho_i((ISIZ3)*(ISIZ2/2*2+1)*(ISIZ1/2*2+1));
std::vector<double> ce((13)*(5));
CountIterator iter(ISIZ3);

/* grid */
static double dxi, deta, dzeta;
static double tx1, tx2, tx3;
static double ty1, ty2, ty3;
static double tz1, tz2, tz3;
static int nx, ny, nz;
static int nx0, ny0, nz0;
static int ist, iend;
static int jst, jend;
static int ii1, ii2;
static int ji1, ji2;
static int ki1, ki2;
/* index calculation constants */
const int ku4_const = (ISIZ2/2*2+1)*(ISIZ1/2*2+1)*5;
const int ku3_const = (ISIZ1/2*2+1)*5, ku2_const = 5;
int ku4, ku3, ku2;
int ku2P, ku2P2, ku2M, ku2M2;
int ku3P, ku3P2, ku3M, ku3M2;
int ku4P, ku4P2, ku4M, ku4M2;

const int ks4_const = (ISIZ1/2*2+1)*5*5;
const int ks3_const = 5*5, ks2_const = 5;
int ks4, ks3, ks2;
int ks2T1=5, ks2T2=10, ks2T3=15, ks2T4=20;

int kfl2, fl_const=5;
int kfl2P, kfl2M;

const int kq3_const=(ISIZ2/2*2+1)*(ISIZ1/2*2+1);
const int kq2_const=(ISIZ1/2*2+1);
int kq3, kq2;
int kq2P, kq2P2, kq2M, kq2M2;
int kq3P, kq3P2, kq3M, kq3M2; 
/* dissipation */
static double dx1, dx2, dx3, dx4, dx5;
static double dy1, dy2, dy3, dy4, dy5;
static double dz1, dz2, dz3, dz4, dz5;
static double dssp;
/* output control parameters */
static int ipr, inorm;
/* newton-raphson iteration control parameters */
static double dt, omega, frc, ttotal;
std::vector<double> tolrsd(5), rsdnm(5), errnm(5);
static int itmax, invert;
/* timer */
static double maxtime;
static boolean timeron;

auto policy = std::execution::par;
const int num_workers = std::thread::hardware_concurrency();
static std::vector<boolean> flag((ISIZ1/2*2+1)*1024);


struct Returnable{
	std::vector<double> sum;
	
	Returnable() : sum(std::vector<double>(5)) {}
	Returnable(const std::vector<double> &sum) : sum(sum) {}

	Returnable operator+(Returnable other) const {
		std::transform(sum.begin(), sum.end(), 	other.sum.begin(), other.sum.begin(), std::plus<double>());
		return other;
	}
};

#include <mutex>
#include <condition_variable>
std::mutex mutex;
std::condition_variable cond[1024];

void WAIT_SIGNAL_BUTS(int k, int worker_id) {
	std::unique_lock<std::mutex> lock(mutex);
	if (worker_id != num_workers-1) cond[worker_id].wait(lock, [&]{return flag[(k*1024)+worker_id+1] == 1;});
}

void GO_SIGNAL_BUTS(int k, int worker_id) {
	std::unique_lock<std::mutex> lock(mutex);
	if (worker_id != num_workers-1) flag[(k*1024)+worker_id+1] = 0;
	if (worker_id != 0) {
		flag[(k*1024)+worker_id] = 1;  
		cond[worker_id-1].notify_one();
	}
}

void WAIT_SIGNAL_BLTS(int k, int worker_id) {
	std::unique_lock<std::mutex> lock(mutex);
	if (worker_id != 0) cond[worker_id].wait(lock, [&]{return flag[(k*1024)+worker_id-1] == 1;});
}

void GO_SIGNAL_BLTS(int k, int worker_id) {
	std::unique_lock<std::mutex> lock(mutex);
	if (worker_id != num_workers-1){
		flag[(k*1024)+worker_id] = 1; 
		cond[worker_id+1].notify_one();
	}
	if (worker_id != 0) flag[(k*1024)+worker_id-1] = 0;
}

/* function prototypes */
void blts(int nx,
		int ny,
		int nz,
		int k,
		double omega,
		std::vector<double> &v,
		std::vector<double> &ldz,
		std::vector<double> &ldy,
		std::vector<double> &ldx,
		std::vector<double> &d,
		int ist,
		int iend,
		int jst,
		int jend,
		int nx0,
		int ny0, 
		int worker_id);
void buts(int nx,
		int ny,
		int nz,
		int k,
		double omega,
		std::vector<double> &v,
		std::vector<double> &d,
		std::vector<double> &udx,
		std::vector<double> &udy,
		std::vector<double> &udz,
		int ist,
		int iend,
		int jst,
		int jend,
		int nx0,
		int ny0, 
		int worker_id);
void domain();
void erhs();
void error();
void exact(int i,
		int j,
		int k,
		std::vector<double> &u000ijk);
void jacld(int k, int worker_id);
void jacu(int k, int worker_id);
void l2norm(int nx0,
		int ny0,
		int nz0,
		int ist,
		int iend,
		int jst,
		int jend,
		std::vector<double> &v,
		double sum[5]);
void pintgr();
void read_input();
void rhs();
void setbv();
void setcoeff();
void setiv();
void ssor(int niter);
void verify(std::vector<double> &xcr,
		std::vector<double> &xce,
		double xci,
		char* class_npb,
		boolean* verified);

/* lu */
int main(int argc, char* argv[]){
	char class_npb;
	boolean verified;
	double mflops;
	double t, tmax, trecs[T_LAST+1];
	int i;
	std::vector<std::string> t_names(T_LAST+1);
	/*
	 * ---------------------------------------------------------------------
	 * setup info for timers
	 * ---------------------------------------------------------------------
	 */
	FILE *fp;
	if((fp=fopen("timer.flag","r"))!= NULL){
		timeron=TRUE;
		t_names[T_TOTAL]="total";
		t_names[T_RHSX]="rhsx";
		t_names[T_RHSY]="rhsy";
		t_names[T_RHSZ]="rhsz";
		t_names[T_RHS]="rhs";
		t_names[T_JACLD]="jacld";
		t_names[T_BLTS]="blts";
		t_names[T_JACU]="jacu";
		t_names[T_BUTS]="buts";
		t_names[T_ADD]="add";
		t_names[T_L2NORM]="l2norm";
		fclose(fp);
	}else{
		timeron=FALSE;
	}
	/*
	 * ---------------------------------------------------------------------
	 * read input data
	 * ---------------------------------------------------------------------
	 */
	read_input();
	/*
	 * ---------------------------------------------------------------------
	 * set up domain sizes
	 * ---------------------------------------------------------------------
	 */
	domain();
	/*
	 * ---------------------------------------------------------------------
	 * set up coefficients
	 * ---------------------------------------------------------------------
	 */
	setcoeff();
	/*
	 * ---------------------------------------------------------------------
	 * set the boundary values for dependent variables
	 * ---------------------------------------------------------------------
	 */
	setbv();
	/*
	 * ---------------------------------------------------------------------
	 * set the initial values for dependent variables
	 * ---------------------------------------------------------------------
	 */
	setiv();
	/*
	 * ---------------------------------------------------------------------
	 * compute the forcing term based on prescribed exact solution
	 * ---------------------------------------------------------------------
	 */
	erhs();
	/*
	 * ---------------------------------------------------------------------
	 * perform one SSOR iteration to touch all pages
	 * ---------------------------------------------------------------------
	 */
	ssor(1);
	/*
	 * ---------------------------------------------------------------------
	 * reset the boundary and initial values
	 * ---------------------------------------------------------------------
	 */
	setbv();
	setiv();
	/*
	 * ---------------------------------------------------------------------
	 * perform the SSOR iterations
	 * ---------------------------------------------------------------------
	 */
	ssor(itmax);
	/*
	 * ---------------------------------------------------------------------
	 * compute the solution error
	 * ---------------------------------------------------------------------
	 */
	error();
	/*
	 * ---------------------------------------------------------------------
	 * compute the surface integral
	 * ---------------------------------------------------------------------
	 */
	pintgr();
	/*
	 * ---------------------------------------------------------------------
	 * verification test
	 * ---------------------------------------------------------------------
	 */
	verify(rsdnm,errnm,frc,&class_npb,&verified);
	mflops=(double)itmax*(1984.77*(double)nx0
			*(double)ny0
			*(double)nz0
			-10923.3*pow(((double)(nx0+ny0+nz0)/3.0),2.0) 
			+27770.9*(double)(nx0+ny0+nz0)/3.0
			-144010.0)
		/(maxtime*1000000.0);
	c_print_results((char*)"LU",
			class_npb,
			nx0,
			ny0,
			nz0,
			itmax,
			maxtime,
			mflops,
			(char*)"          floating point",
			verified,
			(char*)NPBVERSION,
			(char*)COMPILETIME,
			(char*)COMPILERVERSION,
			(char*)CS1,
			(char*)CS2,
			(char*)CS3,
			(char*)CS4,
			(char*)CS5,
			(char*)CS6,
			(char*)"(none)");
	/*
	 * ---------------------------------------------------------------------
	 * more timers
	 * ---------------------------------------------------------------------
	 */
	if(timeron){
		for(int i=1; i<=T_LAST; i++){
			trecs[i]=timer_read(i);
		}
		if(tmax==0.0){tmax=1.0;}
    std::cout
      << "  SECTION   Time (secs)"
      << std::endl;
		for(int i=1; i<=T_LAST; i++){
			std::cout
        << "  "
        << std::left << std::setw(8) << t_names[i]
        << ":"
        << std::setw(9) << std::setprecision(3) << trecs[i]
        << "  ("
        << std::setw(6) << std::setprecision(2) << trecs[i]*100./tmax
        << "%)"
        << std::endl;
			if(i==T_RHS){
				t=trecs[T_RHSX]+trecs[T_RHSY]+trecs[T_RHSZ];
        std::cout
          << "    --> "
          << std::setw(8) << "sub-rhs"
          << ":"
          << std::setw(9) << std::setprecision(3) << t
          << "  ("
          << std::setw(6) << std::setprecision(2) << t*100./tmax
          << "%)"
          << std::endl;
				t=trecs[T_RHS]-t;
        std::cout
          << "    --> "
          << std::setw(8) << "rest-rhs"
          << ":"
          << std::setw(9) << std::setprecision(3) << t
          << "  ("
          << std::setw(6) << std::setprecision(2) << t*100./tmax
          << "%)"
          << std::endl;

			}
		}
	}
	return 0;
}

/*
 * ---------------------------------------------------------------------
 * compute the regular-sparse, block lower triangular solution:
 * v <-- ( L-inv ) * v
 * ---------------------------------------------------------------------
 * to improve cache performance, second two dimensions padded by 1 
 * for even number sizes only. only needed in v.
 * ---------------------------------------------------------------------
 */
void blts(int nx,
		int ny,
		int nz,
		int k,
		double omega,
		std::vector<double> &v, 
		std::vector<double> &ldz,
		std::vector<double> &ldy,
		std::vector<double> &ldx,
		std::vector<double> &d,
		int ist,
		int iend,
		int jst,
		int jend,
		int nx0,
		int ny0,
		int worker_id){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	int ku4, ku3, ku2, ku4M, ku3M, ku2M, ks4, ks3, ks2;
	double tmp, tmp1;
	std::vector<double> tmat(5*5);
  	std::vector<double> tv(5);

	int chunk = (jend-jst)/num_workers;
	int residual = (jend-jst)%num_workers;

	int jst_worker = worker_id*chunk + ((residual>1 && (worker_id<residual)) ? worker_id : residual) + jst;
	int jend_worker = worker_id*chunk + chunk + ((residual>1 && (worker_id<residual)) ? worker_id+1 : residual) + jst;
	
	if(worker_id==0) jst_worker = jst;
	if(worker_id==num_workers-1) jend_worker = jend;

	ku4 = k*ku4_const;
	ku4M = ku4-ku4_const;
	for(int j=jst_worker; j<jend_worker; j++){
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		for(int i=ist; i<iend; i++){
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			for(int m=0; m<5; m++){
				v[ku4+ku3+ku2+m]-=omega*(ldz[ks4+ks3+m]*v[ku4M+ku3+ku2]
							+ldz[ks4+ks3+ks2T1+m]*v[ku4M+ku3+ku2+1]
							+ldz[ks4+ks3+ks2T2+m]*v[ku4M+ku3+ku2+2]
							+ldz[ks4+ks3+ks2T3+m]*v[ku4M+ku3+ku2+3]
							+ldz[ks4+ks3+ks2T4+m]*v[ku4M+ku3+ku2+4]);
			}
		}
	}
	
	WAIT_SIGNAL_BLTS(k,worker_id);
	
	for(int j=jst_worker; j<jend_worker; j++){
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		ku3M = ku3-ku3_const;
		for(int i=ist; i<iend; i++){
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			ku2M = ku2-ku2_const;
			for(int m=0; m<5; m++){
				tv[m]=v[ku4+ku3+ku2+m]
					-omega*(ldy[ks4+ks3+m]*v[ku4+ku3M+ku2]
							+ldx[ks4+ks3+m]*v[ku4+ku3+ku2M]
							+ldy[ks4+ks3+ks2T1+m]*v[ku4+ku3M+ku2+1]
							+ldx[ks4+ks3+ks2T1+m]*v[ku4+ku3+ku2M+1]
							+ldy[ks4+ks3+ks2T2+m]*v[ku4+ku3M+ku2+2]
							+ldx[ks4+ks3+ks2T2+m]*v[ku4+ku3+ku2M+2]
							+ldy[ks4+ks3+ks2T3+m]*v[ku4+ku3M+ku2+3]
							+ldx[ks4+ks3+ks2T3+m]*v[ku4+ku3+ku2M+3]
							+ldy[ks4+ks3+ks2T4+m]*v[ku4+ku3M+ku2+4]
							+ldx[ks4+ks3+ks2T4+m]*v[ku4+ku3+ku2M+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * diagonal block inversion
			 * 
			 * forward elimination
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				tmat[m]=d[ks4+ks3+m];
				tmat[5+m]=d[ks4+ks3+ks2T1+m];
				tmat[10+m]=d[ks4+ks3+ks2T2+m];
				tmat[15+m]=d[ks4+ks3+ks2T3+m];
				tmat[20+m]=d[ks4+ks3+ks2T4+m];
			}
			/* */
			tmp1=1.0/tmat[0];
			tmp=tmp1*tmat[1];
			tmat[5+1]=tmat[5+1]-tmp*tmat[5];
			tmat[10+1]=tmat[10+1]-tmp*tmat[10];
			tmat[15+1]=tmat[15+1]-tmp*tmat[15];
			tmat[20+1]=tmat[20+1]-tmp*tmat[20];
			tv[1]=tv[1]-tv[0]*tmp;
			/* */
			tmp=tmp1*tmat[2];
			tmat[5+2]=tmat[5+2]-tmp*tmat[5];
			tmat[10+2]=tmat[10+2]-tmp*tmat[10];
			tmat[15+2]=tmat[15+2]-tmp*tmat[15];
			tmat[20+2]=tmat[20+2]-tmp*tmat[20];
			tv[2]=tv[2]-tv[0]*tmp;
			/* */
			tmp=tmp1*tmat[3];
			tmat[5+3]=tmat[5+3]-tmp*tmat[5];
			tmat[10+3]=tmat[10+3]-tmp*tmat[10];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20];
			tv[3]=tv[3]-tv[0]*tmp;
			/* */
			tmp=tmp1*tmat[4];
			tmat[5+4]=tmat[5+4]-tmp*tmat[5];
			tmat[10+4]=tmat[10+4]-tmp*tmat[10];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20];
			tv[4]=tv[4]-tv[0]*tmp;
			/* */
			tmp1=1.0/tmat[5+1];
			tmp=tmp1*tmat[5+2];
			tmat[10+2]=tmat[10+2]-tmp*tmat[10+1];
			tmat[15+2]=tmat[15+2]-tmp*tmat[15+1];
			tmat[20+2]=tmat[20+2]-tmp*tmat[20+1];
			tv[2]=tv[2]-tv[1]*tmp;
			/* */
			tmp=tmp1*tmat[5+3];
			tmat[10+3]=tmat[10+3]-tmp*tmat[10+1];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15+1];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20+1];
			tv[3]=tv[3]-tv[1]*tmp;
			/* */
			tmp=tmp1*tmat[5+4];
			tmat[10+4]=tmat[10+4]-tmp*tmat[10+1];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15+1];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+1];
			tv[4]=tv[4]-tv[1]*tmp;
			/* */
			tmp1=1.0/tmat[10+2];
			tmp=tmp1*tmat[10+3];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15+2];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20+2];
			tv[3]=tv[3]-tv[2]*tmp;
			/* */
			tmp=tmp1*tmat[10+4];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15+2];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+2];
			tv[4]=tv[4]-tv[2]*tmp;
			/* */
			tmp1=1.0/tmat[15+3];
			tmp=tmp1*tmat[15+4];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+3];
			tv[4]=tv[4]-tv[3]*tmp;
			/*
			 * ---------------------------------------------------------------------
			 * back substitution
			 * ---------------------------------------------------------------------
			 */
			v[ku4+ku3+ku2+4]=tv[4]/tmat[20+4];
			tv[3]=tv[3] 
				-tmat[20+3]*v[ku4+ku3+ku2+4];
			v[ku4+ku3+ku2+3]=tv[3]/tmat[15+3];
			tv[2]=tv[2]
				-tmat[15+2]*v[ku4+ku3+ku2+3]
				-tmat[20+2]*v[ku4+ku3+ku2+4];
			v[ku4+ku3+ku2+2]=tv[2]/tmat[10+2];
			tv[1]=tv[1]
				-tmat[10+1]*v[ku4+ku3+ku2+2]
				-tmat[15+1]*v[ku4+ku3+ku2+3]
				-tmat[20+1]*v[ku4+ku3+ku2+4];
			v[ku4+ku3+ku2+1]=tv[1]/tmat[5+1];
			tv[0]=tv[0]
				-tmat[5]*v[ku4+ku3+ku2+1]
				-tmat[10]*v[ku4+ku3+ku2+2]
				-tmat[15]*v[ku4+ku3+ku2+3]
				-tmat[20]*v[ku4+ku3+ku2+4];
			v[ku4+ku3+ku2]=tv[0]/tmat[0];
		}
	}

	GO_SIGNAL_BLTS(k,worker_id);
}

/*
 * ---------------------------------------------------------------------
 * compute the regular-sparse, block upper triangular solution:
 * v <-- ( U-inv ) * v
 * ---------------------------------------------------------------------
 * to improve cache performance, second two dimensions padded by 1 
 * for even number sizes only. only needed in v.
 * ---------------------------------------------------------------------
 */
void buts(int nx,
		int ny,
		int nz,
		int k,
		double omega,
		std::vector<double> &v,
		std::vector<double> &d,
		std::vector<double> &udx,
		std::vector<double> &udy,
		std::vector<double> &udz,
		int ist,
		int iend,
		int jst,
		int jend,
		int nx0,
		int ny0,
		int worker_id){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	int ku4, ku3, ku2, ku4P, ku3P, ku2P, ks4, ks3, ks2;
	std::vector<double> tv(ISIZ2*(ISIZ1/2*2+1)*5);
	double tmp, tmp1;
	std::vector<double> tmat(5*5);
	worker_id = num_workers-1-worker_id;
	int chunk = (jend-jst)/num_workers;
	int residual = (jend-jst)%num_workers;

	int jst_worker = worker_id*chunk + ((residual>1 && (worker_id<residual)) ? worker_id : residual) + jst-1;
	int jend_worker = worker_id*chunk + chunk + ((residual>1 && (worker_id<residual)) ? worker_id+1 : residual) + jst-1;


	if(worker_id==0) jst_worker = jst-1;
	if(worker_id==num_workers-1) jend_worker = jend-1;

	ku4 = k*ku4_const;
	ku4P = ku4+ku4_const;
	for(int j=jend_worker; j>jst_worker; j--){
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		int ju3 = (ISIZ1/2*2+1)*5*j;
		for(int i=iend-1; i>=ist; i--){
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			int iu2 = 5*i;
			for(int m=0; m<5; m++){
				tv[ju3+iu2+m]= 
					omega*(udz[ks4+ks3+m]*v[ku4P+ku3+ku2]
							+udz[ks4+ks3+ks2T1+m]*v[ku4P+ku3+ku2+1]
							+udz[ks4+ks3+ks2T2+m]*v[ku4P+ku3+ku2+2]
							+udz[ks4+ks3+ks2T3+m]*v[ku4P+ku3+ku2+3]
							+udz[ks4+ks3+ks2T4+m]*v[ku4P+ku3+ku2+4]);
			}
		}
	}
	
	WAIT_SIGNAL_BUTS(k,worker_id);

	for(int j=jend_worker; j>jst_worker; j--){
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		ku3P = ku3+ku3_const;
		int ju3 = (ISIZ1/2*2+1)*5*j;
		for(int i=iend-1; i>=ist; i--){
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			ku2P = ku2+ku2_const;
			int iu2 = 5*i;
			for(int m=0; m<5; m++){
				tv[ju3+iu2+m]=tv[ju3+iu2+m]
					+omega*(udy[ks4+ks3+m]*v[ku4+ku3P+ku2]
							+udx[ks4+ks3+m]*v[ku4+ku3+ku2P]
							+udy[ks4+ks3+ks2T1+m]*v[ku4+ku3P+ku2+1]
							+udx[ks4+ks3+ks2T1+m]*v[ku4+ku3+ku2P+1]
							+udy[ks4+ks3+ks2T2+m]*v[ku4+ku3P+ku2+2]
							+udx[ks4+ks3+ks2T2+m]*v[ku4+ku3+ku2P+2]
							+udy[ks4+ks3+ks2T3+m]*v[ku4+ku3P+ku2+3]
							+udx[ks4+ks3+ks2T3+m]*v[ku4+ku3+ku2P+3]
							+udy[ks4+ks3+ks2T4+m]*v[ku4+ku3P+ku2+4]
							+udx[ks4+ks3+ks2T4+m]*v[ku4+ku3+ku2P+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * diagonal block inversion
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				tmat[m]=d[ks4+ks3+m];
				tmat[5+m]=d[ks4+ks3+ks2T1+m];
				tmat[10+m]=d[ks4+ks3+ks2T2+m];
				tmat[15+m]=d[ks4+ks3+ks2T3+m];
				tmat[20+m]=d[ks4+ks3+ks2T4+m];
			}
			/* */
			tmp1=1.0/tmat[0];
			tmp=tmp1*tmat[1];
			tmat[5+1]=tmat[5+1]-tmp*tmat[5];
			tmat[10+1]=tmat[10+1]-tmp*tmat[10];
			tmat[15+1]=tmat[15+1]-tmp*tmat[15];
			tmat[20+1]=tmat[20+1]-tmp*tmat[20];
			tv[ju3+iu2+1]=tv[ju3+iu2+1]-tv[ju3+iu2]*tmp;
			/* */
			tmp=tmp1*tmat[2];
			tmat[5+2]=tmat[5+2]-tmp*tmat[5];
			tmat[10+2]=tmat[10+2]-tmp*tmat[10];
			tmat[15+2]=tmat[15+2]-tmp*tmat[15];
			tmat[20+2]=tmat[20+2]-tmp*tmat[20];
			tv[ju3+iu2+2]=tv[ju3+iu2+2]-tv[ju3+iu2]*tmp;
			/* */
			tmp=tmp1*tmat[3];
			tmat[5+3]=tmat[5+3]-tmp*tmat[5];
			tmat[10+3]=tmat[10+3]-tmp*tmat[10];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20];
			tv[ju3+iu2+3]=tv[ju3+iu2+3]-tv[ju3+iu2]*tmp;
			/* */
			tmp=tmp1*tmat[4];
			tmat[5+4]=tmat[5+4]-tmp*tmat[5];
			tmat[10+4]=tmat[10+4]-tmp*tmat[10];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20];
			tv[ju3+iu2+4]=tv[ju3+iu2+4]-tv[ju3+iu2]*tmp;
			/* */
			tmp1=1.0/tmat[5+1];
			tmp=tmp1*tmat[5+2];
			tmat[10+2]=tmat[10+2]-tmp*tmat[10+1];
			tmat[15+2]=tmat[15+2]-tmp*tmat[15+1];
			tmat[20+2]=tmat[20+2]-tmp*tmat[20+1];
			tv[ju3+iu2+2]=tv[ju3+iu2+2]-tv[ju3+iu2+1]*tmp;
			/* */
			tmp=tmp1*tmat[5+3];
			tmat[10+3]=tmat[10+3]-tmp*tmat[10+1];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15+1];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20+1];
			tv[ju3+iu2+3]=tv[ju3+iu2+3]-tv[ju3+iu2+1]*tmp;
			/* */
			tmp=tmp1*tmat[5+4];
			tmat[10+4]=tmat[10+4]-tmp*tmat[10+1];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15+1];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+1];
			tv[ju3+iu2+4]=tv[ju3+iu2+4]-tv[ju3+iu2+1]*tmp;
			/* */
			tmp1=1.0/tmat[10+2];
			tmp=tmp1*tmat[10+3];
			tmat[15+3]=tmat[15+3]-tmp*tmat[15+2];
			tmat[20+3]=tmat[20+3]-tmp*tmat[20+2];
			tv[ju3+iu2+3]=tv[ju3+iu2+3]-tv[ju3+iu2+2]*tmp;
			/* */
			tmp=tmp1*tmat[10+4];
			tmat[15+4]=tmat[15+4]-tmp*tmat[15+2];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+2];
			tv[ju3+iu2+4]=tv[ju3+iu2+4]-tv[ju3+iu2+2]*tmp;
			/* */
			tmp1=1.0/tmat[15+3];
			tmp=tmp1*tmat[15+4];
			tmat[20+4]=tmat[20+4]-tmp*tmat[20+3];
			tv[ju3+iu2+4]=tv[ju3+iu2+4]-tv[ju3+iu2+3]*tmp;
			/*
			 * ---------------------------------------------------------------------
			 * back substitution
			 * ---------------------------------------------------------------------
			 */
			tv[ju3+iu2+4]=tv[ju3+iu2+4]/tmat[20+4];
			tv[ju3+iu2+3]=tv[ju3+iu2+3]-tmat[20+3]*tv[ju3+iu2+4];
			tv[ju3+iu2+3]=tv[ju3+iu2+3]/tmat[15+3];
			tv[ju3+iu2+2]=tv[ju3+iu2+2]
				-tmat[15+2]*tv[ju3+iu2+3]
				-tmat[20+2]*tv[ju3+iu2+4];
			tv[ju3+iu2+2]=tv[ju3+iu2+2]/tmat[10+2];
			tv[ju3+iu2+1]=tv[ju3+iu2+1]
				-tmat[10+1]*tv[ju3+iu2+2]
				-tmat[15+1]*tv[ju3+iu2+3]
				-tmat[20+1]*tv[ju3+iu2+4];
			tv[ju3+iu2+1]=tv[ju3+iu2+1]/tmat[5+1];
			tv[ju3+iu2]=tv[ju3+iu2]
				-tmat[5]*tv[ju3+iu2+1]
				-tmat[10]*tv[ju3+iu2+2]
				-tmat[15]*tv[ju3+iu2+3]
				-tmat[20]*tv[ju3+iu2+4];
			tv[ju3+iu2]=tv[ju3+iu2]/tmat[0];
			v[ku4+ku3+ku2]-=tv[ju3+iu2];
			v[ku4+ku3+ku2+1]-=tv[ju3+iu2+1];
			v[ku4+ku3+ku2+2]-=tv[ju3+iu2+2];
			v[ku4+ku3+ku2+3]-=tv[ju3+iu2+3];
			v[ku4+ku3+ku2+4]-=tv[ju3+iu2+4];
		}
	}

	GO_SIGNAL_BUTS(k,worker_id);
}

void domain(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	nx=nx0;
	ny=ny0;
	nz=nz0;
	/*
	 * ---------------------------------------------------------------------
	 * check the sub-domain size
	 * ---------------------------------------------------------------------
	 */
	if((nx<4)||(ny<4)||(nz<4)){
		std::cout
      << "     SUBDOMAIN SIZE IS TOO SMALL - "
      << std::endl
			<< "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS"
      << std::endl
			<< "     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL"
      << std::endl
			<< "     TO 4 THEY ARE CURRENTLY"
      << std::setw(3) << nx << ny << nz
      << std::endl;
		exit(EXIT_FAILURE);
	}
	if((nx>ISIZ1)||(ny>ISIZ2)||(nz>ISIZ3)){
		std::cout
      << "     SUBDOMAIN SIZE IS TOO LARGE - "
      << std::endl
			<< "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS"
      << std::endl
			<< "     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO "
      << std::endl
			<< "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE"
      << std::endl
			<< "     CURRENTLY"
      << std::setw(4) << nx << ny << nz
      << std::endl;
		exit(EXIT_FAILURE);
	}
	/*
	 * ---------------------------------------------------------------------
	 * set up the start and end in i and j extents for all processors
	 * ---------------------------------------------------------------------
	 */
	ist=1;
	iend=nx-1;
	jst=1;
	jend=ny-1;
	ii1=1;
	ii2=nx0-1;
	ji1=1;
	ji2=ny0-2;
	ki1=2;
	ki2=nz0-1;
}

/*
 * ---------------------------------------------------------------------
 * compute the right hand side based on exact solution
 * ---------------------------------------------------------------------
 */
void erhs(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	double xi, eta, zeta;
	double q;
	double u21, u31, u41;
	double tmp;
	double u21i, u31i, u41i, u51i;
	double u21j, u31j, u41j, u51j;
	double u21k, u31k, u41k, u51k;
	double u21im1, u31im1, u41im1, u51im1;
	double u21jm1, u31jm1, u41jm1, u51jm1;
	double u21km1, u31km1, u41km1, u51km1;
	std::fill(policy, frct.begin(), frct.end(), 0.0);
	std::for_each_n(policy, iter.front(), nz, [&](int k){
		double xi, eta, zeta;
		int ku4, ku3, ku2;

		ku4 = k*ku4_const;
		zeta=((double)k)/(nz-1);
		for(int j=0; j<ny; j++){
			ku3 = j*ku3_const;
			eta=((double)j)/(ny0-1 );
			for(int i=0; i<nx; i++){
				ku2 = i*ku2_const;
				xi=((double)i)/(nx0-1);
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]=ce[(0*5)+m]+
						(ce[(1*5)+m]+
						 (ce[(4*5)+m]+
						  (ce[(7*5)+m]+
						   ce[(10*5)+m]*xi)*xi)*xi)*xi+
						(ce[(2*5)+m]+
						 (ce[(5*5)+m]+
						  (ce[(8*5)+m]+
						   ce[(11*5)+m]*eta)*eta)*eta)*eta+
						(ce[(3*5)+m]+
						 (ce[(6*5)+m]+
						  (ce[(9*5)+m]+
						   ce[(12*5)+m]*zeta)*zeta)*zeta)*zeta;
				}
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * xi-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, nz-2, [&] (int k) {
		double q;
		double u21, u31, u41;
		double tmp;
		double u21i, u31i, u41i, u51i;
		double u21im1, u31im1, u41im1, u51im1;
		int ku4, ku3, ku2;
		int kfl2, kfl2P, kfl2M;
		int ku2P, ku2P2, ku2M, ku2M2;
		std::vector<double> flux(ISIZ1*5, 0);

		ku4 = k*ku4_const;
		for(int j=jst; j<jend; j++){
			ku3 = j*ku3_const;
			for(int i=0; i<nx; i++){
				kfl2 = i*fl_const;
				ku2 = i*ku2_const;
				flux[kfl2]=rsd[ku4+ku3+ku2+1];
				u21=rsd[ku4+ku3+ku2+1]/rsd[ku4+ku3+ku2];
				q=0.50*(rsd[ku4+ku3+ku2+1]*rsd[ku4+ku3+ku2+1]
						+rsd[ku4+ku3+ku2+2]*rsd[ku4+ku3+ku2+2]
						+rsd[ku4+ku3+ku2+3]*rsd[ku4+ku3+ku2+3])
					/rsd[ku4+ku3+ku2];
				flux[kfl2+1]=rsd[ku4+ku3+ku2+1]*u21+C2*(rsd[ku4+ku3+ku2+4]-q);
				flux[kfl2+2]=rsd[ku4+ku3+ku2+2]*u21;
				flux[kfl2+3]=rsd[ku4+ku3+ku2+3]*u21;
				flux[kfl2+4]=(C1*rsd[ku4+ku3+ku2+4]-C2*q)*u21;
			}
			for(int i=ist; i<iend; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=tx2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int i=ist; i<nx; i++){
				kfl2 = i*fl_const;
				ku2 = i*ku2_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21i=tmp*rsd[ku4+ku3+ku2+1];
				u31i=tmp*rsd[ku4+ku3+ku2+2];
				u41i=tmp*rsd[ku4+ku3+ku2+3];
				u51i=tmp*rsd[ku4+ku3+ku2+4];
				ku2 -= ku2_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21im1=tmp*rsd[ku4+ku3+ku2+1];
				u31im1=tmp*rsd[ku4+ku3+ku2+2];
				u41im1=tmp*rsd[ku4+ku3+ku2+3];
				u51im1=tmp*rsd[ku4+ku3+ku2+4];
				flux[kfl2+1]=(4.0/3.0)*tx3*(u21i-u21im1);
				flux[kfl2+2]=tx3*(u31i-u31im1);
				flux[kfl2+3]=tx3*(u41i-u41im1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*tx3*((u21i*u21i+u31i*u31i+u41i*u41i)
							-(u21im1*u21im1+u31im1*u31im1+u41im1*u41im1))
					+(1.0/6.0)
					*tx3*(u21i*u21i-u21im1*u21im1)
					+C1*C5*tx3*(u51i-u51im1);
			}
			for(int i=ist; i<iend; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				int ku2P = ku2+ku2_const;
				int ku2M = ku2-ku2_const;
				frct[ku4+ku3+ku2]+=dx1*tx1*(rsd[ku4+ku3+ku2M]
							-2.0*rsd[ku4+ku3+ku2]
							+rsd[ku4+ku3+ku2P]);
				frct[ku4+ku3+ku2+1]+=tx3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dx2*tx1*(rsd[ku4+ku3+ku2M+1]
							-2.0*rsd[ku4+ku3+ku2+1]
							+rsd[ku4+ku3+ku2P+1]);
				frct[ku4+ku3+ku2+2]+=tx3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
					+dx3*tx1*(rsd[ku4+ku3+ku2M+2]
							-2.0*rsd[ku4+ku3+ku2+2]
							+rsd[ku4+ku3+ku2P+2]);
				frct[ku4+ku3+ku2+3]+=tx3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
					+dx4*tx1*(rsd[ku4+ku3+ku2M+3]
							-2.0*rsd[ku4+ku3+ku2+3]
							+rsd[ku4+ku3+ku2P+3]);
				frct[ku4+ku3+ku2+4]+=tx3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
					+dx5*tx1*(rsd[ku4+ku3+ku2M+4]
							-2.0*rsd[ku4+ku3+ku2+4]
							+rsd[ku4+ku3+ku2P+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			ku2 = 2*ku2_const;
			int ku2P = ku2+ku2_const;
			int ku2P2 = ku2P+ku2_const;
			int ku2M = ku2-ku2_const;
			int ku2M2 = ku2M-ku2_const;
			for(int m=0; m<5; m++){
				frct[ku4+ku3+ku2M+m]-=dssp*(+5.0*rsd[ku4+ku3+ku2M+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+rsd[ku4+ku3+ku2P+m]);
				frct[ku4+ku3+ku2+m]-=dssp*(-4.0*rsd[ku4+ku3+ku2M+m]
							+6.0*rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4+ku3+ku2P+m]
							+rsd[ku4+ku3+ku2P2+m]);
			}
			for(int i=3; i<nx-3; i++){
				ku2 = i*ku2_const;
				ku2P = ku2+ku2_const;
				ku2P2 = ku2P+ku2_const;
				ku2M = ku2-ku2_const;
				ku2M2 = ku2M-ku2_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=dssp*(rsd[ku4+ku3+ku2M2+m]
								-4.0*rsd[ku4+ku3+ku2M+m]
								+6.0*rsd[ku4+ku3+ku2+m]
								-4.0*rsd[ku4+ku3+ku2P+m]
								+rsd[ku4+ku3+ku2P2+m] );
				}
			}
			ku2 = (nx-4)*ku2_const;
			ku2P = ku2+ku2_const;
			ku2P2 = ku2P+ku2_const;
			ku2M = ku2-ku2_const;
			ku2M2 = ku2M-ku2_const;
			for(int m=0; m<5; m++){
				frct[ku4+ku3+ku2P+m]-=dssp*(rsd[ku4+ku3+ku2M+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+6.0*rsd[ku4+ku3+ku2P+m]
							-4.0*rsd[ku4+ku3+ku2P2+m] );
				frct[ku4+ku3+ku2P2+m]-=dssp*(rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4+ku3+ku2P+m]
							+5.0*rsd[ku4+ku3+ku2P2+m] );
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * eta-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, nz-2, [&] (int k) {
		double q;
		double u21, u31, u41;
		double tmp;
		double u21j, u31j, u41j, u51j;
		double u21jm1, u31jm1, u41jm1, u51jm1;
		int ku4, ku3, ku2;
		int kfl2, kfl2P, kfl2M;
		int ku3P, ku3P2, ku3M, ku3M2;
		std::vector<double> flux(ISIZ1*5, 0);

		ku4 = k*ku4_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int j=0; j<ny; j++){
				kfl2 = j*fl_const;
				ku3 = j*ku3_const;
				flux[kfl2]=rsd[ku4+ku3+ku2+2];
				u31=rsd[ku4+ku3+ku2+2]/rsd[ku4+ku3+ku2];
				q=0.50*(rsd[ku4+ku3+ku2+1]*rsd[ku4+ku3+ku2+1]
						+rsd[ku4+ku3+ku2+2]*rsd[ku4+ku3+ku2+2]
						+rsd[ku4+ku3+ku2+3]*rsd[ku4+ku3+ku2+3])
					/rsd[ku4+ku3+ku2];
				flux[kfl2+1]=rsd[ku4+ku3+ku2+1]*u31;
				flux[kfl2+2]=rsd[ku4+ku3+ku2+2]*u31+C2*(rsd[ku4+ku3+ku2+4]-q);
				flux[kfl2+3]=rsd[ku4+ku3+ku2+3]*u31;
				flux[kfl2+4]=(C1*rsd[ku4+ku3+ku2+4]-C2*q)*u31;
			}
			for(int j=jst; j<jend;j++){
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=ty2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int j=jst; j<ny; j++){
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21j=tmp*rsd[ku4+ku3+ku2+1];
				u31j=tmp*rsd[ku4+ku3+ku2+2];
				u41j=tmp*rsd[ku4+ku3+ku2+3];
				u51j=tmp*rsd[ku4+ku3+ku2+4];
				ku3 -= ku3_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21jm1=tmp*rsd[ku4+ku3+ku2+1];
				u31jm1=tmp*rsd[ku4+ku3+ku2+2];
				u41jm1=tmp*rsd[ku4+ku3+ku2+3];
				u51jm1=tmp*rsd[ku4+ku3+ku2+4];
				flux[kfl2+1]=ty3*(u21j-u21jm1);
				flux[kfl2+2]=(4.0/3.0)*ty3*(u31j-u31jm1);
				flux[kfl2+3]=ty3*(u41j-u41jm1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*ty3*((u21j*u21j+u31j*u31j+u41j*u41j)
							-(u21jm1*u21jm1+u31jm1*u31jm1+u41jm1*u41jm1))
					+(1.0/6.0)
					*ty3*(u31j*u31j-u31jm1*u31jm1)
					+C1*C5*ty3*(u51j-u51jm1);
			}
			for(int j=jst; j<jend; j++){
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				ku3P = ku3+ku3_const, ku3M = ku3-ku3_const;
				frct[ku4+ku3+ku2]+=dy1*ty1*(rsd[ku4+ku3M+ku2]
							-2.0*rsd[ku4+ku3+ku2]
							+rsd[ku4+ku3P+ku2]);
				frct[ku4+ku3+ku2+1]+=ty3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dy2*ty1*(rsd[ku4+ku3M+ku2+1]
							-2.0*rsd[ku4+ku3+ku2+1]
							+rsd[ku4+ku3P+ku2+1]);
				frct[ku4+ku3+ku2+2]+=ty3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
					+dy3*ty1*(rsd[ku4+ku3M+ku2+2]
							-2.0*rsd[ku4+ku3+ku2+2]
							+rsd[ku4+ku3P+ku2+2]);
				frct[ku4+ku3+ku2+3]+=ty3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
					+dy4*ty1*(rsd[ku4+ku3M+ku2+3]
							-2.0*rsd[ku4+ku3+ku2+3]
							+rsd[ku4+ku3P+ku2+3]);
				frct[ku4+ku3+ku2+4]+=ty3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
					+dy5*ty1*(rsd[ku4+ku3M+ku2+4]
							-2.0*rsd[ku4+ku3+ku2+4]
							+rsd[ku4+ku3P+ku2+4] );
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			ku3 = 2*ku3_const;
			ku3P = ku3+ku3_const, ku3P2 = ku3P+ku3_const, ku3M = ku3-ku3_const, ku3M2 = ku3M-ku3_const;
			for(int m=0; m<5; m++){
				frct[ku4+ku3M+ku2+m]-=dssp*(+5.0*rsd[ku4+ku3M+ku2+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+rsd[ku4+ku3P+ku2+m]);
				frct[ku4+ku3+ku2+m]-=dssp*(-4.0*rsd[ku4+ku3M+ku2+m]
							+6.0*rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4+ku3P+ku2+m]
							+rsd[ku4+ku3P2+ku2+m]);
			}
			for(int j=3; j<ny-3; j++){
				ku3 = j*ku3_const;
				ku3P = ku3+ku3_const, ku3P2 = ku3P+ku3_const, ku3M = ku3-ku3_const, ku3M2 = ku3M-ku3_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=dssp*(rsd[ku4+ku3M2+ku2+m]
								-4.0*rsd[ku4+ku3M+ku2+m]
								+6.0*rsd[ku4+ku3+ku2+m]
								-4.0*rsd[ku4+ku3P+ku2+m]
								+rsd[ku4+ku3P2+ku2+m]);
				}
			}
			ku3 = (ny-4)*ku3_const;
			ku3P = ku3+ku3_const, ku3P2 = ku3P+ku3_const, ku3M = ku3-ku3_const, ku3M2 = ku3M-ku3_const;
			for(int m=0; m<5; m++){
				frct[ku4+ku3P+ku2+m]-=dssp*(rsd[ku4+ku3M+ku2+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+6.0*rsd[ku4+ku3P+ku2+m]
							-4.0*rsd[ku4+ku3P2+ku2+m]);
				frct[ku4+ku3P2+ku2+m]-=dssp*(rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4+ku3P+ku2+m]
							+5.0*rsd[ku4+ku3P2+ku2+m]);
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * zeta-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+jst, jend-jst, [&] (int j) {
		double q;
		double u21, u31, u41;
		double tmp;
		double u21k, u31k, u41k, u51k;
		double u21km1, u31km1, u41km1, u51km1;
		int ku4, ku3, ku2;
		int kfl2, kfl2P, kfl2M;
		int ku4P, ku4P2, ku4M, ku4M2;
		std::vector<double> flux(ISIZ1*5, 0);
		
		ku3 = j*ku3_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int k=0; k<nz; k++){
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku4 = k*ku4_const;
				flux[kfl2]=rsd[ku4+ku3+ku2+3];
				u41=rsd[ku4+ku3+ku2+3]/rsd[ku4+ku3+ku2];
				q=0.50*(rsd[ku4+ku3+ku2+1]*rsd[ku4+ku3+ku2+1]
						+rsd[ku4+ku3+ku2+2]*rsd[ku4+ku3+ku2+2]
						+rsd[ku4+ku3+ku2+3]*rsd[ku4+ku3+ku2+3])
					/rsd[ku4+ku3+ku2];
				flux[kfl2+1]=rsd[ku4+ku3+ku2+1]*u41;
				flux[kfl2+2]=rsd[ku4+ku3+ku2+2]*u41; 
				flux[kfl2+3]=rsd[ku4+ku3+ku2+3]*u41+C2*(rsd[ku4+ku3+ku2+4]-q);
				flux[kfl2+4]=(C1*rsd[ku4+ku3+ku2+4]-C2*q)*u41;
			}
			for(int k=1; k<nz-1; k++){
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku4 = k*ku4_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=tz2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int k=1; k<nz; k++){
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku4 = k*ku4_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21k=tmp*rsd[ku4+ku3+ku2+1];
				u31k=tmp*rsd[ku4+ku3+ku2+2];
				u41k=tmp*rsd[ku4+ku3+ku2+3];
				u51k=tmp*rsd[ku4+ku3+ku2+4];
				ku4 -= ku4_const;
				tmp=1.0/rsd[ku4+ku3+ku2];
				u21km1=tmp*rsd[ku4+ku3+ku2+1];
				u31km1=tmp*rsd[ku4+ku3+ku2+2];
				u41km1=tmp*rsd[ku4+ku3+ku2+3];
				u51km1=tmp*rsd[ku4+ku3+ku2+4];
				flux[kfl2+1]=tz3*(u21k-u21km1);
				flux[kfl2+2]=tz3*(u31k-u31km1);
				flux[kfl2+3]=(4.0/3.0)*tz3*(u41k-u41km1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*tz3*((u21k*u21k+u31k*u31k+u41k*u41k)
							-(u21km1*u21km1+u31km1*u31km1+u41km1*u41km1))
					+(1.0/6.0)
					*tz3*(u41k*u41k-u41km1*u41km1)
					+C1*C5*tz3*(u51k-u51km1);
			}
			for(int k=1; k<nz-1; k++){
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku4 = k*ku4_const;
				ku4P = ku4+ku4_const, ku4M = ku4-ku4_const;
				frct[ku4+ku3+ku2]+=dz1*tz1*(rsd[ku4P+ku3+ku2]
							-2.0*rsd[ku4+ku3+ku2]
							+rsd[ku4M+ku3+ku2]);
				frct[ku4+ku3+ku2+1]+=tz3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dz2*tz1*(rsd[ku4P+ku3+ku2+1]
							-2.0*rsd[ku4+ku3+ku2+1]
							+rsd[ku4M+ku3+ku2+1]);
				frct[ku4+ku3+ku2+2]+=tz3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
					+dz3*tz1*(rsd[ku4P+ku3+ku2+2]
							-2.0*rsd[ku4+ku3+ku2+2]
							+rsd[ku4M+ku3+ku2+2]);
				frct[ku4+ku3+ku2+3]+=tz3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
					+dz4*tz1*(rsd[ku4P+ku3+ku2+3]
							-2.0*rsd[ku4+ku3+ku2+3]
							+rsd[ku4M+ku3+ku2+3]);
				frct[ku4+ku3+ku2+4]+=tz3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
					+dz5*tz1*(rsd[ku4P+ku3+ku2+4]
							-2.0*rsd[ku4+ku3+ku2+4]
							+rsd[ku4M+ku3+ku2+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			ku4 = 2*ku4_const;
			ku4P = ku4+ku4_const, ku4P2 = ku4P+ku4_const, ku4M = ku4-ku4_const, ku4M2 = ku4M-ku4_const;
			for(int m=0; m<5; m++){
				frct[ku4M+ku3+ku2+m]-=dssp*(+5.0*rsd[ku4M+ku3+ku2+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+rsd[ku4P+ku3+ku2+m]);
				frct[ku4+ku3+ku2+m]-=dssp*(-4.0*rsd[ku4M+ku3+ku2+m]
							+6.0*rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4P+ku3+ku2+m]
							+rsd[ku4P2+ku3+ku2+m]);
			}
			for(int k=3; k<nz-3; k++){
				ku4 = k*ku4_const;
				ku4P = ku4+ku4_const, ku4P2 = ku4P+ku4_const, ku4M = ku4-ku4_const, ku4M2 = ku4M-ku4_const;
				for(int m=0; m<5; m++){
					frct[ku4+ku3+ku2+m]-=dssp*(rsd[ku4M2+ku3+ku2+m]
								-4.0*rsd[ku4M+ku3+ku2+m]
								+6.0*rsd[ku4+ku3+ku2+m]
								-4.0*rsd[ku4P+ku3+ku2+m]
								+rsd[ku4P2+ku3+ku2+m]);
				}
			}
			ku4 = (nz-4)*ku4_const;
			ku4P = ku4+ku4_const, ku4P2 = ku4P+ku4_const, ku4M = ku4-ku4_const, ku4M2 = ku4M-ku4_const;
			for(int m=0; m<5; m++){
				frct[ku4P+ku3+ku2+m]-=dssp*(rsd[ku4M+ku3+ku2+m]
							-4.0*rsd[ku4+ku3+ku2+m]
							+6.0*rsd[ku4P+ku3+ku2+m]
							-4.0*rsd[ku4P2+ku3+ku2+m]);
				frct[ku4P2+ku3+ku2+m]-=dssp*(rsd[ku4+ku3+ku2+m]
							-4.0*rsd[ku4P+ku3+ku2+m]
							+5.0*rsd[ku4P2+ku3+ku2+m]);
			}
		}
	});
}

/*
 * ---------------------------------------------------------------------
 * compute the solution error
 * ---------------------------------------------------------------------
 */
void error(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	double tmp;
	std::vector<double> u000ijk(5);
	std::fill(errnm.begin(), errnm.end(), 0.0);
	for(int k=1; k<nz-1; k++){
		ku4 = k*ku4_const;
		for(int j=jst; j<jend; j++){
			ku3 = j*ku3_const;
			for(int i=ist; i<iend; i++){
				ku2 = i*ku2_const;
				exact(i, j, k, u000ijk);
				for(int m=0; m<5; m++){
					tmp=(u000ijk[m]-u[ku4+ku3+ku2+m]);
					errnm[m]=errnm[m]+tmp*tmp;
				}
			}
		}
	}
	std::transform(errnm.begin(), errnm.end(), errnm.begin(), [&](double x){
		return sqrt(x/((nx0-2)*(ny0-2)*(nz0-2)));
	});
}

/*
 * ---------------------------------------------------------------------
 * compute the exact solution at (i,j,k)
 * ---------------------------------------------------------------------
 */
void exact(int i, int j, int k, std::vector<double> &u000ijk){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	double xi, eta, zeta;
	xi=((double)i)/(nx0-1);
	eta=((double)j)/(ny0-1);
	zeta=((double)k)/(nz-1);
	for(int m=0; m<5; m++){
		u000ijk[m]=ce[(0*5)+m]+
			(ce[(1*5)+m]+
			 (ce[(4*5)+m]+
			  (ce[(7*5)+m]+
			   ce[(10*5)+m]*xi)*xi)*xi)*xi+
			(ce[(2*5)+m]+
			 (ce[(5*5)+m]+
			  (ce[(8*5)+m]+
			   ce[(11*5)+m]*eta)*eta)*eta)*eta+
			(ce[(3*5)+m]+
			 (ce[(6*5)+m]+
			  (ce[(9*5)+m]+
			   ce[(12*5)+m]*zeta)*zeta)*zeta)*zeta;
	}
}

/*
 * ---------------------------------------------------------------------
 * compute the lower triangular part of the jacobian matrix
 * ---------------------------------------------------------------------
 */
void jacld(int k, int worker_id){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	int ku2, ku3, ku4;
	int kq2, kq3;
	int ks2, ks3, ks4;
	double r43;
	double c1345;
	double c34;
	double tmp1, tmp2, tmp3;
	r43=(4.0/3.0);
	c1345=C1*C3*C4*C5;
	c34=C3*C4;
	int chunk = (jend-jst)/num_workers;
	int residual = (jend-jst)%num_workers;

	int jst_worker = worker_id*chunk + ((residual>1 && (worker_id<residual)) ? worker_id : residual) + jst;
	int jend_worker = worker_id*chunk + chunk + ((residual>1 && (worker_id<residual)) ? worker_id+1 : residual) + jst;
	

	if(worker_id==0) jst_worker = jst;
	if(worker_id==num_workers-1) jend_worker = jend;

	ku4 = k*ku4_const;
	int ku4M = (k-1)*ku4_const;
	kq3 = k*kq3_const;
	int kq3M = (k-1)*kq3_const;
	for(int j=jst_worker; j<jend_worker; j++){
		kq2 = j*kq2_const;
		int kq2M = (j-1)*kq2_const;
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		int ku3M = (j-1)*ku3_const;
		for(int i=ist; i<iend; i++){
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			int ku2M = (i-1)*ku2_const;
			/*
			 * ---------------------------------------------------------------------
			 * form the block daigonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			d[ks4+ks3]=1.0+dt*2.0*(tx1*dx1+ty1*dy1+tz1*dz1);
			d[ks4+ks3+ks2T1]=0.0;
			d[ks4+ks3+ks2T2]=0.0;
			d[ks4+ks3+ks2T3]=0.0;
			d[ks4+ks3+ks2T4]=0.0;
			d[ks4+ks3+1]=-dt*2.0
				*(tx1*r43+ty1+tz1)*c34*tmp2*u[ku4+ku3+ku2+1];
			d[ks4+ks3+ks2T1+1]=1.0
				+dt*2.0*c34*tmp1*(tx1*r43+ty1+tz1)
				+dt*2.0*(tx1*dx2+ty1*dy2+tz1*dz2);
			d[ks4+ks3+ks2T2+1]=0.0;
			d[ks4+ks3+ks2T3+1]=0.0;
			d[ks4+ks3+ks2T4+1]=0.0;
			d[ks4+ks3+2]=-dt*2.0 
				*(tx1+ty1*r43+tz1)*c34*tmp2*u[ku4+ku3+ku2+2];
			d[ks4+ks3+ks2T1+2]=0.0;
			d[ks4+ks3+ks2T2+2]=1.0
				+dt*2.0*c34*tmp1*(tx1+ty1*r43+tz1)
				+dt*2.0*(tx1*dx3+ty1*dy3+tz1*dz3);
			d[ks4+ks3+ks2T3+2]=0.0;
			d[ks4+ks3+ks2T4+2]=0.0;
			d[ks4+ks3+3]=-dt*2.0
				*(tx1+ty1+tz1*r43)*c34*tmp2*u[ku4+ku3+ku2+3];
			d[ks4+ks3+ks2T1+3]=0.0;
			d[ks4+ks3+ks2T2+3]=0.0;
			d[ks4+ks3+ks2T3+3]=1.0
				+dt*2.0*c34*tmp1*(tx1+ty1+tz1*r43)
				+dt*2.0*(tx1*dx4+ty1*dy4+tz1*dz4);
			d[ks4+ks3+ks2T4+3]=0.0;
			d[ks4+ks3+4]=-dt*2.0
				*(((tx1*(r43*c34-c1345)
								+ty1*(c34-c1345)
								+tz1*(c34-c1345))*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1])
							+(tx1*(c34-c1345)
								+ty1*(r43*c34-c1345)
								+tz1*(c34-c1345))*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2])
							+(tx1*(c34-c1345)
								+ty1*(c34-c1345)
								+tz1*(r43*c34-c1345))*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
				  )*tmp3
						+(tx1+ty1+tz1)*c1345*tmp2*u[ku4+ku3+ku2+4]);
			d[ks4+ks3+ks2T1+4]=dt*2.0*tmp2*u[ku4+ku3+ku2+1]
				*(tx1*(r43*c34-c1345)
						+ty1*(c34-c1345)
						+tz1*(c34-c1345));
			d[ks4+ks3+ks2T2+4]=dt*2.0*tmp2*u[ku4+ku3+ku2+2]
				*(tx1*(c34-c1345)
						+ty1*(r43*c34-c1345)
						+tz1*(c34-c1345));
			d[ks4+ks3+ks2T3+4]=dt*2.0*tmp2*u[ku4+ku3+ku2+3]
				*(tx1*(c34-c1345)
						+ty1*(c34-c1345)
						+tz1*(r43*c34-c1345));
			d[ks4+ks3+ks2T4+4]=1.0
				+dt*2.0*(tx1+ty1+tz1)*c1345*tmp1
				+dt*2.0*(tx1*dx5+ty1*dy5+tz1*dz5);
			/*
			 * ---------------------------------------------------------------------
			 * form the first block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3M+kq2+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			a[ks4+ks3]=-dt*tz1*dz1;
			a[ks4+ks3+ks2T1]=0.0;
			a[ks4+ks3+ks2T2]=0.0;
			a[ks4+ks3+ks2T3]=-dt*tz2;
			a[ks4+ks3+ks2T4]=0.0;
			a[ks4+ks3+1]=-dt*tz2
				*(-(u[ku4M+ku3+ku2+1]*u[ku4M+ku3+ku2+3])*tmp2)
				-dt*tz1*(-c34*tmp2*u[ku4M+ku3+ku2+1]);
			a[ks4+ks3+ks2T1+1]=-dt*tz2*(u[ku4M+ku3+ku2+3]*tmp1)
				-dt*tz1*c34*tmp1
				-dt*tz1*dz2;
			a[ks4+ks3+ks2T2+1]=0.0;
			a[ks4+ks3+ks2T3+1]=-dt*tz2*(u[ku4M+ku3+ku2+1]*tmp1);
			a[ks4+ks3+ks2T4+1]=0.0;
			a[ks4+ks3+2]=-dt*tz2
				*(-(u[ku4M+ku3+ku2+2]*u[ku4M+ku3+ku2+3])*tmp2)
				-dt*tz1*(-c34*tmp2*u[ku4M+ku3+ku2+2]);
			a[ks4+ks3+ks2T1+2]=0.0;
			a[ks4+ks3+ks2T2+2]=-dt*tz2*(u[ku4M+ku3+ku2+3]*tmp1)
				-dt*tz1*(c34*tmp1)
				-dt*tz1*dz3;
			a[ks4+ks3+ks2T3+2]=-dt*tz2*(u[ku4M+ku3+ku2+2]*tmp1);
			a[ks4+ks3+ks2T4+2]=0.0;
			a[ks4+ks3+3]=-dt*tz2
				*(-(u[ku4M+ku3+ku2+3]*tmp1)*(u[ku4M+ku3+ku2+3]*tmp1)
						+C2*qs[kq3M+kq2+i]*tmp1)
				-dt*tz1*(-r43*c34*tmp2*u[ku4M+ku3+ku2+3]);
			a[ks4+ks3+ks2T1+3]=-dt*tz2
				*(-C2*(u[ku4M+ku3+ku2+1]*tmp1));
			a[ks4+ks3+ks2T2+3]=-dt*tz2
				*(-C2*(u[ku4M+ku3+ku2+2]*tmp1));
			a[ks4+ks3+ks2T3+3]=-dt*tz2*(2.0-C2)
				*(u[ku4M+ku3+ku2+3]*tmp1)
				-dt*tz1*(r43*c34*tmp1)
				-dt*tz1*dz4;
			a[ks4+ks3+ks2T4+3]=-dt*tz2*C2;
			a[ks4+ks3+4]=-dt*tz2
				*((C2*2.0*qs[kq3M+kq2+i]-C1*u[ku4M+ku3+ku2+4])
						*u[ku4M+ku3+ku2+3]*tmp2)
				-dt*tz1
				*(-(c34-c1345)*tmp3*(u[ku4M+ku3+ku2+1]*u[ku4M+ku3+ku2+1])
						-(c34-c1345)*tmp3*(u[ku4M+ku3+ku2+2]*u[ku4M+ku3+ku2+2])
						-(r43*c34-c1345)*tmp3*(u[ku4M+ku3+ku2+3]*u[ku4M+ku3+ku2+3])
						-c1345*tmp2*u[ku4M+ku3+ku2+4]);
			a[ks4+ks3+ks2T1+4]=-dt*tz2
				*(-C2*(u[ku4M+ku3+ku2+1]*u[ku4M+ku3+ku2+3])*tmp2)
				-dt*tz1*(c34-c1345)*tmp2*u[ku4M+ku3+ku2+1];
			a[ks4+ks3+ks2T2+4]=-dt*tz2
				*(-C2*(u[ku4M+ku3+ku2+2]*u[ku4M+ku3+ku2+3])*tmp2)
				-dt*tz1*(c34-c1345)*tmp2*u[ku4M+ku3+ku2+2];
			a[ks4+ks3+ks2T3+4]=-dt*tz2
				*(C1*(u[ku4M+ku3+ku2+4]*tmp1)
						-C2*(qs[kq3M+kq2+i]*tmp1
							+u[ku4M+ku3+ku2+3]*u[ku4M+ku3+ku2+3]*tmp2))
				-dt*tz1*(r43*c34-c1345)*tmp2*u[ku4M+ku3+ku2+3];
			a[ks4+ks3+ks2T4+4]=-dt*tz2
				*(C1*(u[ku4M+ku3+ku2+3]*tmp1))
				-dt*tz1*c1345*tmp1
				-dt*tz1*dz5;
			/*
			 * ---------------------------------------------------------------------
			 * form the second block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2M+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			b[ks4+ks3]=-dt*ty1*dy1;
			b[ks4+ks3+ks2T1]=0.0;
			b[ks4+ks3+ks2T2]=-dt*ty2;
			b[ks4+ks3+ks2T3]=0.0;
			b[ks4+ks3+ks2T4]=0.0;
			b[ks4+ks3+1]=-dt*ty2
				*(-(u[ku4+ku3M+ku2+1]*u[ku4+ku3M+ku2+2])*tmp2)
				-dt*ty1*(-c34*tmp2*u[ku4+ku3M+ku2+1]);
			b[ks4+ks3+ks2T1+1]=-dt*ty2*(u[ku4+ku3M+ku2+2]*tmp1)
				-dt*ty1*(c34*tmp1)
				-dt*ty1*dy2;
			b[ks4+ks3+ks2T2+1]=-dt*ty2*(u[ku4+ku3M+ku2+1]*tmp1);
			b[ks4+ks3+ks2T3+1]=0.0;
			b[ks4+ks3+ks2T4+1]=0.0;
			b[ks4+ks3+2]=-dt*ty2
				*(-(u[ku4+ku3M+ku2+2]*tmp1)*(u[ku4+ku3M+ku2+2]*tmp1)
						+C2*(qs[kq3+kq2M+i]*tmp1))
				-dt*ty1*(-r43*c34*tmp2*u[ku4+ku3M+ku2+2]);
			b[ks4+ks3+ks2T1+2]=-dt*ty2
				*(-C2*(u[ku4+ku3M+ku2+1]*tmp1));
			b[ks4+ks3+ks2T2+2]=-dt*ty2*((2.0-C2)*(u[ku4+ku3M+ku2+2]*tmp1))
				-dt*ty1*(r43*c34*tmp1)
				-dt*ty1*dy3;
			b[ks4+ks3+ks2T3+2]=-dt*ty2*(-C2*(u[ku4+ku3M+ku2+3]*tmp1));
			b[ks4+ks3+ks2T4+2]=-dt*ty2*C2;
			b[ks4+ks3+3]=-dt*ty2
				*(-(u[ku4+ku3M+ku2+2]*u[ku4+ku3M+ku2+3])*tmp2)
				-dt*ty1*(-c34*tmp2*u[ku4+ku3M+ku2+3]);
			b[ks4+ks3+ks2T1+3]=0.0;
			b[ks4+ks3+ks2T2+3]=-dt*ty2*(u[ku4+ku3M+ku2+3]*tmp1);
			b[ks4+ks3+ks2T3+3]=-dt*ty2*(u[ku4+ku3M+ku2+2]*tmp1)
				-dt*ty1*(c34*tmp1)
				-dt*ty1*dy4;
			b[ks4+ks3+ks2T4+3]=0.0;
			b[ks4+ks3+4]=-dt*ty2
				*((C2*2.0*qs[kq3+kq2M+i]-C1*u[ku4+ku3M+ku2+4])
						*(u[ku4+ku3M+ku2+2]*tmp2))
				-dt*ty1
				*(-(c34-c1345)*tmp3*(u[ku4+ku3M+ku2+1]*u[ku4+ku3M+ku2+1])
						-(r43*c34-c1345)*tmp3*(u[ku4+ku3M+ku2+2]*u[ku4+ku3M+ku2+2])
						-(c34-c1345)*tmp3*(u[ku4+ku3M+ku2+3]*u[ku4+ku3M+ku2+3])
						-c1345*tmp2*u[ku4+ku3M+ku2+4]);
			b[ks4+ks3+ks2T1+4]=-dt*ty2
				*(-C2*(u[ku4+ku3M+ku2+1]*u[ku4+ku3M+ku2+2])*tmp2)
				-dt*ty1*(c34-c1345)*tmp2*u[ku4+ku3M+ku2+1];
			b[ks4+ks3+ks2T2+4]=-dt*ty2
				*(C1*(u[ku4+ku3M+ku2+4]*tmp1)
						-C2*(qs[kq3+kq2M+i]*tmp1
							+u[ku4+ku3M+ku2+2]*u[ku4+ku3M+ku2+2]*tmp2))
				-dt*ty1*(r43*c34-c1345)*tmp2*u[ku4+ku3M+ku2+2];
			b[ks4+ks3+ks2T3+4]=-dt*ty2
				*(-C2*(u[ku4+ku3M+ku2+2]*u[ku4+ku3M+ku2+3])*tmp2)
				-dt*ty1*(c34-c1345)*tmp2*u[ku4+ku3M+ku2+3];
			b[ks4+ks3+ks2T4+4]=-dt*ty2
				*(C1*(u[ku4+ku3M+ku2+2]*tmp1))
				-dt*ty1*c1345*tmp1
				-dt*ty1*dy5;
			/*
			 * ---------------------------------------------------------------------
			 * form the third block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2+i-1];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			c[ks4+ks3]=-dt*tx1*dx1;
			c[ks4+ks3+ks2T1]=-dt*tx2;
			c[ks4+ks3+ks2T2]=0.0;
			c[ks4+ks3+ks2T3]=0.0;
			c[ks4+ks3+ks2T4]=0.0;
			c[ks4+ks3+1]=-dt*tx2
				*(-(u[ku4+ku3+ku2M+1]*tmp1)*(u[ku4+ku3+ku2M+1]*tmp1)
						+C2*qs[kq3+kq2+i-1]*tmp1)
				-dt*tx1*(-r43*c34*tmp2*u[ku4+ku3+ku2M+1]);
			c[ks4+ks3+ks2T1+1]=-dt*tx2
				*((2.0-C2)*(u[ku4+ku3+ku2M+1]*tmp1))
				-dt*tx1*(r43*c34*tmp1)
				-dt*tx1*dx2;
			c[ks4+ks3+ks2T2+1]=-dt*tx2
				*(-C2*(u[ku4+ku3+ku2M+2]*tmp1));
			c[ks4+ks3+ks2T3+1]=-dt*tx2
				*(-C2*(u[ku4+ku3+ku2M+3]*tmp1));
			c[ks4+ks3+ks2T4+1]=-dt*tx2*C2;
			c[ks4+ks3+2]=-dt*tx2
				*(-(u[ku4+ku3+ku2M+1]*u[ku4+ku3+ku2M+2])*tmp2)
				-dt*tx1*(-c34*tmp2*u[ku4+ku3+ku2M+2]);
			c[ks4+ks3+ks2T1+2]=-dt*tx2*(u[ku4+ku3+ku2M+2]*tmp1);
			c[ks4+ks3+ks2T2+2]=-dt*tx2*(u[ku4+ku3+ku2M+1]*tmp1)
				-dt*tx1*(c34*tmp1)
				-dt*tx1*dx3;
			c[ks4+ks3+ks2T3+2]=0.0;
			c[ks4+ks3+ks2T4+2]=0.0;
			c[ks4+ks3+3]=-dt*tx2
				*(-(u[ku4+ku3+ku2M+1]*u[ku4+ku3+ku2M+3])*tmp2)
				-dt*tx1*(-c34*tmp2*u[ku4+ku3+ku2M+3]);
			c[ks4+ks3+ks2T1+3]=-dt*tx2*(u[ku4+ku3+ku2M+3]*tmp1);
			c[ks4+ks3+ks2T2+3]=0.0;
			c[ks4+ks3+ks2T3+3]=-dt*tx2*(u[ku4+ku3+ku2M+1]*tmp1)
				-dt*tx1*(c34*tmp1)-dt*tx1*dx4;
			c[ks4+ks3+ks2T4+3]=0.0;
			c[ks4+ks3+4]=-dt*tx2
				*((C2*2.0*qs[kq3+kq2+i-1]-C1*u[ku4+ku3+ku2M+4])
						*u[ku4+ku3+ku2M+1]*tmp2)
				-dt*tx1
				*(-(r43*c34-c1345)*tmp3*(u[ku4+ku3+ku2M+1]*u[ku4+ku3+ku2M+1])
						-(c34-c1345)*tmp3*(u[ku4+ku3+ku2M+2]*u[ku4+ku3+ku2M+2])
						-(c34-c1345)*tmp3*(u[ku4+ku3+ku2M+3]*u[ku4+ku3+ku2M+3])
						-c1345*tmp2*u[ku4+ku3+ku2M+4]);
			c[ks4+ks3+ks2T1+4]=-dt*tx2
				*(C1*(u[ku4+ku3+ku2M+4]*tmp1)
						-C2*(u[ku4+ku3+ku2M+1]*u[ku4+ku3+ku2M+1]*tmp2
							+qs[kq3+kq2+i-1]*tmp1))
				-dt*tx1*(r43*c34-c1345)*tmp2*u[ku4+ku3+ku2M+1];
			c[ks4+ks3+ks2T2+4]=-dt*tx2
				*(-C2*(u[ku4+ku3+ku2M+2]*u[ku4+ku3+ku2M+1])*tmp2)
				-dt*tx1*(c34-c1345)*tmp2*u[ku4+ku3+ku2M+2];
			c[ks4+ks3+ks2T3+4]=-dt*tx2
				*(-C2*(u[ku4+ku3+ku2M+3]*u[ku4+ku3+ku2M+1])*tmp2)
				-dt*tx1*(c34-c1345)*tmp2*u[ku4+ku3+ku2M+3];
			c[ks4+ks3+ks2T4+4]=-dt*tx2
				*(C1*(u[ku4+ku3+ku2M+1]*tmp1))
				-dt*tx1*c1345*tmp1
				-dt*tx1*dx5;
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * compute the upper triangular part of the jacobian matrix
 * ---------------------------------------------------------------------
 */
void jacu(int k, int worker_id){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	double r43;
	double c1345;
	double c34;
	double tmp1, tmp2, tmp3;
	int ku2, ku3, ku4;
	int kq2, kq3;
	int ks2, ks3, ks4;
	r43=(4.0/3.0);
	c1345=C1*C3*C4*C5;
	c34=C3*C4;
	worker_id = num_workers-1-worker_id;
	int chunk = (jend-jst)/num_workers;
	int residual = (jend-jst)%num_workers;

	int jst_worker = worker_id*chunk + ((residual>1 && (worker_id<residual)) ? worker_id : residual) + jst-1;
	int jend_worker = worker_id*chunk + chunk + ((residual>1 && (worker_id<residual)) ? worker_id+1 : residual) + jst-1;

	if(worker_id==0) jst_worker = jst-1;
	if(worker_id==num_workers-1) jend_worker = jend-1;

	ku4 = k*ku4_const;
	int ku4P = ku4+ku4_const;
	kq3 = k*kq3_const;
	int kq3P = kq3+kq3_const;
	for(int j=jend_worker; j>jst_worker; j--){
		kq2 = j*kq2_const;
		int kq2P = (j+1)*kq2_const;
		ks4 = j*ks4_const;
		ku3 = j*ku3_const;
		int ku3P = (j+1)*ku3_const;
		for (int i=iend-1; i>=ist; i--) {
			ks3 = i*ks3_const;
			ku2 = i*ku2_const;
			int ku2P = (i+1)*ku2_const;
			/*
			 * ---------------------------------------------------------------------
			 * form the block daigonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			d[ks4+ks3]=1.0+dt*2.0*(tx1*dx1+ty1*dy1+tz1*dz1);
			d[ks4+ks3+ks2T1]=0.0;
			d[ks4+ks3+ks2T2]=0.0;
			d[ks4+ks3+ks2T3]=0.0;
			d[ks4+ks3+ks2T4]=0.0;
			d[ks4+ks3+1]=dt*2.0
				*(-tx1*r43-ty1-tz1)
				*(c34*tmp2*u[ku4+ku3+ku2+1]);
			d[ks4+ks3+ks2T1+1]=1.0
				+dt*2.0*c34*tmp1 
				*(tx1*r43+ty1+tz1)
				+dt*2.0*(tx1*dx2+ty1*dy2+tz1*dz2);
			d[ks4+ks3+ks2T2+1]=0.0;
			d[ks4+ks3+ks2T3+1]=0.0;
			d[ks4+ks3+ks2T4+1]=0.0;
			d[ks4+ks3+2]=dt*2.0
				*(-tx1-ty1*r43-tz1)
				*(c34*tmp2*u[ku4+ku3+ku2+2]);
			d[ks4+ks3+ks2T1+2]=0.0;
			d[ks4+ks3+ks2T2+2]=1.0
				+dt*2.0*c34*tmp1
				*(tx1+ty1*r43+tz1)
				+dt*2.0*(tx1*dx3+ty1*dy3+tz1*dz3);
			d[ks4+ks3+ks2T3+2]=0.0;
			d[ks4+ks3+ks2T4+2]=0.0;
			d[ks4+ks3+3]=dt*2.0
				*(-tx1-ty1-tz1*r43)
				*(c34*tmp2*u[ku4+ku3+ku2+3]);
			d[ks4+ks3+ks2T1+3]=0.0;
			d[ks4+ks3+ks2T2+3]=0.0;
			d[ks4+ks3+ks2T3+3]=1.0
				+dt*2.0*c34*tmp1
				*(tx1+ty1+tz1*r43)
				+dt*2.0*(tx1*dx4+ty1*dy4+tz1*dz4);
			d[ks4+ks3+ks2T4+3]=0.0;
			d[ks4+ks3+4]=-dt*2.0
				*(((tx1*(r43*c34-c1345)
								+ty1*(c34-c1345)
								+tz1*(c34-c1345))*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1])
							+(tx1*(c34-c1345)
								+ty1*(r43*c34-c1345)
								+tz1*(c34-c1345))*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2])
							+(tx1*(c34-c1345)
								+ty1*(c34-c1345)
								+tz1*(r43*c34-c1345))*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
				  )*tmp3
						+(tx1+ty1+tz1)*c1345*tmp2*u[ku4+ku3+ku2+4]);
			d[ks4+ks3+ks2T1+4]=dt*2.0
				*(tx1*(r43*c34-c1345)
						+ty1*(c34-c1345)
						+tz1*(c34-c1345))*tmp2*u[ku4+ku3+ku2+1];
			d[ks4+ks3+ks2T2+4]=dt*2.0
				*(tx1*(c34-c1345)
						+ty1*(r43*c34-c1345)
						+tz1*(c34-c1345))*tmp2*u[ku4+ku3+ku2+2];
			d[ks4+ks3+ks2T3+4]=dt*2.0
				*(tx1*(c34-c1345)
						+ty1*(c34-c1345)
						+tz1*(r43*c34-c1345))*tmp2*u[ku4+ku3+ku2+3];
			d[ks4+ks3+ks2T4+4]=1.0
				+dt*2.0*(tx1+ty1+tz1)*c1345*tmp1
				+dt*2.0*(tx1*dx5+ty1*dy5+tz1*dz5);
			/*
			 * ---------------------------------------------------------------------
			 * form the first block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2+i+1];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			a[ks4+ks3]=-dt*tx1*dx1;
			a[ks4+ks3+ks2T1]=dt*tx2;
			a[ks4+ks3+ks2T2]=0.0;
			a[ks4+ks3+ks2T3]=0.0;
			a[ks4+ks3+ks2T4]=0.0;
			a[ks4+ks3+1]=dt*tx2
				*(-(u[ku4+ku3+ku2P+1]*tmp1)*(u[ku4+ku3+ku2P+1]*tmp1)
						+C2*qs[kq3+kq2+i+1]*tmp1)
				-dt*tx1*(-r43*c34*tmp2*u[ku4+ku3+ku2P+1]);
			a[ks4+ks3+ks2T1+1]=dt*tx2
				*((2.0-C2)*(u[ku4+ku3+ku2P+1]*tmp1))
				-dt*tx1*(r43*c34*tmp1)
				-dt*tx1*dx2;
			a[ks4+ks3+ks2T2+1]=dt*tx2
				*(-C2*(u[ku4+ku3+ku2P+2]*tmp1));
			a[ks4+ks3+ks2T3+1]=dt*tx2
				*(-C2*(u[ku4+ku3+ku2P+3]*tmp1));
			a[ks4+ks3+ks2T4+1]=dt*tx2*C2;
			a[ks4+ks3+2]=dt*tx2
				*(-(u[ku4+ku3+ku2P+1]*u[ku4+ku3+ku2P+2])*tmp2)
				-dt*tx1*(-c34*tmp2*u[ku4+ku3+ku2P+2]);
			a[ks4+ks3+ks2T1+2]=dt*tx2*(u[ku4+ku3+ku2P+2]*tmp1);
			a[ks4+ks3+ks2T2+2]=dt*tx2*(u[ku4+ku3+ku2P+1]*tmp1)
				-dt*tx1*(c34*tmp1)
				-dt*tx1*dx3;
			a[ks4+ks3+ks2T3+2]=0.0;
			a[ks4+ks3+ks2T4+2]=0.0;
			a[ks4+ks3+3]=dt*tx2
				*(-(u[ku4+ku3+ku2P+1]*u[ku4+ku3+ku2P+3])*tmp2)
				-dt*tx1*(-c34*tmp2*u[ku4+ku3+ku2P+3]);
			a[ks4+ks3+ks2T1+3]=dt*tx2*(u[ku4+ku3+ku2P+3]*tmp1);
			a[ks4+ks3+ks2T2+3]=0.0;
			a[ks4+ks3+ks2T3+3]=dt*tx2*(u[ku4+ku3+ku2P+1]*tmp1)
				-dt*tx1*(c34*tmp1)
				-dt*tx1*dx4;
			a[ks4+ks3+ks2T4+3]=0.0;
			a[ks4+ks3+4]=dt*tx2
				*((C2*2.0*qs[kq3+kq2+i+1]
							-C1*u[ku4+ku3+ku2P+4])
						*(u[ku4+ku3+ku2P+1]*tmp2))
				-dt*tx1
				*(-(r43*c34-c1345)*tmp3*(u[ku4+ku3+ku2P+1]*u[ku4+ku3+ku2P+1])
						-(c34-c1345)*tmp3*(u[ku4+ku3+ku2P+2]*u[ku4+ku3+ku2P+2])
						-(c34-c1345)*tmp3*( u[ku4+ku3+ku2P+3]*u[ku4+ku3+ku2P+3])
						-c1345*tmp2*u[ku4+ku3+ku2P+4]);
			a[ks4+ks3+ks2T1+4]=dt*tx2
				*(C1*(u[ku4+ku3+ku2P+4]*tmp1)
						-C2
						*(u[ku4+ku3+ku2P+1]*u[ku4+ku3+ku2P+1]*tmp2
							+qs[kq3+kq2+i+1]*tmp1))
				-dt*tx1
				*(r43*c34-c1345)*tmp2*u[ku4+ku3+ku2P+1];
			a[ks4+ks3+ks2T2+4]=dt*tx2
				*(-C2*(u[ku4+ku3+ku2P+2]*u[ku4+ku3+ku2P+1])*tmp2)
				-dt*tx1
				*(c34-c1345)*tmp2*u[ku4+ku3+ku2P+2];
			a[ks4+ks3+ks2T3+4]=dt*tx2
				*(-C2*(u[ku4+ku3+ku2P+3]*u[ku4+ku3+ku2P+1])*tmp2)
				-dt*tx1
				*(c34-c1345)*tmp2*u[ku4+ku3+ku2P+3];
			a[ks4+ks3+ks2T4+4]=dt*tx2
				*(C1*(u[ku4+ku3+ku2P+1]*tmp1))
				-dt*tx1*c1345*tmp1
				-dt*tx1*dx5;
			/*
			 * ---------------------------------------------------------------------
			 * form the second block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3+kq2P+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			b[ks4+ks3]=-dt*ty1*dy1;
			b[ks4+ks3+ks2T1]=0.0;
			b[ks4+ks3+ks2T2]=dt*ty2;
			b[ks4+ks3+ks2T3]=0.0;
			b[ks4+ks3+ks2T4]=0.0;
			b[ks4+ks3+1]=dt*ty2
				*(-(u[ku4+ku3P+ku2+1]*u[ku4+ku3P+ku2+2])*tmp2)
				-dt*ty1*(-c34*tmp2*u[ku4+ku3P+ku2+1]);
			b[ks4+ks3+ks2T1+1]=dt*ty2*(u[ku4+ku3P+ku2+2]*tmp1)
				-dt*ty1*(c34*tmp1)
				-dt*ty1*dy2;
			b[ks4+ks3+ks2T2+1]=dt*ty2*(u[ku4+ku3P+ku2+1]*tmp1);
			b[ks4+ks3+ks2T3+1]=0.0;
			b[ks4+ks3+ks2T4+1]=0.0;
			b[ks4+ks3+2]=dt*ty2
				*(-(u[ku4+ku3P+ku2+2]*tmp1)*(u[ku4+ku3P+ku2+2]*tmp1)
						+C2*(qs[kq3+kq2P+i]*tmp1))
				-dt*ty1*(-r43*c34*tmp2*u[ku4+ku3P+ku2+2]);
			b[ks4+ks3+ks2T1+2]=dt*ty2
				*(-C2*(u[ku4+ku3P+ku2+1]*tmp1));
			b[ks4+ks3+ks2T2+2]=dt*ty2*((2.0-C2)
					*(u[ku4+ku3P+ku2+2]*tmp1))
				-dt*ty1*(r43*c34*tmp1)
				-dt*ty1*dy3;
			b[ks4+ks3+ks2T3+2]=dt*ty2
				*(-C2*(u[ku4+ku3P+ku2+3]*tmp1));
			b[ks4+ks3+ks2T4+2]=dt*ty2*C2;
			b[ks4+ks3+3]=dt*ty2
				*(-(u[ku4+ku3P+ku2+2]*u[ku4+ku3P+ku2+3])*tmp2)
				-dt*ty1*(-c34*tmp2*u[ku4+ku3P+ku2+3]);
			b[ks4+ks3+ks2T1+3]=0.0;
			b[ks4+ks3+ks2T2+3]=dt*ty2*(u[ku4+ku3P+ku2+3]*tmp1);
			b[ks4+ks3+ks2T3+3]=dt*ty2*(u[ku4+ku3P+ku2+2]*tmp1)
				-dt*ty1*(c34*tmp1)
				-dt*ty1*dy4;
			b[ks4+ks3+ks2T4+3]=0.0;
			b[ks4+ks3+4]=dt*ty2
				*((C2*2.0*qs[kq3+kq2P+i]
							-C1*u[ku4+ku3P+ku2+4])
						*(u[ku4+ku3P+ku2+2]*tmp2))
				-dt*ty1
				*(-(c34-c1345)*tmp3*(u[ku4+ku3P+ku2+1]*u[ku4+ku3P+ku2+1])
						-(r43*c34-c1345)*tmp3*(u[ku4+ku3P+ku2+2]*u[ku4+ku3P+ku2+2])
						-(c34-c1345)*tmp3*(u[ku4+ku3P+ku2+3]*u[ku4+ku3P+ku2+3])
						-c1345*tmp2*u[ku4+ku3P+ku2+4]);
			b[ks4+ks3+ks2T1+4]=dt*ty2
				*(-C2*(u[ku4+ku3P+ku2+1]*u[ku4+ku3P+ku2+2])*tmp2)
				-dt*ty1
				*(c34-c1345)*tmp2*u[ku4+ku3P+ku2+1];
			b[ks4+ks3+ks2T2+4]=dt*ty2
				*(C1*(u[ku4+ku3P+ku2+4]*tmp1)
						-C2 
						*(qs[kq3+kq2P+i]*tmp1
							+u[ku4+ku3P+ku2+2]*u[ku4+ku3P+ku2+2]*tmp2))
				-dt*ty1
				*(r43*c34-c1345)*tmp2*u[ku4+ku3P+ku2+2];
			b[ks4+ks3+ks2T3+4]=dt*ty2
				*(-C2*(u[ku4+ku3P+ku2+2]*u[ku4+ku3P+ku2+3])*tmp2)
				-dt*ty1*(c34-c1345)*tmp2*u[ku4+ku3P+ku2+3];
			b[ks4+ks3+ks2T4+4]=dt*ty2
				*(C1*(u[ku4+ku3P+ku2+2]*tmp1))
				-dt*ty1*c1345*tmp1
				-dt*ty1*dy5;
			/*
			 * ---------------------------------------------------------------------
			 * form the third block sub-diagonal
			 * ---------------------------------------------------------------------
			 */
			tmp1=rho_i[kq3P+kq2+i];
			tmp2=tmp1*tmp1;
			tmp3=tmp1*tmp2;
			c[ks4+ks3]=-dt*tz1*dz1;
			c[ks4+ks3+ks2T1]=0.0;
			c[ks4+ks3+ks2T2]=0.0;
			c[ks4+ks3+ks2T3]=dt*tz2;
			c[ks4+ks3+ks2T4]=0.0;
			c[ks4+ks3+1]=dt*tz2
				*(-(u[ku4P+ku3+ku2+1]*u[ku4P+ku3+ku2+3])*tmp2)
				-dt*tz1*(-c34*tmp2*u[ku4P+ku3+ku2+1]);
			c[ks4+ks3+ks2T1+1]=dt*tz2*(u[ku4P+ku3+ku2+3]*tmp1)
				-dt*tz1*c34*tmp1
				-dt*tz1*dz2;
			c[ks4+ks3+ks2T2+1]=0.0;
			c[ks4+ks3+ks2T3+1]=dt*tz2*(u[ku4P+ku3+ku2+1]*tmp1);
			c[ks4+ks3+ks2T4+1]=0.0;
			c[ks4+ks3+2]=dt*tz2
				*(-(u[ku4P+ku3+ku2+2]*u[ku4P+ku3+ku2+3])*tmp2)
				-dt*tz1*(-c34*tmp2*u[ku4P+ku3+ku2+2]);
			c[ks4+ks3+ks2T1+2]=0.0;
			c[ks4+ks3+ks2T2+2]=dt*tz2*(u[ku4P+ku3+ku2+3]*tmp1)
				-dt*tz1*(c34*tmp1)
				-dt*tz1*dz3;
			c[ks4+ks3+ks2T3+2]=dt*tz2*(u[ku4P+ku3+ku2+2]*tmp1);
			c[ks4+ks3+ks2T4+2]=0.0;
			c[ks4+ks3+3]=dt*tz2
				*(-(u[ku4P+ku3+ku2+3]*tmp1)*(u[ku4P+ku3+ku2+3]*tmp1)
						+C2*(qs[kq3P+kq2+i]*tmp1))
				-dt*tz1*(-r43*c34*tmp2*u[ku4P+ku3+ku2+3]);
			c[ks4+ks3+ks2T1+3]=dt*tz2
				*(-C2*(u[ku4P+ku3+ku2+1]*tmp1));
			c[ks4+ks3+ks2T2+3]=dt*tz2
				*(-C2*(u[ku4P+ku3+ku2+2]*tmp1));
			c[ks4+ks3+ks2T3+3]=dt*tz2*(2.0-C2)
				*(u[ku4P+ku3+ku2+3]*tmp1)
				-dt*tz1*(r43*c34*tmp1)
				-dt*tz1*dz4;
			c[ks4+ks3+ks2T4+3]=dt*tz2*C2;
			c[ks4+ks3+4]=dt*tz2
				*((C2*2.0*qs[kq3P+kq2+i]
							-C1*u[ku4P+ku3+ku2+4])
						*(u[ku4P+ku3+ku2+3]*tmp2))
				-dt*tz1
				*(-(c34-c1345)*tmp3*(u[ku4P+ku3+ku2+1]*u[ku4P+ku3+ku2+1])
						-(c34-c1345)*tmp3*(u[ku4P+ku3+ku2+2]*u[ku4P+ku3+ku2+2])
						-(r43*c34-c1345)*tmp3*(u[ku4P+ku3+ku2+3]*u[ku4P+ku3+ku2+3])
						-c1345*tmp2*u[ku4P+ku3+ku2+4]);
			c[ks4+ks3+ks2T1+4]=dt*tz2
				*(-C2*(u[ku4P+ku3+ku2+1]*u[ku4P+ku3+ku2+3])*tmp2)
				-dt*tz1*(c34-c1345)*tmp2*u[ku4P+ku3+ku2+1];
			c[ks4+ks3+ks2T2+4]=dt*tz2
				*(-C2*(u[ku4P+ku3+ku2+2]*u[ku4P+ku3+ku2+3])*tmp2)
				-dt*tz1*(c34-c1345)*tmp2*u[ku4P+ku3+ku2+2];
			c[ks4+ks3+ks2T3+4]=dt*tz2
				*(C1*(u[ku4P+ku3+ku2+4]*tmp1)
						-C2
						*(qs[kq3P+kq2+i]*tmp1
							+u[ku4P+ku3+ku2+3]*u[ku4P+ku3+ku2+3]*tmp2))
				-dt*tz1*(r43*c34-c1345)*tmp2*u[ku4P+ku3+ku2+3];
			c[ks4+ks3+ks2T4+4]=dt*tz2
				*(C1*(u[ku4P+ku3+ku2+3]*tmp1))
				-dt*tz1*c1345*tmp1
				-dt*tz1*dz5;
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * to compute the l2-norm of vector v.
 * ---------------------------------------------------------------------
 * to improve cache performance, second two dimensions padded by 1 
 * for even number sizes only.  Only needed in v.
 * ---------------------------------------------------------------------
 */
void l2norm(int nx0,
		int ny0,
		int nz0,
		int ist,
		int iend,
		int jst,
		int jend,
		std::vector<double> &v,
		std::vector<double> &sum){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	std::fill(sum.begin(), sum.end(), 0.0);
	Returnable result = std::transform_reduce(policy, iter.front()+1, iter.front()+1+nz0-1, Returnable(), std::plus<Returnable>(), [&] (int k) {
		std::vector<double> sum(5, 0.0);
		int ku4 = k*ku4_const;
		for(int j=jst; j<jend; j++){
			int ku3 = j*ku3_const;
			for(int i=ist; i<iend; i++){
				int ku2 = i*ku2_const;
				sum[0] += v[ku4+ku3+ku2]*v[ku4+ku3+ku2];				
				sum[1] += v[ku4+ku3+ku2+1]*v[ku4+ku3+ku2+1];
				sum[2] += v[ku4+ku3+ku2+2]*v[ku4+ku3+ku2+2];
				sum[3] += v[ku4+ku3+ku2+3]*v[ku4+ku3+ku2+3];
				sum[4] += v[ku4+ku3+ku2+4]*v[ku4+ku3+ku2+4];
			}
		}

		return Returnable(sum);
	});

	sum = result.sum;

	std::transform(sum.begin(), sum.end(), sum.begin(), [&](double x){return sqrt(x/((nx0-2)*(ny0-2)*(nz0-2)));});
}

void pintgr(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
  	int k;
	int ibeg, ifin, ifin1;
	int jbeg, jfin, jfin1;
	std::vector<double> phi1((ISIZ3+2)*(ISIZ2+2), 0.0);
	std::vector<double> phi2((ISIZ3+2)*(ISIZ2+2), 0.0);
	int phi_const = ISIZ2+2, phi=0, phip=0;
	double frc1, frc2, frc3;
	/*
	 * ---------------------------------------------------------------------
	 * set up the sub-domains for integeration in each processor
	 * ---------------------------------------------------------------------
	 */
	ibeg=ii1;
	ifin=ii2;
	jbeg=ji1;
	jfin=ji2;
	ifin1=ifin-1;
	jfin1=jfin-1;
	/*
	 * ---------------------------------------------------------------------
	 * initialize
	 * ---------------------------------------------------------------------
	 */
	std::fill(phi1.begin(), phi1.end(), 0.0);
	std::fill(phi2.begin(), phi2.end(), 0.0);
	for(int j=jbeg; j<jfin; j++){
		ku3 = j*ku3_const;
		phi = j*phi_const;
		for(int i=ibeg; i<ifin; i++){
			ku2 = i*ku2_const;
			k=ki1;
			ku4 = k*ku4_const;
			phi1[phi+i]=C2*(u[ku4+ku3+ku2+4]
					-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
						+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
						+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					/u[ku4+ku3+ku2]);
			k=ki2-1;
			ku4 = k*ku4_const;
			phi2[phi+i]=C2*(u[ku4+ku3+ku2+4]
					-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
						+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
						+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					/u[ku4+ku3+ku2]);
		}
	}
	frc1=0.0;
	for(int j=jbeg; j<jfin1; j++){
		phi = j*phi_const;
		phip = (j+1)*phi_const;
		for(int i=ibeg; i<ifin1; i++){
			frc1=frc1+(phi1[phi+i]
					+phi1[phi+i+1]
					+phi1[phip+i]
					+phi1[phip+i+1]
					+phi2[phi+i]
					+phi2[phi+i+1]
					+phi2[phip+i]
					+phi2[phip+i+1]);
		}
	}
	frc1=dxi*deta*frc1;
	/*
	 * ---------------------------------------------------------------------
	 * initialize
	 * ---------------------------------------------------------------------
	 */
	std::fill(phi1.begin(), phi1.end(), 0.0);
	std::fill(phi2.begin(), phi2.end(), 0.0);
	if(jbeg==ji1){
		ku3 = jbeg*ku3_const;
		for(int k=ki1; k<ki2; k++){
			ku4 = k*ku4_const;
			phi = k*phi_const;
			for(int i=ibeg; i<ifin; i++){
				ku2 = i*ku2_const;
				phi1[phi+i]=C2*(u[ku4+ku3+ku2+4]
						-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
							+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
							+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
						/u[ku4+ku3+ku2]);
			}
		}
	}
	if(jfin==ji2){
		ku3 = (jfin-1)*ku3_const;
		for(int k=ki1; k<ki2; k++){
			ku4 = k*ku4_const;
			phi = k*phi_const;
			for(int i=ibeg; i<ifin; i++){
				ku2 = i*ku2_const;
				phi2[phi+i]=C2*(u[ku4+ku3+ku2+4]
						-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
							+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
							+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
						/u[ku4+ku3+ku2]);
			}
		}
	}
	frc2=0.0;
	for(int k=ki1; k<ki2-1; k++){
		phi = k*phi_const;
		phip = (k+1)*phi_const;
		for(int i=ibeg; i<ifin1; i++){
			frc2=frc2+(phi1[phi+i]
					+phi1[phi+i+1]
					+phi1[phip+i]
					+phi1[phip+i+1]
					+phi2[phi+i]
					+phi2[phi+i+1]
					+phi2[phip+i]
					+phi2[phip+i+1]);
		}
	}
	frc2=dxi*dzeta*frc2;
	/*
	 * ---------------------------------------------------------------------
	 * initialize
	 * ---------------------------------------------------------------------
	 */
	std::fill(phi1.begin(), phi1.end(), 0.0);
	std::fill(phi2.begin(), phi2.end(), 0.0);
	if(ibeg==ii1){
		ku2 = ibeg*ku2_const;
		for(int k=ki1; k<ki2; k++){
			ku4 = k*ku4_const;
			phi = k*phi_const;
			for(int j=jbeg; j<jfin; j++){
				ku3 = j*ku3_const;
				phi1[phi+j]=C2*(u[ku4+ku3+ku2+4]
						-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
							+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
							+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
						/u[ku4+ku3+ku2]);
			}
		}
	}
	if(ifin==ii2){
		ku2 = (ifin-1)*ku2_const;
		for(int k=ki1; k<ki2; k++){
			ku4 = k*ku4_const;
			phi = k*phi_const;
			for(int j=jbeg; j<jfin; j++){
				ku3 = j*ku3_const;
				phi2[phi+j]=C2*(u[ku4+ku3+ku2+4]
						-0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
							+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
							+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
						/u[ku4+ku3+ku2]);
			}
		}
	}
	frc3=0.0;
	for(int k=ki1; k<ki2-1; k++){
		phi = k*phi_const;
		phip = (k+1)*phi_const;
		for(int j=jbeg; j<jfin1; j++){
			frc3=frc3+(phi1[phi+j]
					+phi1[phi+j+1]
					+phi1[phip+j]
					+phi1[phip+j+1]
					+phi2[phi+j]
					+phi2[phi+j+1]
					+phi2[phip+j]
					+phi2[phip+j+1]);
		}
	}
	frc3=deta*dzeta*frc3;
	frc=0.25*(frc1+frc2+frc3);
}

void read_input(){
	/*
	 * ---------------------------------------------------------------------
	 * if input file does not exist, it uses defaults
	 * ipr = 1 for detailed progress output
	 * inorm = how often the norm is printed (once every inorm iterations)
	 * itmax = number of pseudo time steps
	 * dt = time step
	 * omega 1 over-relaxation factor for SSOR
	 * tolrsd = steady state residual tolerance levels
	 * nx, ny, nz = number of grid points in x, y, z directions
	 * ---------------------------------------------------------------------
	 */	
	std::fstream fp;
  fp.open("inputlu.data", std::ios::in);
	if(fp.is_open()){
		std::cout
      << "Reading from input file inputlu.data"
      <<std::endl;
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		fp >> ipr >> inorm;
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		while(fp.get() != '\n');
    fp >> itmax;
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		fp >> dt;
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		while(fp.get() != '\n');
    fp >> omega;
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		fp
      >> tolrsd[0]
      >> tolrsd[1]
      >> tolrsd[2]
      >> tolrsd[3]
      >> tolrsd[4];
		while(fp.get() != '\n');
		while(fp.get() != '\n');
		fp >> nx0 >> ny0 >> nz0;
		fp.close();

	}else{
		ipr=IPR_DEFAULT;
		inorm=INORM_DEFAULT;
		itmax=ITMAX_DEFAULT;
		dt=DT_DEFAULT;
		omega=OMEGA_DEFAULT;
		tolrsd[0]=TOLRSD1_DEF;
		tolrsd[1]=TOLRSD2_DEF;
		tolrsd[2]=TOLRSD3_DEF;
		tolrsd[3]=TOLRSD4_DEF;
		tolrsd[4]=TOLRSD5_DEF;
		nx0=ISIZ1;
		ny0=ISIZ2;
		nz0=ISIZ3;
	}
	/*
	 * ---------------------------------------------------------------------
	 * check problem size
	 * ---------------------------------------------------------------------
	 */
	if((nx0<4)||(ny0<4)||(nz0<4)){
		std::cout
      << "     PROBLEM SIZE IS TOO SMALL - "
      << std::endl
			<< "     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5"
      << std::endl;
		exit(EXIT_FAILURE);
	}
	if((nx0>ISIZ1)||(ny0>ISIZ2)||(nz0>ISIZ3)){
		std::cout
      << "     PROBLEM SIZE IS TOO LARGE - "
      << std::endl
			<< "     NX, NY AND NZ SHOULD BE EQUAL TO "
      << std::endl
			<< "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY"
      << std::endl;
		exit(EXIT_FAILURE);
	}
  std::cout
    << std::endl
    << std::endl;
  
  std::cout
    << " NAS Parallel Benchmarks 4.1 Serial C++ version - LU Benchmark";
  
  std::cout
    << std::endl
    << std::endl;

  std::cout
    << " Size: "
    << std::setw(4) << nx0 << "x"
    << std::setw(4) << ny0 << "x"
    << std::setw(4) << nz0
    << std::endl;

  std::cout
    << " Iterations: "
    << std::setw(4) << itmax
    << std::endl;
  
  std::cout
    << std::endl;
}

/*
 * ---------------------------------------------------------------------
 * compute the right hand sides
 * ---------------------------------------------------------------------
 */
void rhs(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	double q;
	double u21, u31, u41;
	double u21i, u31i, u41i, u51i;
	double u21j, u31j, u41j, u51j;
	double u21k, u31k, u41k, u51k;
	double u21im1, u31im1, u41im1, u51im1;
	double u21jm1, u31jm1, u41jm1, u51jm1;
	double u21km1, u31km1, u41km1, u51km1;
	if(timeron){timer_start(T_RHS);}
	std::for_each_n(policy, iter.front(), nz, [&] (int k){
		double tmp;
		int ku2, ku3, ku4, kq2, kq3;

		kq3 = k*kq3_const;
		ku4 = k*ku4_const;
		for(int j=0; j<ny; j++){
			kq2 = j*kq2_const;
			ku3 = j*ku3_const;
			for(int i=0; i<nx; i++){
				ku2 = i*ku2_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]=-frct[ku4+ku3+ku2+m];
				}
				tmp=1.0/u[ku4+ku3+ku2];
				rho_i[kq3+kq2+i]=tmp;
				qs[kq3+kq2+i]=0.50*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]
						+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]
						+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					*tmp;
			}
		}
	});
	if(timeron){timer_start(T_RHSX);}
	/*
	 * ---------------------------------------------------------------------
	 * xi-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, nz-2, [&] (int k){
		double q;
		double tmp;
		double u21;
		double u21i, u31i, u41i, u51i;
		double u21im1, u31im1, u41im1, u51im1;
		int ku2, ku3, ku4, kq2, kq3;
		int kfl2, kfl2P, kfl2M;
		int ku2P, ku2P2, ku2M, ku2M2;
		std::vector<double> flux(ISIZ1*5, 0);

		kq3 = k*kq3_const;
		ku4 = k*ku4_const;
		for(int j=jst; j<jend; j++){
			kq2 = j*kq2_const;
			ku3 = j*ku3_const;
			for(int i=0; i<nx; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				flux[kfl2]=u[ku4+ku3+ku2+1];
				u21=u[ku4+ku3+ku2+1]*rho_i[kq3+kq2+i];
				q=qs[kq3+kq2+i];
				flux[kfl2+1]=u[ku4+ku3+ku2+1]*u21+C2*(u[ku4+ku3+ku2+4]-q);
				flux[kfl2+2]=u[ku4+ku3+ku2+2]*u21;
				flux[kfl2+3]=u[ku4+ku3+ku2+3]*u21;
				flux[kfl2+4]=(C1*u[ku4+ku3+ku2+4]-C2*q)*u21;
			}
			for(int i=ist; i<iend; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]-=tx2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int i=ist; i<nx; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				tmp=rho_i[kq3+kq2+i];
				u21i=tmp*u[ku4+ku3+ku2+1];
				u31i=tmp*u[ku4+ku3+ku2+2];
				u41i=tmp*u[ku4+ku3+ku2+3];
				u51i=tmp*u[ku4+ku3+ku2+4];
				ku2 -= ku2_const;
				tmp=rho_i[kq3+kq2+i-1];
				u21im1=tmp*u[ku4+ku3+ku2+1];
				u31im1=tmp*u[ku4+ku3+ku2+2];
				u41im1=tmp*u[ku4+ku3+ku2+3];
				u51im1=tmp*u[ku4+ku3+ku2+4];
				flux[kfl2+1]=(4.0/3.0)*tx3*(u21i-u21im1);
				flux[kfl2+2]=tx3*(u31i-u31im1);
				flux[kfl2+3]=tx3*(u41i-u41im1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*tx3*((u21i*u21i+u31i*u31i+u41i*u41i)
							-(u21im1*u21im1+u31im1*u31im1+u41im1*u41im1))
					+(1.0/6.0)
					*tx3*(u21i*u21i-u21im1*u21im1)
					+C1*C5*tx3*(u51i-u51im1);
			}
			for(int i=ist; i<iend; i++){
				kfl2 = i*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku2 = i*ku2_const;
				int ku2P = ku2+ku2_const;
				int ku2M = ku2-ku2_const;
				rsd[ku4+ku3+ku2]+=dx1*tx1*(u[ku4+ku3+ku2M]
							-2.0* u[ku4+ku3+ku2]
							+u[ku4+ku3+ku2P]);
				rsd[ku4+ku3+ku2+1]+=tx3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dx2*tx1*(u[ku4+ku3+ku2M+1]
							-2.0*u[ku4+ku3+ku2+1]
							+u[ku4+ku3+ku2P+1]);
				rsd[ku4+ku3+ku2+2]+=tx3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
							+dx3*tx1*(u[ku4+ku3+ku2M+2]
							-2.0*u[ku4+ku3+ku2+2]
							+u[ku4+ku3+ku2P+2]);
				rsd[ku4+ku3+ku2+3]+=tx3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
							+dx4*tx1*(u[ku4+ku3+ku2M+3]
							-2.0*u[ku4+ku3+ku2+3]
							+u[ku4+ku3+ku2P+3]);
				rsd[ku4+ku3+ku2+4]+=tx3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
							+dx5*tx1*(u[ku4+ku3+ku2M+4]
							-2.0*u[ku4+ku3+ku2+4]
							+u[ku4+ku3+ku2P+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			ku2 = 2*ku2_const;
			int ku2M = ku2-ku2_const;
			int ku2P = ku2+ku2_const;
			int ku2P2 = ku2P+ku2_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3+ku2M+m]-=dssp*(+5.0*u[ku4+ku3+ku2M+m]
							-4.0*u[ku4+ku3+ku2+m]
							+u[ku4+ku3+ku2P+m]);
				rsd[ku4+ku3+ku2+m]-=dssp*(-4.0*u[ku4+ku3+ku2M+m]
							+6.0*u[ku4+ku3+ku2+m]
							-4.0*u[ku4+ku3+ku2P+m]
							+u[ku4+ku3+ku2P2+m]);
			}
			for(int i=3; i<nx-3; i++){
				ku2 = i*ku2_const;
				ku2M = ku2-ku2_const;
				int ku2M2 = ku2M-ku2_const;
				ku2P = ku2+ku2_const;
				ku2P2 = ku2P+ku2_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]-=dssp*(u[ku4+ku3+ku2M2+m]
								-4.0*u[ku4+ku3+ku2M+m]
								+6.0*u[ku4+ku3+ku2+m]
								-4.0*u[ku4+ku3+ku2P+m]
								+u[ku4+ku3+ku2P2+m]);
				}
			}
			ku2 = (nx-4)*ku2_const;
			ku2M = ku2-ku2_const;
			ku2P = ku2+ku2_const;
			ku2P2 = ku2P+ku2_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3+ku2P+m]-=dssp*(u[ku4+ku3+ku2M+m]
							-4.0*u[ku4+ku3+ku2+m]
							+6.0*u[ku4+ku3+ku2P+m]
							-4.0*u[ku4+ku3+ku2P2+m]);
				rsd[ku4+ku3+ku2P2+m]-=dssp*(u[ku4+ku3+ku2+m]
							-4.0*u[ku4+ku3+ku2P+m]
							+5.0*u[ku4+ku3+ku2P2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSX);}
	if(timeron){timer_start(T_RHSY);}
	/*
	 * ---------------------------------------------------------------------
	 * eta-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, nz-2, [&] (int k){
		double q;
		double tmp;
		double u21, u31, u41;
		double u21j, u31j, u41j, u51j;
		double u21jm1, u31jm1, u41jm1, u51jm1;
		int ku2, ku3, ku4, kq2, kq3;
		int kfl2, kfl2P, kfl2M;
		int ku3P, ku3P2, ku3M, ku3M2;
		std::vector<double> flux(ISIZ1*5, 0);

		kq3 = k*kq3_const;
		ku4 = k*ku4_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int j=0; j<ny; j++){
				kq2 = j*kq2_const;
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				flux[kfl2]=u[ku4+ku3+ku2+2];
				u31=u[ku4+ku3+ku2+2]*rho_i[kq3+kq2+i];
				q=qs[kq3+kq2+i];
				flux[kfl2+1]=u[ku4+ku3+ku2+1]*u31;
				flux[kfl2+2]=u[ku4+ku3+ku2+2]*u31+C2*(u[ku4+ku3+ku2+4]-q);
				flux[kfl2+3]=u[ku4+ku3+ku2+3]*u31;
				flux[kfl2+4]=(C1*u[ku4+ku3+ku2+4]-C2*q)*u31;
			}
			for(int j=jst; j<jend; j++){
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]-=ty2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int j=jst; j<ny; j++){
				kq2 = j*kq2_const;
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				tmp=rho_i[kq3+kq2+i];
				u21j=tmp*u[ku4+ku3+ku2+1];
				u31j=tmp*u[ku4+ku3+ku2+2];
				u41j=tmp*u[ku4+ku3+ku2+3];
				u51j=tmp*u[ku4+ku3+ku2+4];
				ku3 -= ku3_const;
				kq2 -= kq2_const;
				tmp=rho_i[kq3+kq2+i];
				u21jm1=tmp*u[ku4+ku3+ku2+1];
				u31jm1=tmp*u[ku4+ku3+ku2+2];
				u41jm1=tmp*u[ku4+ku3+ku2+3];
				u51jm1=tmp*u[ku4+ku3+ku2+4];
				flux[kfl2+1]=ty3*(u21j-u21jm1);
				flux[kfl2+2]=(4.0/3.0)*ty3*(u31j-u31jm1);
				flux[kfl2+3]=ty3*(u41j-u41jm1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*ty3*((u21j*u21j+u31j*u31j+u41j*u41j)
							-(u21jm1*u21jm1+u31jm1*u31jm1+u41jm1*u41jm1))
					+(1.0/6.0)
					*ty3*(u31j*u31j-u31jm1*u31jm1)
					+C1*C5*ty3*(u51j-u51jm1);
			}
			for(int j=jst; j<jend; j++){
				kfl2 = j*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku3 = j*ku3_const;
				int ku3P = ku3+ku3_const;
				int ku3M = ku3-ku3_const;
				rsd[ku4+ku3+ku2]+=dy1*ty1*(u[ku4+ku3M+ku2]
							-2.0*u[ku4+ku3+ku2]
							+u[ku4+ku3P+ku2]);
				rsd[ku4+ku3+ku2+1]+=ty3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dy2*ty1*(u[ku4+ku3M+ku2+1]
							-2.0*u[ku4+ku3+ku2+1]
							+u[ku4+ku3P+ku2+1]);
				rsd[ku4+ku3+ku2+2]+=ty3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
					+dy3*ty1*(u[ku4+ku3M+ku2+2]
							-2.0*u[ku4+ku3+ku2+2]
							+u[ku4+ku3P+ku2+2]);
				rsd[ku4+ku3+ku2+3]+=ty3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
					+dy4*ty1*(u[ku4+ku3M+ku2+3]
							-2.0*u[ku4+ku3+ku2+3]
							+u[ku4+ku3P+ku2+3]);
				rsd[ku4+ku3+ku2+4]+=ty3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
					+dy5*ty1*(u[ku4+ku3M+ku2+4]
							-2.0*u[ku4+ku3+ku2+4]
							+u[ku4+ku3P+ku2+4]);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * fourth-order dissipation
		 * ---------------------------------------------------------------------
		 */
		ku3 = 2*ku3_const;
		ku3M = ku3-ku3_const;
		ku3P = ku3+ku3_const;
		ku3P2 = ku3P+ku3_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3M+ku2+m]-=dssp*(+5.0*u[ku4+ku3M+ku2+m]
							-4.0*u[ku4+ku3+ku2+m]
							+u[ku4+ku3P+ku2+m]);
				rsd[ku4+ku3+ku2+m]-=dssp*(-4.0*u[ku4+ku3M+ku2+m]
							+6.0*u[ku4+ku3+ku2+m]
							-4.0*u[ku4+ku3P+ku2+m]
							+u[ku4+ku3P2+ku2+m]);
			}
		}
		for(int j=3; j<ny-3; j++){
			ku3 = j*ku3_const;
			ku3M = ku3-ku3_const;
			ku3M2 = ku3M-ku3_const;
			ku3P = ku3+ku3_const;
			ku3P2 = ku3P+ku3_const;
			for(int i=ist; i<iend; i++){
				ku2 = i*ku2_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]-=dssp*(u[ku4+ku3M2+ku2+m]
								-4.0*u[ku4+ku3M+ku2+m]
								+6.0*u[ku4+ku3+ku2+m]
								-4.0*u[ku4+ku3P+ku2+m]
								+u[ku4+ku3P2+ku2+m]);
				}
			}
		}
		ku3 = (ny-4)*ku3_const;
		ku3M = ku3-ku3_const;
		ku3P = ku3+ku3_const;
		ku3P2 = ku3P+ku3_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3P+ku2+m]-=dssp*(u[ku4+ku3M+ku2+m]
							-4.0*u[ku4+ku3+ku2+m]
							+6.0*u[ku4+ku3P+ku2+m]
							-4.0*u[ku4+ku3P2+ku2+m]);
				rsd[ku4+ku3P2+ku2+m]-=dssp*(u[ku4+ku3+ku2+m]
							-4.0*u[ku4+ku3P+ku2+m]
							+5.0*u[ku4+ku3P2+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSY);}
	if(timeron){timer_start(T_RHSZ);}
	/*
	 * ---------------------------------------------------------------------
	 * zeta-direction flux differences
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+jst, jend-jst, [&] (int j){
		double q;
		double u21, u31, u41;
		double tmp;
		std::vector<double> utmp(ISIZ3*6, 0);
		std::vector<double> rtmp(ISIZ3*5, 0);
		double u21k, u31k, u41k, u51k;
		double u21km1, u31km1, u41km1, u51km1;
		int ku2, ku3, ku4, kq2, kq3;
		int kfl2, kfl2P, kfl2M;
		int ku4P, ku4M;
		std::vector<double> flux(ISIZ1*5, 0);

		kq2 = j*kq2_const;
		ku3 = j*ku3_const;
		for(int i=ist; i<iend; i++){
			ku2 = i*ku2_const;
			for(int k=0; k<nz; k++){
				kq3 = k*kq3_const;
				ku4 = k*ku4_const;
				int kt2 = k*6;
				utmp[kt2]=u[ku4+ku3+ku2];
				utmp[kt2+1]=u[ku4+ku3+ku2+1];
				utmp[kt2+2]=u[ku4+ku3+ku2+2];
				utmp[kt2+3]=u[ku4+ku3+ku2+3];
				utmp[kt2+4]=u[ku4+ku3+ku2+4];
				utmp[kt2+5]=rho_i[kq3+kq2+i];
			}
			for(int k=0; k<nz; k++){
				kq3 = k*kq3_const;
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				int kt2 = k*6;
				flux[kfl2]=utmp[kt2+3];
				u41=utmp[kt2+3]*utmp[kt2+5];
				q=qs[kq3+kq2+i];
				flux[kfl2+1]=utmp[kt2+1]*u41;
				flux[kfl2+2]=utmp[kt2+2]*u41;
				flux[kfl2+3]=utmp[kt2+3]*u41+C2*(utmp[kt2+4]-q);
				flux[kfl2+4]=(C1*utmp[kt2+4]-C2*q)*u41;
			}
			for(int k=1; k<nz-1; k++){
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				ku4 = k*ku4_const;
				for(int m=0; m<5; m++){
					rtmp[k*5+m]=rsd[ku4+ku3+ku2+m]
						-tz2*(flux[kfl2P+m]-flux[kfl2M+m]);
				}
			}
			for(int k=1; k<nz; k++){
				int kt2 = k*6;
				int kt2M = (k-1)*6;
				kfl2 = k*fl_const;
				tmp=utmp[kt2+5];
				u21k=tmp*utmp[kt2+1];
				u31k=tmp*utmp[kt2+2];
				u41k=tmp*utmp[kt2+3];
				u51k=tmp*utmp[kt2+4];
				tmp=utmp[kt2M+5];
				u21km1=tmp*utmp[kt2M+1];
				u31km1=tmp*utmp[kt2M+2];
				u41km1=tmp*utmp[kt2M+3];
				u51km1=tmp*utmp[kt2M+4];
				flux[kfl2+1]=tz3*(u21k-u21km1);
				flux[kfl2+2]=tz3*(u31k-u31km1);
				flux[kfl2+3]=(4.0/3.0)*tz3*(u41k-u41km1);
				flux[kfl2+4]=0.50*(1.0-C1*C5)
					*tz3*((u21k*u21k+u31k*u31k+u41k*u41k)
							-(u21km1*u21km1+u31km1*u31km1+u41km1*u41km1))
					+(1.0/6.0)
					*tz3*(u41k*u41k-u41km1*u41km1)
					+C1*C5*tz3*(u51k-u51km1);
			}
			for(int k=1; k<nz-1; k++){
				int kt2 = k*6;
				int ktr2 = k*5;
				int kt2M = (k-1)*6;
				int kt2P = (k+1)*6;
				kfl2 = k*fl_const;
				kfl2P = kfl2+fl_const;
				kfl2M = kfl2-fl_const;
				rtmp[ktr2]=rtmp[ktr2]
					+dz1*tz1*(utmp[kt2M]
							-2.0*utmp[kt2]
							+utmp[kt2P]);
				rtmp[ktr2+1]=rtmp[ktr2+1]
					+tz3*C3*C4*(flux[kfl2P+1]-flux[kfl2+1])
					+dz2*tz1*(utmp[kt2M+1]
							-2.0*utmp[kt2+1]
							+utmp[kt2P+1]);
				rtmp[ktr2+2]=rtmp[ktr2+2]
					+tz3*C3*C4*(flux[kfl2P+2]-flux[kfl2+2])
					+dz3*tz1*(utmp[kt2M+2]
							-2.0*utmp[kt2+2]
							+utmp[kt2P+2]);
				rtmp[ktr2+3]=rtmp[ktr2+3]
					+tz3*C3*C4*(flux[kfl2P+3]-flux[kfl2+3])
					+dz4*tz1*(utmp[kt2M+3]
							-2.0*utmp[kt2+3]
							+utmp[kt2P+3]);
				rtmp[ktr2+4]=rtmp[ktr2+4]
					+tz3*C3*C4*(flux[kfl2P+4]-flux[kfl2+4])
					+dz5*tz1*(utmp[kt2M+4]
							-2.0*utmp[kt2+4]
							+utmp[kt2P+4]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			ku4 = ku4_const;
			ku4P = ku4+ku4_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3+ku2+m]=rtmp[5+m]
					-dssp*(+5.0*utmp[6+m]
							-4.0*utmp[12+m]
							+utmp[18+m]);
				rsd[ku4P+ku3+ku2+m]=rtmp[10+m]
					-dssp*(-4.0*utmp[6+m]
							+6.0*utmp[12+m]
							-4.0*utmp[18+m]
							+utmp[24+m] );
			}
			for(int k=3; k<nz-3; k++){
				int kt2 = k*6;
				int kt2M = (k-1)*6;
				int kt2P = (k+1)*6;
				ku4 = k*ku4_const;
				for(int m=0; m<5; m++){
					rsd[ku4+ku3+ku2+m]=rtmp[k*5+m]
						-dssp*(utmp[kt2M-6+m]
								-4.0*utmp[kt2M+m]
								+6.0*utmp[kt2+m]
								-4.0*utmp[kt2P+m]
								+utmp[kt2P+6+m]);
				}
			}
			ku4 = (nz-3)*ku4_const;
			ku4P = ku4+ku4_const;
			for(int m=0; m<5; m++){
				rsd[ku4+ku3+ku2+m]=rtmp[(nz-3)*5+m]
					-dssp*(utmp[(nz-5)*6+m]
							-4.0*utmp[(nz-4)*6+m]
							+6.0*utmp[(nz-3)*6+m]
							-4.0*utmp[(nz-2)*6+m]);
				rsd[ku4P+ku3+ku2+m]=rtmp[(nz-2)*5+m]
					-dssp*(utmp[(nz-4)*6+m]
							-4.0*utmp[(nz-3)*6+m]
							+5.0*utmp[(nz-2)*6+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSZ);}
	if(timeron){timer_stop(T_RHS);}
}

/*
 * ---------------------------------------------------------------------
 * set the boundary values of dependent variables
 * ---------------------------------------------------------------------
 */
void setbv(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	std::vector<double> temp1(5), temp2(5);
	/*
	 * ---------------------------------------------------------------------
	 * set the dependent variable values along the top and bottom faces
	 * ---------------------------------------------------------------------
	 */
	ku4 = (nz-1)*ku4_const;
	std::for_each_n(policy, iter.front(), ny, [&](int j){
		std::vector<double> temp1(5), temp2(5);
		int ku3, ku2;
		
		ku3 = j*ku3_const;
		for(int i=0; i<nx; i++){
			ku2 = i*ku2_const;
			exact(i, j, 0, temp1);
			exact(i, j, nz-1, temp2);
			for(int m=0; m<5; m++){
				u[ku3+ku2+m]=temp1[m];
				u[ku4+ku3+ku2+m]=temp2[m];
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * set the dependent variable values along north and south faces
	 * ---------------------------------------------------------------------
	 */
	ku3 = (ny-1)*ku3_const;
	std::for_each_n(policy, iter.front(), nz, [&](int k){
		std::vector<double> temp1(5), temp2(5);
		int ku4, ku2;

		ku4 = k*ku4_const;
		for(int i=0; i<nx; i++){
			ku2 = i*ku2_const;
			exact(i, 0, k, temp1);
			exact(i, ny-1, k, temp2);
			for(int m=0; m<5; m++){
				u[ku4+ku2+m]=temp1[m];
				u[ku4+ku3+ku2+m]=temp2[m];
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * set the dependent variable values along east and west faces
	 * ---------------------------------------------------------------------
	 */
	ku2 = (nx-1)*ku2_const;
	std::for_each_n(policy, iter.front(), nz, [&](int k){
		std::vector<double> temp1(5), temp2(5);
		int ku4, ku3;

		ku4 = k*ku4_const;
		for(int j=0; j<ny; j++){
			ku3 = j*ku3_const;
			exact(0, j, k, temp1);
			exact(nx-1, j, k, temp2);
			for(int m=0; m<5; m++){
				u[ku4+ku3+m]=temp1[m];
				u[ku4+ku3+ku2+m]=temp2[m];
			}
		}
	});
}

void setcoeff(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 * set up coefficients
	 * ---------------------------------------------------------------------
	 */
	dxi=1.0/(nx0-1);
	deta=1.0/(ny0-1);
	dzeta=1.0/(nz0-1);
	tx1=1.0/(dxi*dxi);
	tx2=1.0/(2.0*dxi);
	tx3=1.0/dxi;
	ty1=1.0/(deta*deta);
	ty2=1.0/(2.0*deta);
	ty3=1.0/deta;
	tz1=1.0/(dzeta*dzeta);
	tz2=1.0/(2.0*dzeta);
	tz3=1.0/dzeta;
	/*
	 * ---------------------------------------------------------------------
	 * diffusion coefficients
	 * ---------------------------------------------------------------------
	 */
	dx1=0.75;
	dx2=dx1;
	dx3=dx1;
	dx4=dx1;
	dx5=dx1;
	dy1=0.75;
	dy2=dy1;
	dy3=dy1;
	dy4=dy1;
	dy5=dy1;
	dz1=1.00;
	dz2=dz1;
	dz3=dz1;
	dz4=dz1;
	dz5=dz1;
	/*
	 * ---------------------------------------------------------------------
	 * fourth difference dissipation
	 * ---------------------------------------------------------------------
	 */
	dssp=(std::max(std::max(dx1,dy1),dz1))/4.0;
	/*
	 * ---------------------------------------------------------------------
	 * coefficients of the exact solution to the first pde
	 * ---------------------------------------------------------------------
	 */
	ce[(0*5)]=2.0;
	ce[(1*5)]=0.0;
	ce[(2*5)]=0.0;
	ce[(3*5)]=4.0;
	ce[(4*5)]=5.0;
	ce[(5*5)]=3.0;
	ce[(6*5)]=5.0e-01;
	ce[(7*5)]=2.0e-02;
	ce[(8*5)]=1.0e-02;
	ce[(9*5)]=3.0e-02;
	ce[(10*5)]=5.0e-01;
	ce[(11*5)]=4.0e-01;
	ce[(12*5)]=3.0e-01;
	/*
	 * ---------------------------------------------------------------------
	 * coefficients of the exact solution to the second pde
	 * ---------------------------------------------------------------------
	 */
	ce[(0*5)+1]=1.0;
	ce[(1*5)+1]=0.0;
	ce[(2*5)+1]=0.0;
	ce[(3*5)+1]=0.0;
	ce[(4*5)+1]=1.0;
	ce[(5*5)+1]=2.0;
	ce[(6*5)+1]=3.0;
	ce[(7*5)+1]=1.0e-02;
	ce[(8*5)+1]=3.0e-02;
	ce[(9*5)+1]=2.0e-02;
	ce[(10*5)+1]=4.0e-01;
	ce[(11*5)+1]=3.0e-01;
	ce[(12*5)+1]=5.0e-01;
	/*
	 * ---------------------------------------------------------------------
	 * coefficients of the exact solution to the third pde
	 * ---------------------------------------------------------------------
	 */
	ce[(0*5)+2]=2.0;
	ce[(1*5)+2]=2.0;
	ce[(2*5)+2]=0.0;
	ce[(3*5)+2]=0.0;
	ce[(4*5)+2]=0.0;
	ce[(5*5)+2]=2.0;
	ce[(6*5)+2]=3.0;
	ce[(7*5)+2]=4.0e-02;
	ce[(8*5)+2]=3.0e-02;
	ce[(9*5)+2]=5.0e-02;
	ce[(10*5)+2]=3.0e-01;
	ce[(11*5)+2]=5.0e-01;
	ce[(12*5)+2]=4.0e-01;
	/*
	 * ---------------------------------------------------------------------
	 * coefficients of the exact solution to the fourth pde
	 * ---------------------------------------------------------------------
	 */
	ce[(0*5)+3]=2.0;
	ce[(1*5)+3]=2.0;
	ce[(2*5)+3]=0.0;
	ce[(3*5)+3]=0.0;
	ce[(4*5)+3]=0.0;
	ce[(5*5)+3]=2.0;
	ce[(6*5)+3]=3.0;
	ce[(7*5)+3]=3.0e-02;
	ce[(8*5)+3]=5.0e-02;
	ce[(9*5)+3]=4.0e-02;
	ce[(10*5)+3]=2.0e-01;
	ce[(11*5)+3]=1.0e-01;
	ce[(12*5)+3]=3.0e-01;
	/*
	 * ---------------------------------------------------------------------
	 * coefficients of the exact solution to the fifth pde
	 * ---------------------------------------------------------------------
	 */
	ce[(0*5)+4]=5.0;
	ce[(1*5)+4]=4.0;
	ce[(2*5)+4]=3.0;
	ce[(3*5)+4]=2.0;
	ce[(4*5)+4]=1.0e-01;
	ce[(5*5)+4]=4.0e-01;
	ce[(6*5)+4]=3.0e-01;
	ce[(7*5)+4]=5.0e-02;
	ce[(8*5)+4]=4.0e-02;
	ce[(9*5)+4]=3.0e-02;
	ce[(10*5)+4]=1.0e-01;
	ce[(11*5)+4]=3.0e-01;
	ce[(12*5)+4]=2.0e-01;
}

/*
 * ---------------------------------------------------------------------
 * set the initial values of independent variables based on tri-linear
 * interpolation of boundary values in the computational space.
 * ---------------------------------------------------------------------
 */
void setiv(){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, nz-2, [&](int k){
		int ku2, ku3, ku4;
		double xi, eta, zeta;
		double pxi, peta, pzeta;
		std::vector<double> ue_1jk(5), ue_nx0jk(5), ue_i1k(5);
		std::vector<double> ue_iny0k(5), ue_ij1(5), ue_ijnz(5);
		ku4 = k*ku4_const;
		zeta=((double)k)/(nz-1);
		for(int j=1; j<ny-1; j++){
			ku3 = j*ku3_const;
			eta=((double)j)/(ny0-1);
			for(int i=1; i<nx-1; i++){
				ku2 = i*ku2_const;
				xi=((double)i)/(nx0-1);
				exact(0, j, k, ue_1jk);
				exact(nx0-1, j, k, ue_nx0jk);
				exact(i, 0, k, ue_i1k);
				exact(i, ny0-1, k, ue_iny0k);
				exact(i, j, 0, ue_ij1);
				exact(i, j, nz-1, ue_ijnz);
				for(int m=0; m<5; m++){
					pxi=(1.0-xi)*ue_1jk[m]
						+xi*ue_nx0jk[m];
					peta=(1.0-eta)*ue_i1k[m]
						+eta*ue_iny0k[m];
					pzeta=(1.0-zeta)*ue_ij1[m]
						+zeta*ue_ijnz[m];
					u[ku4+ku3+ku2+m]=pxi+peta+pzeta
						-pxi*peta-peta*pzeta-pzeta*pxi
						+pxi*peta*pzeta;
				}
			}
		}
	});
}

/*
 * ---------------------------------------------------------------------
 * to perform pseudo-time stepping SSOR iterations
 * for five nonlinear pde's.
 * ---------------------------------------------------------------------
 */
void ssor(int niter){
	/*
	 * ---------------------------------------------------------------------
	 * local variables
	 * ---------------------------------------------------------------------
	 */
	int i, j, k, m, n;
	int istep;
	double tmp, tv[ISIZ2*(ISIZ1/2*2+1)*5];
	std::vector<double> delunm(5);
	/*
	 * ---------------------------------------------------------------------
	 * begin pseudo-time stepping iterations
	 * ---------------------------------------------------------------------
	 */
	tmp=1.0/(omega*(2.0-omega));
	/*
	 * ---------------------------------------------------------------------
	 * initialize a,b,c,d to zero (guarantees that page tables have been
	 * formed, if applicable on given architecture, before timestepping).
	 * ---------------------------------------------------------------------
	 */
	std::fill(policy, a.begin(), a.end(), 0.0);
	std::fill(policy, b.begin(), b.end(), 0.0);
	std::fill(policy, c.begin(), c.end(), 0.0);
	std::fill(policy, d.begin(), d.end(), 0.0);
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	/*
	 * ---------------------------------------------------------------------
	 * compute the steady-state residuals
	 * ---------------------------------------------------------------------
	 */
	rhs();
	/*
	 * ---------------------------------------------------------------------
	 * compute the L2 norms of newton iteration residuals
	 * ---------------------------------------------------------------------
	 */
	l2norm(	nx0,
			ny0,
			nz0,
			ist,
			iend,
			jst,
			jend,
			rsd,
			rsdnm);
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	timer_start(1);
	/*
	 * ---------------------------------------------------------------------
	 * the timestep loop
	 * ---------------------------------------------------------------------
	 */
	for(int istep=1; istep<=niter; istep++){
		if((istep%20)==0||istep==itmax||istep==1){
			if(niter>1){
				std::cout
				<< " Time step "
				<< std::setw(4) << istep
				<< std::endl;
			}	
		}
		/*
		 * ---------------------------------------------------------------------
		 * perform SSOR iteration
		 * ---------------------------------------------------------------------
		 */
		if(timeron){timer_start(T_RHS);}
		std::transform(policy, rsd.begin(), rsd.end(), rsd.begin(), [&] (double rsd) {
			return rsd*dt;
		});
		if(timeron){timer_stop(T_RHS);}

		std::for_each_n(policy, iter.front(), num_workers, [&] (int worker_id) {
			for(int k=1; k<nz-1; k++){
				/*
				* ---------------------------------------------------------------------
				* form the lower triangular part of the jacobian matrix
				* ---------------------------------------------------------------------
				*/
				if(timeron){timer_start(T_JACLD);}
				jacld(k, worker_id);
				if(timeron){timer_stop(T_JACLD);}
				/*
				* ---------------------------------------------------------------------
				* perform the lower triangular solution
				* ---------------------------------------------------------------------
				*/
				if(timeron){timer_start(T_BLTS);}
				blts(	nx,
						ny,
						nz,
						k,
						omega,
						rsd,
						a,
						b,
						c,
						d,
						ist,
						iend,
						jst,
						jend,
						nx0,
						ny0, 
						worker_id);
				if(timeron){timer_stop(T_BLTS);}
			}
		});

		std::for_each_n(policy, iter.front(), num_workers, [&] (int worker_id) {

			for(int k=nz-2; k>0; k--){
				/*
				* ---------------------------------------------------------------------
				* form the strictly upper triangular part of the jacobian matrix
				* ---------------------------------------------------------------------
				*/
				if(timeron){timer_start(T_JACU);}
				jacu(k, worker_id);
				if(timeron){timer_stop(T_JACU);}
				/*
				* ---------------------------------------------------------------------
				* perform the upper triangular solution
				* ---------------------------------------------------------------------
				*/
				if(timeron){timer_start(T_BUTS);}
				buts(	nx,
						ny,
						nz,
						k,
						omega,
						rsd,
						d,
						a,
						b,
						c,
						ist,
						iend,
						jst,
						jend,
						nx0,
						ny0,
						worker_id);
				if(timeron){timer_stop(T_BUTS);}
			}
		});
		if(timeron){timer_start(T_ADD);}
		/*
		 * ---------------------------------------------------------------------
		 * update the variables
		 * ---------------------------------------------------------------------
		 */
		std::transform(policy, u.begin(), u.end(), rsd.begin(), u.begin(), [&] (double u, double rsd) {
			return u+tmp*rsd;
		});
		if(timeron){timer_stop(T_ADD);}
		/*
		 * ---------------------------------------------------------------------
		 * compute the max-norms of newton iteration corrections
		 * ---------------------------------------------------------------------
		 */
		if((istep%inorm)==0){
			if(timeron){timer_start(T_L2NORM);}
			l2norm(	nx0,
					ny0,
					nz0,
					ist,
					iend,
					jst,
					jend,
					rsd,
					delunm);
			if(timeron){timer_stop(T_L2NORM);}
		}
		/*
		 * ---------------------------------------------------------------------
		 * compute the steady-state residuals
		 * ---------------------------------------------------------------------
		 */
		rhs();
		/*
		 * ---------------------------------------------------------------------
		 * compute the max-norms of newton iteration residuals
		 * ---------------------------------------------------------------------
		 */
		if(((istep%inorm)==0)||( istep == itmax)){
			if(timeron){timer_start(T_L2NORM);}
			l2norm(	nx0,
					ny0,
					nz0,
					ist,
					iend,
					jst,
					jend,
					rsd,
					rsdnm);
			if(timeron){timer_stop(T_L2NORM);}
		}
		/*
		 * ---------------------------------------------------------------------
		 * check the newton-iteration residuals against the tolerance levels
		 * ---------------------------------------------------------------------
		 */
		if((rsdnm[0]<tolrsd[0])&&
				(rsdnm[1]<tolrsd[1])&&
				(rsdnm[2]<tolrsd[2])&&
				(rsdnm[3]<tolrsd[3])&&
				(rsdnm[4]<tolrsd[4])){
			std::cout
        << std::endl
        << " convergence was achieved after "
        << std::setw(4) << istep
        << " pseudo-time steps"
        << std::endl;
			break;
		}
	}
	timer_stop(1);
	maxtime=timer_read(1);
}

/*
 * ---------------------------------------------------------------------
 * verification routine                         
 * ---------------------------------------------------------------------
 */
void verify(std::vector<double> &xcr,
		std::vector<double> &xce,
		double xci,
		char* class_npb,
		boolean* verified){
	double xciref, xcidif;
	std::vector<double> xcrref(5), xceref(5);
  std::vector<double> xcrdif(5), xcedif(5);
	double epsilon, dtref=0.0;
	int m;
	/*
	 * ---------------------------------------------------------------------
	 * tolerance level
	 * ---------------------------------------------------------------------
	 */
	epsilon=1.0e-08;
	*class_npb='U';
	*verified=TRUE;
	std::fill_n(xcrref.begin(), 5, 1.0);
	std::fill_n(xceref.begin(), 5, 1.0);
	xciref=1.0;
	if((nx0==12)&&(ny0==12)&&(nz0==12)&&(itmax==50)){
		*class_npb='S';
		dtref=5.0e-1;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (12X12X12) grid,
		 * after 50 time steps, with DT = 5.0d-01
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=1.6196343210976702e-02;
		xcrref[1]=2.1976745164821318e-03;
		xcrref[2]=1.5179927653399185e-03;
		xcrref[3]=1.5029584435994323e-03;
		xcrref[4]=3.4264073155896461e-02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (12X12X12) grid,
		 * after 50 time steps, with DT = 5.0d-01
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=6.4223319957960924e-04;
		xceref[1]=8.4144342047347926e-05;
		xceref[2]=5.8588269616485186e-05;
		xceref[3]=5.8474222595157350e-05;
		xceref[4]=1.3103347914111294e-03;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (12X12X12) grid,
		 * after 50 time steps, with DT = 5.0d-01
		 * ---------------------------------------------------------------------
		 */
		xciref=7.8418928865937083e+00;
	}else if((nx0==33)&&(ny0==33)&&(nz0==33)&&(itmax==300)){
		*class_npb='W'; /* SPEC95fp size */
		dtref=1.5e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (33x33x33) grid,
		 * after 300 time steps, with DT = 1.5d-3
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.1236511638192e+02;
		xcrref[1]=0.1317228477799e+01;
		xcrref[2]=0.2550120713095e+01;
		xcrref[3]=0.2326187750252e+01;
		xcrref[4]=0.2826799444189e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (33X33X33) grid,
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.4867877144216e+00;
		xceref[1]=0.5064652880982e-01;
		xceref[2]=0.9281818101960e-01;
		xceref[3]=0.8570126542733e-01;
		xceref[4]=0.1084277417792e+01;
		/*
		 * ---------------------------------------------------------------------
		 * rReference value of surface integral, for the (33X33X33) grid,
		 * after 300 time steps, with DT = 1.5d-3
		 * ---------------------------------------------------------------------
		 */
		xciref=0.1161399311023e+02;
	}else if((nx0==64)&&(ny0==64)&&(nz0==64)&&(itmax==250)){
		*class_npb='A';
		dtref=2.0e+0;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (64X64X64) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=7.7902107606689367e+02;
		xcrref[1]=6.3402765259692870e+01;
		xcrref[2]=1.9499249727292479e+02;
		xcrref[3]=1.7845301160418537e+02;
		xcrref[4]=1.8384760349464247e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (64X64X64) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=2.9964085685471943e+01;
		xceref[1]=2.8194576365003349e+00;
		xceref[2]=7.3473412698774742e+00;
		xceref[3]=6.7139225687777051e+00;
		xceref[4]=7.0715315688392578e+01;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (64X64X64) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=2.6030925604886277e+01;
	}else if((nx0==102)&&(ny0==102)&&(nz0==102)&&(itmax==250)){
		*class_npb='B';
		dtref=2.0e+0;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (102X102X102) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=3.5532672969982736e+03;
		xcrref[1]=2.6214750795310692e+02;
		xcrref[2]=8.8333721850952190e+02;
		xcrref[3]=7.7812774739425265e+02;
		xcrref[4]=7.3087969592545314e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (102X102X102) 
		 * grid, after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=1.1401176380212709e+02;
		xceref[1]=8.1098963655421574e+00;
		xceref[2]=2.8480597317698308e+01;
		xceref[3]=2.5905394567832939e+01;
		xceref[4]=2.6054907504857413e+02;
		/*
		   c---------------------------------------------------------------------
		 * reference value of surface integral, for the (102X102X102) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=4.7887162703308227e+01;
	}else if((nx0==162)&&(ny0==162)&&(nz0==162)&&(itmax==250)){
		*class_npb='C';
		dtref=2.0e+0;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (162X162X162) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=1.03766980323537846e+04;
		xcrref[1]=8.92212458801008552e+02;
		xcrref[2]=2.56238814582660871e+03;
		xcrref[3]=2.19194343857831427e+03;
		xcrref[4]=1.78078057261061185e+04;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (162X162X162) 
		 * grid, after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=2.15986399716949279e+02;
		xceref[1]=1.55789559239863600e+01;
		xceref[2]=5.41318863077207766e+01;
		xceref[3]=4.82262643154045421e+01;
		xceref[4]=4.55902910043250358e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (162X162X162) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=6.66404553572181300e+01;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (162X162X162) grid,
		 * after 250 time steps, with DT = 2.0d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=6.66404553572181300e+01;
	}else if((nx0==408)&&(ny0==408)&&(nz0==408)&&(itmax== 300)){
		*class_npb='D';
		dtref=1.0e+0;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (408X408X408) grid,
		 * after 300 time steps, with DT = 1.0d+00
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.4868417937025e+05;
		xcrref[1]=0.4696371050071e+04;
		xcrref[2]=0.1218114549776e+05;
		xcrref[3]=0.1033801493461e+05;
		xcrref[4]=0.7142398413817e+05;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (408X408X408) 
		 * grid, after 300 time steps, with DT = 1.0d+00
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.3752393004482e+03;
		xceref[1]=0.3084128893659e+02;
		xceref[2]=0.9434276905469e+02;
		xceref[3]=0.8230686681928e+02;
		xceref[4]=0.7002620636210e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (408X408X408) grid,
		 * after 300 time steps, with DT = 1.0d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=0.8334101392503e+02;
	}else if((nx0==1020)&&(ny0==1020)&&(nz0==1020)&&(itmax==300)){
		*class_npb='E';
		dtref=0.5e+0;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual, for the (1020X1020X1020) grid,
		 * after 300 time steps, with DT = 0.5d+00
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.2099641687874e+06;
		xcrref[1]=0.2130403143165e+05;
		xcrref[2]=0.5319228789371e+05;
		xcrref[3]=0.4509761639833e+05;
		xcrref[4]=0.2932360006590e+06;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error, for the (1020X1020X1020) 
		 * grid, after 300 time steps, with DT = 0.5d+00
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.4800572578333e+03;
		xceref[1]=0.4221993400184e+02;
		xceref[2]=0.1210851906824e+03;
		xceref[3]=0.1047888986770e+03;
		xceref[4]=0.8363028257389e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference value of surface integral, for the (1020X1020X1020) grid,
		 * after 300 time steps, with DT = 0.5d+00
		 * ---------------------------------------------------------------------
		 */
		xciref=0.9512163272273e+02;
	}else{
		*verified=FALSE;
	}
	/*
	 * ---------------------------------------------------------------------
	 * verification test for residuals if gridsize is one of 
	 * the defined grid sizes above (class .ne. 'U')
	 * ---------------------------------------------------------------------
	 * compute the difference of solution values and the known reference values.
	 * ---------------------------------------------------------------------
	 */
	std::transform(xcr.begin(), xcr.end(), xcrref.begin(), xcrdif.begin(), [&] (double x, double y) {return fabs((x-y)/y);});
	std::transform(xce.begin(), xce.end(), xceref.begin(), xcedif.begin(), [&] (double x, double y) {return fabs((x-y)/y);});
	xcidif=fabs((xci-xciref)/xciref);
	/*
	 * ---------------------------------------------------------------------
	 * output the comparison of computed results to known cases.
	 * ---------------------------------------------------------------------
	 */
	if(*class_npb!='U'){
    std::cout
      << " Verification being performed for class_npb "
      << *class_npb
      << std::endl;
    
    std::cout
      << " accuracy setting for epsilon = "
      << std::setw(20) << std::setprecision(13) << std::scientific
      << epsilon
      << std::endl;
		*verified=(fabs(dt-dtref)<=epsilon);
		if(!(*verified)){  
			*class_npb='U';
			std::cout
        << " DT does not match the reference value of "
        << std::setw(15) << std::setprecision(8) << std::scientific
        << dtref
        << std::endl;
		} 
	}else{
		std::cout
      << " Unknown class_npb"
      << std::endl;
	}
	if(*class_npb!='U'){
		std::cout
      << " Comparison of RMS-norms of residual"
      << std::endl;
	}else{
    std::cout
      << " RMS-norms of residual"
      << std::endl;
	}
	for(int m=0; m<5; m++){
		if(*class_npb=='U'){
			std::cout
        << "          "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcr[m]
        << std::endl;
		}else if(xcrdif[m]<=epsilon){
      std::cout
        << "          "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcr[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcrref[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcrdif[m]
        << std::endl;
		}else {
			*verified=FALSE;
      std::cout
        << " FAILURE: "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcr[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcrref[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcrdif[m]
        << std::endl;
		}
	}
	if(*class_npb!='U'){
		std::cout
      << " Comparison of RMS-norms of solution error"
      << std::endl;
	}else{
    std::cout
      << " RMS-norms of solution error"
      << std::endl;
	}
	for(int m=0; m<5; m++){
		if(*class_npb=='U'){
			std::cout
        << "          "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xce[m]
        << std::endl;
		}else if(xcedif[m]<=epsilon){
			std::cout
        << "          "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xce[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xceref[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcedif[m]
        << std::endl;
		}else{
			*verified = FALSE;
      std::cout
        << " FAILURE: "
        << std::setw(2) << m+1
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xce[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xceref[m]
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcedif[m]
        << std::endl;
		}
	}
  if(*class_npb!='U'){
		std::cout
      << " Comparison of surface integral"
      << std::endl;
	}else{
		std::cout
      << " Surface integral"
      << std::endl;
	}
  if(*class_npb=='U'){
    std::cout
        << "              "
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xci
        << std::endl;
  }else if(xcidif<=epsilon){
		std::cout
        << "              "
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xci
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xciref
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcidif
        << std::endl;
  }else{
		*verified=FALSE;
    std::cout
        << " FAILURE:     "
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xci
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xciref
        << std::setw(20) << std::setprecision(13) << std::scientific
        << xcidif
        << std::endl;
  }
	if(*class_npb=='U'){
		std::cout
      << " No reference values provided"
      << std::endl;
    
    std::cout
      << " No verification performed"
      << std::endl;
	}else if(*verified){
    std::cout
      << " Verification Successful"
      << std::endl;
	}else{
    std::cout
      << " Verification failed"
      << std::endl;
	}
}
