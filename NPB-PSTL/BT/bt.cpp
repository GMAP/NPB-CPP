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
	R. Van der Wijngaart
	T. Harris
	M. Yarrow
	H. Jin

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

#define IMAX PROBLEM_SIZE
#define JMAX PROBLEM_SIZE
#define KMAX PROBLEM_SIZE
#define IMAXP (IMAX/2*2)
#define JMAXP (JMAX/2*2)
#define AA 0
#define BB 1
#define CC 2
#define BLOCK_SIZE 5
#define T_TOTAL 1
#define T_RHSX 2
#define T_RHSY 3
#define T_RHSZ 4
#define T_RHS 5
#define T_XSOLVE 6
#define T_YSOLVE 7
#define T_ZSOLVE 8
#define T_RDIS1 9
#define T_RDIS2 10
#define T_ADD 11
#define T_LAST 11

/* global variables */
std::vector<double> u((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> rhs((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> forcing((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> us((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> vs((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> ws((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> qs((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> rho_i((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> square((KMAX)*(JMAXP+1)*(IMAXP+1));
std::vector<double> ue((PROBLEM_SIZE+1)*5);
std::vector<double> buf((PROBLEM_SIZE+1)*5);
std::vector<double> lhs((PROBLEM_SIZE+1)*(3)*(5)*(5));
std::vector<double> cuf(PROBLEM_SIZE+1);
std::vector<double> q(PROBLEM_SIZE+1);
std::vector<double> fjac((PROBLEM_SIZE+1)*5*5);
std::vector<double> njac((PROBLEM_SIZE+1)*5*5);
std::vector<double> ce(13*5);

CountIterator iter(0);
static auto policy = std::execution::par;

static double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, 
	      dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, 
	      dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, 
	      dxmax, dymax, dzmax, xxcon1, xxcon2, 
	      xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
	      dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
	      yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
	      zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, 
	      dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, 
	      dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, 
	      c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1,
	      dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, 
	      c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, 
	      c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16,
	      elapsed_time, tmp1, tmp2, tmp3;
std::vector<int> grid_points(3, 0);
static boolean timeron;

/* function prototypes */
static void add();
static void adi();
static void binvcrhs(double tLhs[], double c[], double r[5]);
static void binvrhs(double lhs[], double r[5]);
static void compute_rhs();
static void error_norm(std::vector<double> &rms);
static void exact_rhs();
static void exact_solution(double xi, double eta, double zeta, std::span<double> dtemp);
static void initialize();
static void lhsinit(double* lhs, int size);
static void matmul_sub(double ablock[], double bblock[], double cblock[]);
static void matvec_sub(double ablock[], double avec[5], double bvec[5]);
static void rhs_norm(double rms[5]);
static void set_constants();
static void verify(int no_time_steps, char* class_npb, boolean* verified);
static void x_solve();
static void y_solve();
static void z_solve();

/* bt */
int main(int argc, char* argv[]){
	int niter;
	double navg, mflops, n3;
	double tmax, t;
  std::vector<double> trecs(T_LAST+1, 0);
	boolean verified;
	char class_npb;
	std::vector<std::string> t_names(T_LAST+1);
	/*
	 * ---------------------------------------------------------------------
	 * root node reads input file (if it exists) else takes
	 * defaults from parameters
	 * ---------------------------------------------------------------------
	 */
	std::fstream fp;
  fp.open("inputbt.data", std::ios::in);
	if(fp.is_open()){
		int avoid_warning;
		std::cout
      << " Reading from input file inputbt.data"
      << std::endl;
		fp >> niter;
		while(fp.get() != '\n');
		fp >> dt;
		while(fp.get() != '\n');
		fp >> grid_points[0] >> grid_points[1] >> grid_points[2];
		fp.close();
	}else{
		std::cout
      << " No input file inputbt.data. Using compiled defaults"
      << std::endl;
		niter=NITER_DEFAULT;
		dt=DT_DEFAULT;
		grid_points[0]=PROBLEM_SIZE;
		grid_points[1]=PROBLEM_SIZE;
		grid_points[2]=PROBLEM_SIZE;
	}
  fp.open("timer.flag", std::ios::in);
	if(fp.is_open()){
		timeron=TRUE;
		t_names[T_TOTAL]="total";
		t_names[T_RHSX]="rhsx";
		t_names[T_RHSY]="rhsy";
		t_names[T_RHSZ]="rhsz";
		t_names[T_RHS]="rhs";
		t_names[T_XSOLVE]="xsolve";
		t_names[T_YSOLVE]="ysolve";
		t_names[T_ZSOLVE]="zsolve";
		t_names[T_RDIS1]="redist1";
		t_names[T_RDIS2]="redist2";
		t_names[T_ADD]="add";
		fp.close();
	}else{
		timeron=FALSE;
	}

	iter = CountIterator(std::max({	T_LAST,
									niter+1,
									grid_points[0], 
									grid_points[1], 
									grid_points[2]}));

	std::cout
    << std::endl
    << std::endl;
  
  std::cout
    << " NAS Parallel Benchmarks 4.1 Serial C++ version - BT Benchmark";
  
  std::cout
    << std::endl
    << std::endl;

  std::cout
    << " Size: "
    << std::setw(4) << grid_points[0] << "x"
    << std::setw(4) << grid_points[1] << "x"
    << std::setw(4) << grid_points[2]
    << std::endl;

  std::cout
    << " Iterations: "
    << std::setw(4) << niter
    << "    dt: "
    << std::setw(10) << std::setprecision(6) << dt
    << std::endl;

  std::cout
    << std::endl;
  
	if((grid_points[0]>IMAX)||(grid_points[1]>JMAX)||(grid_points[2]>KMAX)){
		std::cout
      << " " << grid_points[0] << ","
      << " " << grid_points[1] << ","
      << " " << grid_points[2]
      << std::endl;
		
    std::cout
      << " Problem size too big for compiled array sizes"
      << std::endl;
    
		return 0;
	}

	set_constants();
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	initialize();
	exact_rhs();
	/*
	 * ---------------------------------------------------------------------
	 * do one time step to touch all code, and reinitialize
	 * ---------------------------------------------------------------------
	 */
	adi();
	initialize();
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	timer_start(1);
	for(int step=1; step<=niter; step++){
		if((step%20)==0||step==1){
			std::cout
        << " Time step "
        << std::setw(4) << step
        << std::endl;
		}
		adi();
	}
	timer_stop(1);
	tmax=timer_read(1);
	verify(niter, &class_npb, &verified);
	n3=1.0*grid_points[0]*grid_points[1]*grid_points[2];
	navg=(grid_points[0]+grid_points[1]+grid_points[2])/3.0;
	if(tmax!=0.0){
		mflops=1.0e-6*(double)niter*
			(3478.8*n3-17655.7*(navg*navg)+28023.7*navg)
			/tmax;
	}else{
		mflops=0.0;
	}
	c_print_results("BT",
			class_npb,
			grid_points[0],
			grid_points[1],
			grid_points[2],
			niter,
			tmax,
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
			"(none)");
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
          << "  "
          << std::setw(8) << "rest-rhs"
          << ":"
          << std::setw(9) << std::setprecision(3) << t
          << "  ("
          << std::setw(6) << std::setprecision(2) << t*100./tmax
          << "%)"
          << std::endl;
			}else if(i==T_ZSOLVE){
				t=trecs[T_ZSOLVE]-trecs[T_RDIS1]-trecs[T_RDIS2];
        std::cout
          << "  "
          << std::setw(8) << "sub-zsol"
          << ":"
          << std::setw(9) << std::setprecision(3) << t
          << "  ("
          << std::setw(6) << std::setprecision(2) << t*100./tmax
          << "%)"
          << std::endl;
			}else if(i==T_RDIS2){
				t=trecs[T_RDIS1]+trecs[T_RDIS2];
        std::cout
          << "  "
          << std::setw(8) << "redist"
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
 * addition of update to the vector u
 * ---------------------------------------------------------------------
 */
void add(){
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	if(timeron){timer_start(T_ADD);}
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			int ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				int ku2 = i*5;
				for(int m=0; m<5; m++){
					u[ku4+ku3+ku2+m]=u[ku4+ku3+ku2+m]+rhs[ku4+ku3+ku2+m];
				}
			}
		}
	});
	if(timeron){timer_stop(T_ADD);}
}

void adi(){
	compute_rhs();
	x_solve();
	y_solve();
	z_solve();
	add();
}

void binvcrhs(double lhs[], double c[], double r[5]){
	double pivot, coeff;
	pivot=1.00/lhs[0*5+0];
	lhs[1*5+0]=lhs[1*5+0]*pivot;
	lhs[2*5+0]=lhs[2*5+0]*pivot;
	lhs[3*5+0]=lhs[3*5+0]*pivot;
	lhs[4*5+0]=lhs[4*5+0]*pivot;
	c[0*5+0]=c[0*5+0]*pivot;
	c[1*5+0]=c[1*5+0]*pivot;
	c[2*5+0]=c[2*5+0]*pivot;
	c[3*5+0]=c[3*5+0]*pivot;
	c[4*5+0]=c[4*5+0]*pivot;
	r[0]=r[0]*pivot;
	/* */
	coeff=lhs[0*5+1];
	lhs[1*5+1]=lhs[1*5+1]-coeff*lhs[1*5+0];
	lhs[2*5+1]=lhs[2*5+1]-coeff*lhs[2*5+0];
	lhs[3*5+1]=lhs[3*5+1]-coeff*lhs[3*5+0];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+0];
	c[0*5+1]=c[0*5+1]-coeff*c[0*5+0];
	c[1*5+1]=c[1*5+1]-coeff*c[1*5+0];
	c[2*5+1]=c[2*5+1]-coeff*c[2*5+0];
	c[3*5+1]=c[3*5+1]-coeff*c[3*5+0];
	c[4*5+1]=c[4*5+1]-coeff*c[4*5+0];
	r[1]=r[1]-coeff*r[0];
	/* */
	coeff=lhs[0*5+2];
	lhs[1*5+2]=lhs[1*5+2]-coeff*lhs[1*5+0];
	lhs[2*5+2]=lhs[2*5+2]-coeff*lhs[2*5+0];
	lhs[3*5+2]=lhs[3*5+2]-coeff*lhs[3*5+0];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+0];
	c[0*5+2]=c[0*5+2]-coeff*c[0*5+0];
	c[1*5+2]=c[1*5+2]-coeff*c[1*5+0];
	c[2*5+2]=c[2*5+2]-coeff*c[2*5+0];
	c[3*5+2]=c[3*5+2]-coeff*c[3*5+0];
	c[4*5+2]=c[4*5+2]-coeff*c[4*5+0];
	r[2]=r[2]-coeff*r[0];
	/* */
	coeff=lhs[0*5+3];
	lhs[1*5+3]=lhs[1*5+3]-coeff*lhs[1*5+0];
	lhs[2*5+3]=lhs[2*5+3]-coeff*lhs[2*5+0];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+0];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+0];
	c[0*5+3]=c[0*5+3]-coeff*c[0*5+0];
	c[1*5+3]=c[1*5+3]-coeff*c[1*5+0];
	c[2*5+3]=c[2*5+3]-coeff*c[2*5+0];
	c[3*5+3]=c[3*5+3]-coeff*c[3*5+0];
	c[4*5+3]=c[4*5+3]-coeff*c[4*5+0];
	r[3]=r[3]-coeff*r[0];
	/* */
	coeff=lhs[0*5+4];
	lhs[1*5+4]=lhs[1*5+4]-coeff*lhs[1*5+0];
	lhs[2*5+4]=lhs[2*5+4]-coeff*lhs[2*5+0];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+0];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+0];
	c[0*5+4]=c[0*5+4]-coeff*c[0*5+0];
	c[1*5+4]=c[1*5+4]-coeff*c[1*5+0];
	c[2*5+4]=c[2*5+4]-coeff*c[2*5+0];
	c[3*5+4]=c[3*5+4]-coeff*c[3*5+0];
	c[4*5+4]=c[4*5+4]-coeff*c[4*5+0];
	r[4]=r[4]-coeff*r[0];
	/* */
	pivot=1.00/lhs[1*5+1];
	lhs[2*5+1]=lhs[2*5+1]*pivot;
	lhs[3*5+1]=lhs[3*5+1]*pivot;
	lhs[4*5+1]=lhs[4*5+1]*pivot;
	c[0*5+1]=c[0*5+1]*pivot;
	c[1*5+1]=c[1*5+1]*pivot;
	c[2*5+1]=c[2*5+1]*pivot;
	c[3*5+1]=c[3*5+1]*pivot;
	c[4*5+1]=c[4*5+1]*pivot;
	r[1]=r[1]*pivot;
	/* */
	coeff=lhs[1*5+0];
	lhs[2*5+0]=lhs[2*5+0]-coeff*lhs[2*5+1];
	lhs[3*5+0]=lhs[3*5+0]-coeff*lhs[3*5+1];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+1];
	c[0*5+0]=c[0*5+0]-coeff*c[0*5+1];
	c[1*5+0]=c[1*5+0]-coeff*c[1*5+1];
	c[2*5+0]=c[2*5+0]-coeff*c[2*5+1];
	c[3*5+0]=c[3*5+0]-coeff*c[3*5+1];
	c[4*5+0]=c[4*5+0]-coeff*c[4*5+1];
	r[0]=r[0]-coeff*r[1];
	/* */
	coeff = lhs[1*5+2];
	lhs[2*5+2]=lhs[2*5+2]-coeff*lhs[2*5+1];
	lhs[3*5+2]=lhs[3*5+2]-coeff*lhs[3*5+1];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+1];
	c[0*5+2]=c[0*5+2]-coeff*c[0*5+1];
	c[1*5+2]=c[1*5+2]-coeff*c[1*5+1];
	c[2*5+2]=c[2*5+2]-coeff*c[2*5+1];
	c[3*5+2]=c[3*5+2]-coeff*c[3*5+1];
	c[4*5+2]=c[4*5+2]-coeff*c[4*5+1];
	r[2]=r[2]-coeff*r[1];
	/* */
	coeff=lhs[1*5+3];
	lhs[2*5+3]=lhs[2*5+3]-coeff*lhs[2*5+1];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+1];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+1];
	c[0*5+3]=c[0*5+3]-coeff*c[0*5+1];
	c[1*5+3]=c[1*5+3]-coeff*c[1*5+1];
	c[2*5+3]=c[2*5+3]-coeff*c[2*5+1];
	c[3*5+3]=c[3*5+3]-coeff*c[3*5+1];
	c[4*5+3]=c[4*5+3]-coeff*c[4*5+1];
	r[3]=r[3]-coeff*r[1];
	/* */
	coeff=lhs[1*5+4];
	lhs[2*5+4]=lhs[2*5+4]-coeff*lhs[2*5+1];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+1];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+1];
	c[0*5+4]=c[0*5+4]-coeff*c[0*5+1];
	c[1*5+4]=c[1*5+4]-coeff*c[1*5+1];
	c[2*5+4]=c[2*5+4]-coeff*c[2*5+1];
	c[3*5+4]=c[3*5+4]-coeff*c[3*5+1];
	c[4*5+4]=c[4*5+4]-coeff*c[4*5+1];
	r[4]=r[4]-coeff*r[1];
	/* */
	pivot = 1.00/lhs[2*5+2];
	lhs[3*5+2]=lhs[3*5+2]*pivot;
	lhs[4*5+2]=lhs[4*5+2]*pivot;
	c[0*5+2]=c[0*5+2]*pivot;
	c[1*5+2]=c[1*5+2]*pivot;
	c[2*5+2]=c[2*5+2]*pivot;
	c[3*5+2]=c[3*5+2]*pivot;
	c[4*5+2]=c[4*5+2]*pivot;
	r[2]=r[2]*pivot;
	/* */
	coeff=lhs[2*5+0];
	lhs[3*5+0]=lhs[3*5+0]-coeff*lhs[3*5+2];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+2];
	c[0*5+0]=c[0*5+0]-coeff*c[0*5+2];
	c[1*5+0]=c[1*5+0]-coeff*c[1*5+2];
	c[2*5+0]=c[2*5+0]-coeff*c[2*5+2];
	c[3*5+0]=c[3*5+0]-coeff*c[3*5+2];
	c[4*5+0]=c[4*5+0]-coeff*c[4*5+2];
	r[0]=r[0]-coeff*r[2];
	/* */
	coeff=lhs[2*5+1];
	lhs[3*5+1]=lhs[3*5+1]-coeff*lhs[3*5+2];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+2];
	c[0*5+1]=c[0*5+1]-coeff*c[0*5+2];
	c[1*5+1]=c[1*5+1]-coeff*c[1*5+2];
	c[2*5+1]=c[2*5+1]-coeff*c[2*5+2];
	c[3*5+1]=c[3*5+1]-coeff*c[3*5+2];
	c[4*5+1]=c[4*5+1]-coeff*c[4*5+2];
	r[1]=r[1]-coeff*r[2];
	/* */
	coeff=lhs[2*5+3];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+2];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+2];
	c[0*5+3]=c[0*5+3]-coeff*c[0*5+2];
	c[1*5+3]=c[1*5+3]-coeff*c[1*5+2];
	c[2*5+3]=c[2*5+3]-coeff*c[2*5+2];
	c[3*5+3]=c[3*5+3]-coeff*c[3*5+2];
	c[4*5+3]=c[4*5+3]-coeff*c[4*5+2];
	r[3]=r[3]-coeff*r[2];
	/* */
	coeff=lhs[2*5+4];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+2];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+2];
	c[0*5+4]=c[0*5+4]-coeff*c[0*5+2];
	c[1*5+4]=c[1*5+4]-coeff*c[1*5+2];
	c[2*5+4]=c[2*5+4]-coeff*c[2*5+2];
	c[3*5+4]=c[3*5+4]-coeff*c[3*5+2];
	c[4*5+4]=c[4*5+4]-coeff*c[4*5+2];
	r[4]=r[4]-coeff*r[2];
	/* */
	pivot=1.00/lhs[3*5+3];
	lhs[4*5+3]=lhs[4*5+3]*pivot;
	c[0*5+3]=c[0*5+3]*pivot;
	c[1*5+3]=c[1*5+3]*pivot;
	c[2*5+3]=c[2*5+3]*pivot;
	c[3*5+3]=c[3*5+3]*pivot;
	c[4*5+3]=c[4*5+3]*pivot;
	r[3]=r[3] *pivot;
	/* */
	coeff=lhs[3*5+0];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+3];
	c[0*5+0]=c[0*5+0]-coeff*c[0*5+3];
	c[1*5+0]=c[1*5+0]-coeff*c[1*5+3];
	c[2*5+0]=c[2*5+0]-coeff*c[2*5+3];
	c[3*5+0]=c[3*5+0]-coeff*c[3*5+3];
	c[4*5+0]=c[4*5+0]-coeff*c[4*5+3];
	r[0]=r[0]-coeff*r[3];
	/* */
	coeff=lhs[3*5+1];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+3];
	c[0*5+1]=c[0*5+1]-coeff*c[0*5+3];
	c[1*5+1]=c[1*5+1]-coeff*c[1*5+3];
	c[2*5+1]=c[2*5+1]-coeff*c[2*5+3];
	c[3*5+1]=c[3*5+1]-coeff*c[3*5+3];
	c[4*5+1]=c[4*5+1]-coeff*c[4*5+3];
	r[1]=r[1]-coeff*r[3];
	/* */
	coeff=lhs[3*5+2];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+3];
	c[0*5+2]=c[0*5+2]-coeff*c[0*5+3];
	c[1*5+2]=c[1*5+2]-coeff*c[1*5+3];
	c[2*5+2]=c[2*5+2]-coeff*c[2*5+3];
	c[3*5+2]=c[3*5+2]-coeff*c[3*5+3];
	c[4*5+2]=c[4*5+2]-coeff*c[4*5+3];
	r[2]=r[2]-coeff*r[3];
	/* */
	coeff=lhs[3*5+4];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+3];
	c[0*5+4]=c[0*5+4]-coeff*c[0*5+3];
	c[1*5+4]=c[1*5+4]-coeff*c[1*5+3];
	c[2*5+4]=c[2*5+4]-coeff*c[2*5+3];
	c[3*5+4]=c[3*5+4]-coeff*c[3*5+3];
	c[4*5+4]=c[4*5+4]-coeff*c[4*5+3];
	r[4]=r[4]-coeff*r[3];
	/* */
	pivot=1.00/lhs[4*5+4];
	c[0*5+4]=c[0*5+4]*pivot;
	c[1*5+4]=c[1*5+4]*pivot;
	c[2*5+4]=c[2*5+4]*pivot;
	c[3*5+4]=c[3*5+4]*pivot;
	c[4*5+4]=c[4*5+4]*pivot;
	r[4]=r[4]*pivot;
	/* */
	coeff=lhs[4*5+0];
	c[0*5+0]=c[0*5+0]-coeff*c[0*5+4];
	c[1*5+0]=c[1*5+0]-coeff*c[1*5+4];
	c[2*5+0]=c[2*5+0]-coeff*c[2*5+4];
	c[3*5+0]=c[3*5+0]-coeff*c[3*5+4];
	c[4*5+0]=c[4*5+0]-coeff*c[4*5+4];
	r[0]=r[0]-coeff*r[4];
	/* */
	coeff=lhs[4*5+1];
	c[0*5+1]=c[0*5+1]-coeff*c[0*5+4];
	c[1*5+1]=c[1*5+1]-coeff*c[1*5+4];
	c[2*5+1]=c[2*5+1]-coeff*c[2*5+4];
	c[3*5+1]=c[3*5+1]-coeff*c[3*5+4];
	c[4*5+1]=c[4*5+1]-coeff*c[4*5+4];
	r[1]=r[1]-coeff*r[4];
	/* */
	coeff=lhs[4*5+2];
	c[0*5+2]=c[0*5+2]-coeff*c[0*5+4];
	c[1*5+2]=c[1*5+2]-coeff*c[1*5+4];
	c[2*5+2]=c[2*5+2]-coeff*c[2*5+4];
	c[3*5+2]=c[3*5+2]-coeff*c[3*5+4];
	c[4*5+2]=c[4*5+2]-coeff*c[4*5+4];
	r[2]=r[2]-coeff*r[4];
	/* */
	coeff=lhs[4*5+3];
	c[0*5+3]=c[0*5+3]-coeff*c[0*5+4];
	c[1*5+3]=c[1*5+3]-coeff*c[1*5+4];
	c[2*5+3]=c[2*5+3]-coeff*c[2*5+4];
	c[3*5+3]=c[3*5+3]-coeff*c[3*5+4];
	c[4*5+3]=c[4*5+3]-coeff*c[4*5+4];
	r[3]=r[3]-coeff*r[4];
}

void binvrhs(double lhs[], double r[5]){
	double pivot, coeff;
	pivot=1.00/lhs[0*5+0];
	lhs[1*5+0]=lhs[1*5+0]*pivot;
	lhs[2*5+0]=lhs[2*5+0]*pivot;
	lhs[3*5+0]=lhs[3*5+0]*pivot;
	lhs[4*5+0]=lhs[4*5+0]*pivot;
	r[0]=r[0]*pivot;
	/* */
	coeff=lhs[0*5+1];
	lhs[1*5+1]=lhs[1*5+1]-coeff*lhs[1*5+0];
	lhs[2*5+1]=lhs[2*5+1]-coeff*lhs[2*5+0];
	lhs[3*5+1]=lhs[3*5+1]-coeff*lhs[3*5+0];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+0];
	r[1]=r[1]-coeff*r[0];
	/* */
	coeff=lhs[0*5+2];
	lhs[1*5+2]=lhs[1*5+2]-coeff*lhs[1*5+0];
	lhs[2*5+2]=lhs[2*5+2]-coeff*lhs[2*5+0];
	lhs[3*5+2]=lhs[3*5+2]-coeff*lhs[3*5+0];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+0];
	r[2]=r[2]-coeff*r[0];
	/* */
	coeff=lhs[0*5+3];
	lhs[1*5+3]=lhs[1*5+3]-coeff*lhs[1*5+0];
	lhs[2*5+3]=lhs[2*5+3]-coeff*lhs[2*5+0];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+0];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+0];
	r[3]=r[3]-coeff*r[0];
	/* */
	coeff=lhs[0*5+4];
	lhs[1*5+4]=lhs[1*5+4]-coeff*lhs[1*5+0];
	lhs[2*5+4]=lhs[2*5+4]-coeff*lhs[2*5+0];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+0];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+0];
	r[4]=r[4]-coeff*r[0];
	/* */
	pivot=1.00/lhs[1*5+1];
	lhs[2*5+1]=lhs[2*5+1]*pivot;
	lhs[3*5+1]=lhs[3*5+1]*pivot;
	lhs[4*5+1]=lhs[4*5+1]*pivot;
	r[1]=r[1]*pivot;
	/* */
	coeff=lhs[1*5+0];
	lhs[2*5+0]=lhs[2*5+0]-coeff*lhs[2*5+1];
	lhs[3*5+0]=lhs[3*5+0]-coeff*lhs[3*5+1];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+1];
	r[0]=r[0]-coeff*r[1];
	/* */
	coeff=lhs[1*5+2];
	lhs[2*5+2]=lhs[2*5+2]-coeff*lhs[2*5+1];
	lhs[3*5+2]=lhs[3*5+2]-coeff*lhs[3*5+1];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+1];
	r[2]=r[2]-coeff*r[1];
	/* */
	coeff=lhs[1*5+3];
	lhs[2*5+3]=lhs[2*5+3]-coeff*lhs[2*5+1];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+1];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+1];
	r[3]=r[3]-coeff*r[1];
	/* */
	coeff=lhs[1*5+4];
	lhs[2*5+4]=lhs[2*5+4]-coeff*lhs[2*5+1];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+1];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+1];
	r[4]=r[4]-coeff*r[1];
	/* */
	pivot=1.00/lhs[2*5+2];
	lhs[3*5+2]=lhs[3*5+2]*pivot;
	lhs[4*5+2]=lhs[4*5+2]*pivot;
	r[2]=r[2]*pivot;
	/* */
	coeff=lhs[2*5+0];
	lhs[3*5+0]=lhs[3*5+0]-coeff*lhs[3*5+2];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+2];
	r[0]=r[0]-coeff*r[2];
	/* */
	coeff=lhs[2*5+1];
	lhs[3*5+1]=lhs[3*5+1]-coeff*lhs[3*5+2];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+2];
	r[1]=r[1]-coeff*r[2];
	/* */
	coeff=lhs[2*5+3];
	lhs[3*5+3]=lhs[3*5+3]-coeff*lhs[3*5+2];
	lhs[4*5+3]=lhs[4*5+3]-coeff*lhs[4*5+2];
	r[3]=r[3]-coeff*r[2];
	/* */
	coeff=lhs[2*5+4];
	lhs[3*5+4]=lhs[3*5+4]-coeff*lhs[3*5+2];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+2];
	r[4]=r[4]-coeff*r[2];
	/* */
	pivot=1.00/lhs[3*5+3];
	lhs[4*5+3]=lhs[4*5+3]*pivot;
	r[3]=r[3]*pivot;
	/* */
	coeff=lhs[3*5+0];
	lhs[4*5+0]=lhs[4*5+0]-coeff*lhs[4*5+3];
	r[0]=r[0]-coeff*r[3];
	/* */
	coeff=lhs[3*5+1];
	lhs[4*5+1]=lhs[4*5+1]-coeff*lhs[4*5+3];
	r[1]=r[1]-coeff*r[3];
	/* */
	coeff=lhs[3*5+2];
	lhs[4*5+2]=lhs[4*5+2]-coeff*lhs[4*5+3];
	r[2]=r[2]-coeff*r[3];
	/* */
	coeff=lhs[3*5+4];
	lhs[4*5+4]=lhs[4*5+4]-coeff*lhs[4*5+3];
	r[4]=r[4]-coeff*r[3];
	/* */
	pivot=1.00/lhs[4*5+4];
	r[4]=r[4]*pivot;
	/* */
	coeff=lhs[4*5+0];
	r[0]=r[0]-coeff*r[4];
	/* */
	coeff=lhs[4*5+1];
	r[1]=r[1]-coeff*r[4];
	/* */
	coeff=lhs[4*5+2];
	r[2]=r[2]-coeff*r[4];
	/* */
	coeff=lhs[4*5+3];
	r[3]=r[3]-coeff*r[4];
}

void compute_rhs(){
	int i, j, k, m;
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	int ks3, ks2;
	int ks3_const=(JMAXP+1)*(IMAXP+1), ks2_const=(IMAXP+1);
	double rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;

	if(timeron){timer_start(T_RHS);}
	/*
	 * ---------------------------------------------------------------------
	 * compute the reciprocal of density, and the kinetic energy, 
	 * and the speed of sound.
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4 = k*ku4_const;
		int ks3 = k*ks3_const;
		for(int j=0; j<=grid_points[1]-1; j++){
			int ks2 = j*ks2_const;
			int ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-1; i++){
				int ku2 = i*5;
				double rho_inv=1.0/u[ku4+ku3+ku2];
				rho_i[ks3+ks2+i]=rho_inv;
				us[ks3+ks2+i]=u[ku4+ku3+ku2+1]*rho_inv;
				vs[ks3+ks2+i]=u[ku4+ku3+ku2+2]*rho_inv;
				ws[ks3+ks2+i]=u[ku4+ku3+ku2+3]*rho_inv;
				square[ks3+ks2+i]=0.5*(
						u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]+ 
						u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]+
						u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])*rho_inv;
				qs[ks3+ks2+i]=square[ks3+ks2+i]*rho_inv;
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * copy the exact forcing term to the right hand side; because 
	 * this forcing term is known, we can store it on the whole grid
	 * including the boundary                   
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4 = k*ku4_const;
		for(int j=0; j<=grid_points[1]-1; j++){
			int ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-1; i++){
				int ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m];
				}
			}
		}
	});
	if(timeron){timer_start(T_RHSX);}
	/*
	 * ---------------------------------------------------------------------
	 * compute xi-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2;
		double uijk, up1, um1;
		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			ks2 = j*ks2_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				uijk=us[ks3+ks2+i];
				up1=us[ks3+ks2+i+1];
				um1=us[ks3+ks2+i-1];
				rhs[ku4+ku3+ku2]=rhs[ku4+ku3+ku2]+dx1tx1* 
					(u[ku4+ku3+ku2+5]-2.0*u[ku4+ku3+ku2]+ 
					 u[ku4+ku3+ku2-5+0])-
					tx2*(u[ku4+ku3+ku2+5+1]-u[ku4+ku3+ku2-5+1]);
				rhs[ku4+ku3+ku2+1]=rhs[ku4+ku3+ku2+1]+dx2tx1* 
					(u[ku4+ku3+ku2+5+1]-2.0*u[ku4+ku3+ku2+1]+ 
					 u[ku4+ku3+ku2-5+1])+
					xxcon2*con43*(up1-2.0*uijk+um1)-
					tx2*(u[ku4+ku3+ku2+5+1]*up1- 
							u[ku4+ku3+ku2-5+1]*um1+
							(u[ku4+ku3+ku2+5+4]- square[ks3+ks2+i+1]-
							 u[ku4+ku3+ku2-5+4]+ square[ks3+ks2+i-1])*
							c2);
				rhs[ku4+ku3+ku2+2]=rhs[ku4+ku3+ku2+2]+dx3tx1* 
					(u[ku4+ku3+ku2+5+2]-2.0*u[ku4+ku3+ku2+2]+
					 u[ku4+ku3+ku2-5+2])+
					xxcon2*(vs[ks3+ks2+i+1]-2.0*vs[ks3+ks2+i]+
							vs[ks3+ks2+i-1])-
					tx2*(u[ku4+ku3+ku2+5+2]*up1- 
							u[ku4+ku3+ku2-5+2]*um1);
				rhs[ku4+ku3+ku2+3]=rhs[ku4+ku3+ku2+3]+dx4tx1* 
					(u[ku4+ku3+ku2+5+3]-2.0*u[ku4+ku3+ku2+3]+
					 u[ku4+ku3+ku2-5+3])+
					xxcon2*(ws[ks3+ks2+i+1]-2.0*ws[ks3+ks2+i]+
							ws[ks3+ks2+i-1])-
					tx2*(u[ku4+ku3+ku2+5+3]*up1- 
							u[ku4+ku3+ku2-5+3]*um1);
				rhs[ku4+ku3+ku2+4]=rhs[ku4+ku3+ku2+4]+dx5tx1* 
					(u[ku4+ku3+ku2+5+4]-2.0*u[ku4+ku3+ku2+4]+
					 u[ku4+ku3+ku2-5+4])+
					xxcon3*(qs[ks3+ks2+i+1]-2.0*qs[ks3+ks2+i]+
							qs[ks3+ks2+i-1])+
					xxcon4*(up1*up1-2.0*uijk*uijk+ 
							um1*um1)+
					xxcon5*(u[ku4+ku3+ku2+5+4]*rho_i[ks3+ks2+i+1]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4+ku3+ku2-5+4]*rho_i[ks3+ks2+i-1])-
					tx2*((c1*u[ku4+ku3+ku2+5+4]- 
								c2*square[ks3+ks2+i+1])*up1-
							(c1*u[ku4+ku3+ku2-5+4]- 
							 c2*square[ks3+ks2+i-1])*um1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order xi-direction dissipation               
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			i=1;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+5+m]=rhs[ku4+ku3+5+m]-dssp* 
					(5.0*u[ku4+ku3+5+m]-4.0*u[ku4+ku3+10+m]+
					 u[ku4+ku3+15+m]);
			}
			i=2;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+10+m]=rhs[ku4+ku3+10+m]-dssp* 
					(-4.0*u[ku4+ku3+5+m]+6.0*u[ku4+ku3+10+m]-
					 4.0*u[ku4+ku3+15+m]+u[ku4+ku3+20+m]);
			}
		}
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			for(int i=3; i<=grid_points[0]-4; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp * 
						(u[ku4+ku3+ku2-10+m]-4.0*u[ku4+ku3+ku2-5+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3+ku2+5+m]+ 
						 u[ku4+ku3+ku2+10+m]);
				}
			}
		}
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			i=grid_points[0]-3;
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp*
					(u[ku4+ku3+ku2-10+m]-4.0*u[ku4+ku3+ku2-5+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3+ku2+5+m]);
			}
			i=grid_points[0]-2;
			ku2 += 5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp*
					(u[ku4+ku3+ku2-10+m]-4.*u[ku4+ku3+ku2-5+m]+
					 5.*u[ku4+ku3+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSX);}
	if(timeron){timer_start(T_RHSY);}
	/*
	 * ---------------------------------------------------------------------
	 * compute eta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2;
    	double vijk, vp1, vm1;
		ku4 = k*ku4_const;
		ks3 = k*ks3_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			ks2 = j*ks2_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				vijk=vs[ks3+ks2+i];
				vp1=vs[ks3+ks2+ks2_const+i];
				vm1=vs[ks3+ks2-ks2_const+i];
				rhs[ku4+ku3+ku2]=rhs[ku4+ku3+ku2]+dy1ty1* 
					(u[ku4+ku3+ku3_const+ku2]-2.0*u[ku4+ku3+ku2]+ 
					 u[ku4+ku3-ku3_const+ku2])-
					ty2*(u[ku4+ku3+ku3_const+ku2+2]-u[ku4+ku3-ku3_const+ku2+2]);
				rhs[ku4+ku3+ku2+1]=rhs[ku4+ku3+ku2+1]+dy2ty1* 
					(u[ku4+ku3+ku3_const+ku2+1]-2.0*u[ku4+ku3+ku2+1]+ 
					 u[ku4+ku3-ku3_const+ku2+1])+
					yycon2*(us[ks3+ks2+ks2_const+i]-2.0*us[ks3+ks2+i]+ 
							us[ks3+ks2-ks2_const+i])-
					ty2*(u[ku4+ku3+ku3_const+ku2+1]*vp1- 
							u[ku4+ku3-ku3_const+ku2+1]*vm1);
				rhs[ku4+ku3+ku2+2]=rhs[ku4+ku3+ku2+2]+dy3ty1* 
					(u[ku4+ku3+ku3_const+ku2+2]-2.0*u[ku4+ku3+ku2+2]+ 
					 u[ku4+ku3-ku3_const+ku2+2])+
					yycon2*con43*(vp1-2.0*vijk+vm1)-
					ty2*(u[ku4+ku3+ku3_const+ku2+2]*vp1- 
							u[ku4+ku3-ku3_const+ku2+2]*vm1+
							(u[ku4+ku3+ku3_const+ku2+4]-square[ks3+ks2+ks2_const+i]- 
							 u[ku4+ku3-ku3_const+ku2+4]+square[ks3+ks2-ks2_const+i])
							*c2);
				rhs[ku4+ku3+ku2+3]=rhs[ku4+ku3+ku2+3]+dy4ty1* 
					(u[ku4+ku3+ku3_const+ku2+3]-2.0*u[ku4+ku3+ku2+3]+ 
					 u[ku4+ku3-ku3_const+ku2+3])+
					yycon2*(ws[ks3+ks2+ks2_const+i]-2.0*ws[ks3+ks2+i]+ 
							ws[ks3+ks2-ks2_const+i])-
					ty2*(u[ku4+ku3+ku3_const+ku2+3]*vp1- 
							u[ku4+ku3-ku3_const+ku2+3]*vm1);
				rhs[ku4+ku3+ku2+4]=rhs[ku4+ku3+ku2+4]+dy5ty1* 
					(u[ku4+ku3+ku3_const+ku2+4]-2.0*u[ku4+ku3+ku2+4]+ 
					 u[ku4+ku3-ku3_const+ku2+4])+
					yycon3*(qs[ks3+ks2+ks2_const+i]-2.0*qs[ks3+ks2+i]+ 
							qs[ks3+ks2-ks2_const+i])+
					yycon4*(vp1*vp1-2.0*vijk*vijk+ 
							vm1*vm1)+
					yycon5*(u[ku4+ku3+ku3_const+ku2+4]*rho_i[ks3+ks2+ks2_const+i]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4+ku3-ku3_const+ku2+4]*rho_i[ks3+ks2-ks2_const+i])-
					ty2*((c1*u[ku4+ku3+ku3_const+ku2+4]- 
								c2*square[ks3+ks2+ks2_const+i])*vp1-
							(c1*u[ku4+ku3-ku3_const+ku2+4]- 
							 c2*square[ks3+ks2-ks2_const+i])*vm1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order eta-direction dissipation         
		 * ---------------------------------------------------------------------
		 */
		int ku3P, ku3M;
		j=1;
		ku3 = ku3_const;
		ku3P = ku3+ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
					(5.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3P+ku2+m]+
					 u[ku4+ku3P+ku3_const+ku2+m]);
			}
		}
		j=2;
		ku3 = j*ku3_const;
		ku3P += ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
					(-4.0*u[ku4+ku3_const+ku2+m]+6.0*u[ku4+ku3P-ku3_const+ku2+m]-
					 4.0*u[ku4+ku3P+ku2+m]+u[ku4+ku3P+ku3_const+ku2+m]);
			}
		}
		
		for(int j=3; j<=grid_points[1]-4; j++){
			ku3 = j*ku3_const;
			ku3P = ku3+ku3_const;
			ku3M = ku3-ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
						(u[ku4+ku3M-ku3_const+ku2+m]-4.0*u[ku4+ku3M+ku2+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3P+ku2+m]+ 
						 u[ku4+ku3P+ku3_const+ku2+m]);
				}
			}
		}
		j=grid_points[1]-3;
		ku3 = j*ku3_const;
		ku3M = ku3-ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp*
					(u[ku4+ku3M-ku3_const+ku2+m]-4.0*u[ku4+ku3M+ku2+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3+ku3_const+ku2+m]);
			}
		}
		j=grid_points[1]-2;
		ku3P = ku3+ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3P+ku2+m]=rhs[ku4+ku3P+ku2+m]-dssp*
					(u[ku4+ku3M+ku2+m]-4.*u[ku4+ku3+ku2+m]+
					 5.*u[ku4+ku3P+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSY);}
	if(timeron){timer_start(T_RHSZ);}
	/*
	 * ---------------------------------------------------------------------
	 * compute zeta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2;
    	double wijk, wp1, wm1;
		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ks2 = j*ks2_const;
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				wijk=ws[ks3+ks2+i];
				wp1=ws[ks3+ks3_const+ks2+i];
				wm1=ws[ks3-ks3_const+ks2+i];
				rhs[ku4+ku3+ku2+0]=rhs[ku4+ku3+ku2+0]+dz1tz1* 
					(u[ku4+ku4_const+ku3+ku2]-2.0*u[ku4+ku3+ku2]+ 
					 u[ku4-ku4_const+ku3+ku2])-
					tz2*(u[ku4+ku4_const+ku3+ku2+3]-u[ku4-ku4_const+ku3+ku2+3]);
				rhs[ku4+ku3+ku2+1]=rhs[ku4+ku3+ku2+1]+dz2tz1* 
					(u[ku4+ku4_const+ku3+ku2+1]-2.0*u[ku4+ku3+ku2+1]+ 
					 u[ku4-ku4_const+ku3+ku2+1])+
					zzcon2*(us[ks3+ks3_const+ks2+i]-2.0*us[ks3+ks2+i]+ 
							us[ks3-ks3_const+ks2+i])-
					tz2*(u[ku4+ku4_const+ku3+ku2+1]*wp1- 
							u[ku4-ku4_const+ku3+ku2+1]*wm1);
				rhs[ku4+ku3+ku2+2]=rhs[ku4+ku3+ku2+2]+dz3tz1* 
					(u[ku4+ku4_const+ku3+ku2+2]-2.0*u[ku4+ku3+ku2+2]+ 
					 u[ku4-ku4_const+ku3+ku2+2])+
					zzcon2*(vs[ks3+ks3_const+ks2+i]-2.0*vs[ks3+ks2+i]+ 
							vs[ks3-ks3_const+ks2+i])-
					tz2*(u[ku4+ku4_const+ku3+ku2+2]*wp1- 
							u[ku4-ku4_const+ku3+ku2+2]*wm1);
				rhs[ku4+ku3+ku2+3]=rhs[ku4+ku3+ku2+3]+dz4tz1* 
					(u[ku4+ku4_const+ku3+ku2+3]-2.0*u[ku4+ku3+ku2+3]+ 
					 u[ku4-ku4_const+ku3+ku2+3])+
					zzcon2*con43*(wp1-2.0*wijk+wm1)-
					tz2*(u[ku4+ku4_const+ku3+ku2+3]*wp1- 
							u[ku4-ku4_const+ku3+ku2+3]*wm1+
							(u[ku4+ku4_const+ku3+ku2+4]-square[ks3+ks3_const+ks2+i]- 
							 u[ku4-ku4_const+ku3+ku2+4]+square[ks3-ks3_const+ks2+i])
							*c2);
				rhs[ku4+ku3+ku2+4]=rhs[ku4+ku3+ku2+4]+dz5tz1* 
					(u[ku4+ku4_const+ku3+ku2+4]-2.0*u[ku4+ku3+ku2+4]+ 
					 u[ku4-ku4_const+ku3+ku2+4])+
					zzcon3*(qs[ks3+ks3_const+ks2+i]-2.0*qs[ks3+ks2+i]+ 
							qs[ks3-ks3_const+ks2+i])+
					zzcon4*(wp1*wp1-2.0*wijk*wijk+ 
							wm1*wm1)+
					zzcon5*(u[ku4+ku4_const+ku3+ku2+4]*rho_i[ks3+ks3_const+ks2+i]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4-ku4_const+ku3+ku2+4]*rho_i[ks3-ks3_const+ks2+i])-
					tz2*((c1*u[ku4+ku4_const+ku3+ku2+4]- 
								c2*square[ks3+ks3_const+ks2+i])*wp1-
							(c1*u[ku4-ku4_const+ku3+ku2+4]- 
							 c2*square[ks3-ks3_const+ks2+i])*wm1);
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * add fourth order zeta-direction dissipation                
	 * ---------------------------------------------------------------------
	 */
	int ku4M, ku4P;
	k=1;
	ku4 = ku4_const;
	ku4P = ku4+ku4_const;
	std::for_each_n(policy, iter.front()+1, grid_points[1]-2, [&](int j){
		int ku3, ku2;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
					(5.0*u[ku4+ku3+ku2+m]-4.0*u[ku4P+ku3+ku2+m]+
					 u[ku4P+ku4_const+ku3+ku2+m]);
			}
		}
	});
	k=2;
	ku4 += ku4_const;
	ku4P += ku4_const;
	std::for_each_n(policy, iter.front()+1, grid_points[1]-2, [&](int j){
		int ku3, ku2;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
					(-4.0*u[ku4-ku4_const+ku3+ku2+m]+6.0*u[ku4+ku3+ku2+m]-
					 4.0*u[ku4P+ku3+ku2+m]+u[ku4P+ku4_const+ku3+ku2+m]);
			}
		}
	});
	std::for_each_n(policy, iter.front()+3, grid_points[2]-6, [&](int k){
		int ku4, ku3, ku2, ku4M, ku4P;
		ku4 = k*ku4_const;
		ku4M = ku4-ku4_const;
		ku4P = ku4+ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp* 
						(u[ku4M-ku4_const+ku3+ku2+m]-4.0*u[ku4M+ku3+ku2+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4P+ku3+ku2+m]+ 
						 u[ku4P+ku4_const+ku3+ku2+m]);
				}
			}
		}
	});
	k=grid_points[2]-3;
	ku4 = k*ku4_const;
	ku4M = ku4-ku4_const;
	std::for_each_n(policy, iter.front()+1, grid_points[1]-2, [&](int j){
		int ku3, ku2;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]-dssp*
					(u[ku4M-ku4_const+ku3+ku2+m]-4.0*u[ku4M+ku3+ku2+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku4_const+ku3+ku2+m]);
			}
		}
	});
	k=grid_points[2]-2;
	ku4P = ku4+ku4_const;
	std::for_each_n(policy, iter.front()+1, grid_points[1]-2, [&](int j){
		int ku3, ku2;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4P+ku3+ku2+m]=rhs[ku4P+ku3+ku2+m]-dssp*
					(u[ku4M+ku3+ku2+m]-4.*u[ku4+ku3+ku2+m]+
					 5.*u[ku4P+ku3+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSZ);}
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2;
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]=rhs[ku4+ku3+ku2+m]*dt;
				}
			}
		}
	});
	if(timeron){timer_stop(T_RHS);}
}

/*
 * ---------------------------------------------------------------------
 * this function computes the norm of the difference between the
 * computed solution and the exact solution
 * ---------------------------------------------------------------------
 */
void error_norm(std::vector<double> &rms){
	int ku4, ku3, ku2;
	std::vector<double> u_exact(5, 0.0);
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	double xi, eta, zeta, add;

	std::fill(rms.begin(), rms.end(), 0.0);
	for(int k=0; k<=grid_points[2]-1; k++){
		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
			ku3 = j*ku3_const;
			eta=(double)(j)*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
				ku2 = i*5;
				xi=(double)(i)*dnxm1;
				exact_solution(xi, eta, zeta, u_exact);
				for(int m=0; m<5; m++){
					add=u[ku4+ku3+ku2+m]-u_exact[m];
					rms[m]=rms[m]+add*add;
				}
			}
		}
	}
	for(int m=0; m<5; m++){
		for(int d=0; d<3; d++){
			rms[m]/=(double)(grid_points[d]-2);
		}
		rms[m]=sqrt(rms[m]);
	}
}

/*
 * ---------------------------------------------------------------------
 * compute the right hand side based on exact solution
 * ---------------------------------------------------------------------
 */
void exact_rhs(){
	std::vector<double> dtemp(5, 0.0);
	double xi, eta, zeta, dtpp;
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	int m_calc, m_const=(PROBLEM_SIZE+1);
	int m_constT2=2*m_const, m_constT3=3*m_const, m_constT4=4*m_const;
	int i, j, k, ip1, im1, jp1, jm1, km1, kp1;
	/*
	 * ---------------------------------------------------------------------
	 * initialize                                  
	 * ---------------------------------------------------------------------
	 */
	std::fill(policy, forcing.begin(), forcing.end(), 0.0);
	/*
	 * ---------------------------------------------------------------------
	 * xi-direction flux differences                      
	 * ---------------------------------------------------------------------
	 */
	for(int k=1; k<=grid_points[2]-2; k++){
		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			eta=(double)(j)*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
				xi=(double)(i)*dnxm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
					m_calc = m*(PROBLEM_SIZE+1);
					ue[m_calc+i]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
					m_calc = m*(PROBLEM_SIZE+1);
					buf[m_calc+i]=dtpp*dtemp[m];
				}
				cuf[i]=buf[m_const+i]*buf[m_const+i];
				buf[i]=cuf[i]+buf[m_constT2+i]*buf[m_constT2+i]+buf[m_constT3+i]*buf[m_constT3+i];
				q[i]=0.5*(buf[m_const+i]*ue[m_const+i]+buf[m_constT2+i]*ue[m_constT2+i]+
						buf[m_constT3+i]*ue[m_constT3+i]);
			}
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				im1=i-1;
				ip1=i+1;
				forcing[ku4+ku3+ku2+0]=forcing[ku4+ku3+ku2+0]-
					tx2*(ue[m_const+ip1]-ue[m_const+im1])+
					dx1tx1*(ue[ip1]-2.0*ue[i]+ue[im1]);
				forcing[ku4+ku3+ku2+1]=forcing[ku4+ku3+ku2+1]-tx2*(
						(ue[m_const+ip1]*buf[m_const+ip1]+c2*(ue[m_constT4+ip1]-q[ip1]))-
						(ue[m_const+im1]*buf[m_const+im1]+c2*(ue[m_constT4+im1]-q[im1])))+
					xxcon1*(buf[m_const+ip1]-2.0*buf[m_const+i]+buf[m_const+im1])+
					dx2tx1*(ue[m_const+ip1]-2.0*ue[m_const+i]+ue[m_const+im1]);
				forcing[ku4+ku3+ku2+2]=forcing[ku4+ku3+ku2+2]-tx2*(
						ue[m_constT2+ip1]*buf[m_const+ip1]-ue[m_constT2+im1]*buf[m_const+im1])+
					xxcon2*(buf[m_constT2+ip1]-2.0*buf[m_constT2+i]+buf[m_constT2+im1])+
					dx3tx1*(ue[m_constT2+ip1]-2.0*ue[m_constT2+i] +ue[m_constT2+im1]);
				forcing[ku4+ku3+ku2+3]=forcing[ku4+ku3+ku2+3]-tx2*(
						ue[m_constT3+ip1]*buf[m_const+ip1]-ue[m_constT3+im1]*buf[m_const+im1])+
					xxcon2*(buf[m_constT3+ip1]-2.0*buf[m_constT3+i]+buf[m_constT3+im1])+
					dx4tx1*(ue[m_constT3+ip1]-2.0*ue[m_constT3+i]+ue[m_constT3+im1]);
				forcing[ku4+ku3+ku2+4]=forcing[ku4+ku3+ku2+4]-tx2*(
						buf[m_const+ip1]*(c1*ue[m_constT4+ip1]-c2*q[ip1])-
						buf[m_const+im1]*(c1*ue[m_constT4+im1]-c2*q[im1]))+
					0.5*xxcon3*(buf[ip1]-2.0*buf[i]+
							buf[im1])+
					xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
					xxcon5*(buf[m_constT4+ip1]-2.0*buf[m_constT4+i]+buf[m_constT4+im1])+
					dx5tx1*(ue[m_constT4+ip1]-2.0*ue[m_constT4+i]+ue[m_constT4+im1]);
			}
			/* 
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                         
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				i=1;
				forcing[ku4+ku3+5+m]=forcing[ku4+ku3+5+m]-dssp*
					(5.0*ue[m_calc+i]-4.0*ue[m_calc+i+1]+ue[m_calc+i+2]);
				i=2;
				forcing[ku4+ku3+10+m]=forcing[ku4+ku3+10+m]-dssp*
					(-4.0*ue[m_calc+i-1]+6.0*ue[m_calc+i]-
					 4.0*ue[m_calc+i+1]+ue[m_calc+i+2]);
			}
			for(int i=3; i<=grid_points[0]-4; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					m_calc = m*m_const;
					forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
						(ue[m_calc+i-2]-4.0*ue[m_calc+i-1]+
						 6.0*ue[m_calc+i]-4.0*ue[m_calc+i+1]+ue[m_calc+i+2]);
				}
			}
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				i=grid_points[0]-3;
				ku2 = i*5;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(ue[m_calc+i-2]-4.0*ue[m_calc+i-1]+
					 6.0*ue[m_calc+i]-4.0*ue[m_calc+i+1]);
				i=grid_points[0]-2;
				forcing[ku4+ku3+ku2+5+m]=forcing[ku4+ku3+ku2+5+m]-dssp*
					(ue[m_calc+i-2]-4.0*ue[m_calc+i-1]+5.0*ue[m_calc+i]);
			}
		}
	}
	/* 
	 * ---------------------------------------------------------------------
	 * eta-direction flux differences             
	 * ---------------------------------------------------------------------
	 */
	for(int k=1; k<=grid_points[2]-2; k++){
		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			for(int j=0; j<=grid_points[1]-1; j++){
				eta=(double)(j)*dnym1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
					m_calc = m*m_const;
					ue[m_calc+j]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
					m_calc = m*m_const;
					buf[m_calc+j]=dtpp*dtemp[m];
				}
				cuf[j]=buf[m_constT2+j]*buf[m_constT2+j];
				buf[j]=cuf[j]+buf[m_const+j]*buf[m_const+j]+buf[m_constT3+j]*buf[m_constT3+j];
				q[j]=0.5*(buf[m_const+j]*ue[m_const+j]+buf[m_constT2+j]*ue[m_constT2+j]+
						buf[m_constT3+j]*ue[m_constT3+j]);
			}
			for(int j=1; j<=grid_points[1]-2; j++){
				ku3 = j*ku3_const;
				jm1=j-1;
				jp1=j+1;
				forcing[ku4+ku3+ku2+0]=forcing[ku4+ku3+ku2+0]-
					ty2*(ue[m_constT2+jp1]-ue[m_constT2+jm1])+
					dy1ty1*(ue[jp1]-2.0*ue[j]+ue[jm1]);
				forcing[ku4+ku3+ku2+1]=forcing[ku4+ku3+ku2+1]-ty2*(
						ue[m_const+jp1]*buf[m_constT2+jp1]-ue[m_const+jm1]*buf[m_constT2+jm1])+
					yycon2*(buf[m_const+jp1]-2.0*buf[m_const+j]+buf[m_const+jm1])+
					dy2ty1*(ue[m_const+jp1]-2.0*ue[m_const+j]+ue[m_const+jm1]);
				forcing[ku4+ku3+ku2+2]=forcing[ku4+ku3+ku2+2]-ty2*(
						(ue[m_constT2+jp1]*buf[m_constT2+jp1]+c2*(ue[m_constT4+jp1]-q[jp1]))-
						(ue[m_constT2+jm1]*buf[m_constT2+jm1]+c2*(ue[m_constT4+jm1]-q[jm1])))+
					yycon1*(buf[m_constT2+jp1]-2.0*buf[m_constT2+j]+buf[m_constT2+jm1])+
					dy3ty1*(ue[m_constT2+jp1]-2.0*ue[m_constT2+j]+ue[m_constT2+jm1]);
				forcing[ku4+ku3+ku2+3]=forcing[ku4+ku3+ku2+3]-ty2*(
						ue[m_constT3+jp1]*buf[m_constT2+jp1]-ue[m_constT3+jm1]*buf[m_constT2+jm1])+
					yycon2*(buf[m_constT3+jp1]-2.0*buf[m_constT3+j]+buf[m_constT3+jm1])+
					dy4ty1*(ue[m_constT3+jp1]-2.0*ue[m_constT3+j]+ue[m_constT3+jm1]);
				forcing[ku4+ku3+ku2+4]=forcing[ku4+ku3+ku2+4]-ty2*(
						buf[m_constT2+jp1]*(c1*ue[m_constT4+jp1]-c2*q[jp1])-
						buf[m_constT2+jm1]*(c1*ue[m_constT4+jm1]-c2*q[jm1]))+
					0.5*yycon3*(buf[jp1]-2.0*buf[j]+
							buf[jm1])+
					yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
					yycon5*(buf[m_constT4+jp1]-2.0*buf[m_constT4+j]+buf[m_constT4+jm1])+
					dy5ty1*(ue[m_constT4+jp1]-2.0*ue[m_constT4+j]+ue[m_constT4+jm1]);
			}
			/* 
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                      
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				j=1;
				forcing[ku4+ku3_const+ku2+m]=forcing[ku4+ku3_const+ku2+m]-dssp*
					(5.0*ue[m_calc+j]-4.0*ue[m_calc+j+1] +ue[m_calc+j+2]);
				j=2;
				ku3 = j*ku3_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(-4.0*ue[m_calc+j-1]+6.0*ue[m_calc+j]-
					 4.0*ue[m_calc+j+1]+ue[m_calc+j+2]);
			}
			for(int j=3; j<=grid_points[1]-4; j++){
				ku3 = j*ku3_const;
				for(int m=0; m<5; m++){
					m_calc = m*m_const;
					forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
						(ue[m_calc+j-2]-4.0*ue[m_calc+j-1]+
						 6.0*ue[m_calc+j]-4.0*ue[m_calc+j+1]+ue[m_calc+j+2]);
				}
			}
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				j=grid_points[1]-3;
				ku3 = j*ku3_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(ue[m_calc+j-2]-4.0*ue[m_calc+j-1]+
					 6.0*ue[m_calc+j]-4.0*ue[m_calc+j+1]);
				j=grid_points[1]-2;
				ku3 += ku3_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(ue[m_calc+j-2]-4.0*ue[m_calc+j-1]+5.0*ue[m_calc+j]);
			}
		}
	}
	/* 
	 * ---------------------------------------------------------------------
	 * zeta-direction flux differences                      
	 * ---------------------------------------------------------------------
	 */
	for(int j=1; j<=grid_points[1]-2; j++){
		ku3 = j*ku3_const;
		eta=(double)(j)*dnym1;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			for(int k=0; k<=grid_points[2]-1; k++){
				zeta=(double)(k)*dnzm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
					m_calc = m*m_const;
					ue[m_calc+k]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
					m_calc = m*m_const;
					buf[m_calc+k]=dtpp*dtemp[m];
				}
				cuf[k]=buf[m_constT3+k]*buf[m_constT3+k];
				buf[k]=cuf[k]+buf[m_const+k]*buf[m_const+k]+buf[m_constT2+k]*buf[m_constT2+k];
				q[k]=0.5*(buf[m_const+k]*ue[m_const+k]+buf[m_constT2+k]*ue[m_constT2+k]+
						buf[m_constT3+k]*ue[m_constT3+k]);
			}
			for(int k=1; k<=grid_points[2]-2; k++){
				ku4 = k*ku4_const;
				km1=k-1;
				kp1=k+1;
				forcing[ku4+ku3+ku2+0]=forcing[ku4+ku3+ku2+0]-
					tz2*(ue[m_constT3+kp1]-ue[m_constT3+km1])+
					dz1tz1*(ue[kp1]-2.0*ue[k]+ue[km1]);
				forcing[ku4+ku3+ku2+1]=forcing[ku4+ku3+ku2+1]-tz2*(
						ue[m_const+kp1]*buf[m_constT3+kp1]-ue[m_const+km1]*buf[m_constT3+km1])+
					zzcon2*(buf[m_const+kp1]-2.0*buf[m_const+k]+buf[m_const+km1])+
					dz2tz1*(ue[m_const+kp1]-2.0*ue[m_const+k]+ue[m_const+km1]);
				forcing[ku4+ku3+ku2+2]=forcing[ku4+ku3+ku2+2]-tz2*(
						ue[m_constT2+kp1]*buf[m_constT3+kp1]-ue[m_constT2+km1]*buf[m_constT3+km1])+
					zzcon2*(buf[m_constT2+kp1]-2.0*buf[m_constT2+k]+buf[m_constT2+km1])+
					dz3tz1*(ue[m_constT2+kp1]-2.0*ue[m_constT2+k]+ue[m_constT2+km1]);
				forcing[ku4+ku3+ku2+3]=forcing[ku4+ku3+ku2+3]-tz2*(
						(ue[m_constT3+kp1]*buf[m_constT3+kp1]+c2*(ue[m_constT4+kp1]-q[kp1]))-
						(ue[m_constT3+km1]*buf[m_constT3+km1]+c2*(ue[m_constT4+km1]-q[km1])))+
					zzcon1*(buf[m_constT3+kp1]-2.0*buf[m_constT3+k]+buf[m_constT3+km1])+
					dz4tz1*(ue[m_constT3+kp1]-2.0*ue[m_constT3+k]+ue[m_constT3+km1]);
				forcing[ku4+ku3+ku2+4]=forcing[ku4+ku3+ku2+4]-tz2*(
						buf[m_constT3+kp1]*(c1*ue[m_constT4+kp1]-c2*q[kp1])-
						buf[m_constT3+km1]*(c1*ue[m_constT4+km1]-c2*q[km1]))+
					0.5*zzcon3*(buf[kp1]-2.0*buf[k]
							+buf[km1])+
					zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
					zzcon5*(buf[m_constT4+kp1]-2.0*buf[m_constT4+k]+buf[m_constT4+km1])+
					dz5tz1*(ue[m_constT4+kp1]-2.0*ue[m_constT4+k]+ue[m_constT4+km1]);
			}
			/* 
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                        
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				k=1;
				ku4 = ku4_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(5.0*ue[m_calc+k]-4.0*ue[m_calc+k+1]+ue[m_calc+k+2]);
				k=2;
				ku4 += ku4_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(-4.0*ue[m_calc+k-1]+6.0*ue[m_calc+k]-
					 4.0*ue[m_calc+k+1]+ue[m_calc+k+2]);
			}
			for(int k=3; k<=grid_points[2]-4; k++){
				ku4 = k*ku4_const;
				for(int m=0; m<5; m++){
					m_calc = m*m_const;
					forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
						(ue[m_calc+k-2]-4.0*ue[m_calc+k-1]+
						 6.0*ue[m_calc+k]-4.0*ue[m_calc+k+1]+ue[m_calc+k+2]);
				}
			}
			for(int m=0; m<5; m++){
				m_calc = m*m_const;
				k=grid_points[2]-3;
				ku4 = k*ku4_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(ue[m_calc+k-2]-4.0*ue[m_calc+k-1]+
					 6.0*ue[m_calc+k]-4.0*ue[m_calc+k+1]);
				k=grid_points[2]-2;
				ku4 += ku4_const;
				forcing[ku4+ku3+ku2+m]=forcing[ku4+ku3+ku2+m]-dssp*
					(ue[m_calc+k-2]-4.0*ue[m_calc+k-1]+5.0*ue[m_calc+k]);
			}
		}
	}
	/* 
	 * ---------------------------------------------------------------------
	 * now change the sign of the forcing function
	 * ---------------------------------------------------------------------
	 */
	std::transform(policy, forcing.begin(), forcing.end(), forcing.begin(), std::negate<double>());
}

/* 
 * ---------------------------------------------------------------------
 * this function returns the exact solution at point xi, eta, zeta  
 * ---------------------------------------------------------------------
 */
void exact_solution(double xi, double eta, double zeta, std::span<double> dtemp){
	int m;
	for(m=0; m<5; m++){
		dtemp[m]=ce[(0*5)+m]+
			xi*(ce[(1*5)+m]+
					xi*(ce[(4*5)+m]+
						xi*(ce[(7*5)+m]+
							xi*ce[(10*5)+m])))+
			eta*(ce[(2*5)+m]+
					eta*(ce[(5*5)+m]+
						eta*(ce[(8*5)+m]+
							eta*ce[(11*5)+m])))+
			zeta*(ce[(3*5)+m]+
					zeta*(ce[(6*5)+m]+
						zeta*(ce[(9*5)+m]+ 
							zeta*ce[(12*5)+m])));
	}
}

/* 
 * ---------------------------------------------------------------------
 * this subroutine initializes the field variable u using 
 * tri-linear transfinite interpolation of the boundary values     
 * ---------------------------------------------------------------------
 */
void initialize(){
	int i, j, k;
	double xi, eta, zeta, Pxi, Peta, Pzeta;
	std::vector<double> Pface(30), temp(5);
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	/* 
	 * ---------------------------------------------------------------------
	 * later (in compute_rhs) we compute 1/u for every element. a few of 
	 * the corner elements are not used, but it convenient (and faster) 
	 * to compute the whole thing with a simple loop. make sure those 
	 * values are nonzero by initializing the whole thing here. 
	 * ---------------------------------------------------------------------
	 */
	std::fill(policy, u.begin(), u.end(), 1.0);
	/* 
	 * ---------------------------------------------------------------------
	 * first store the "interpolated" values everywhere on the grid    
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku3, ku2;
		std::vector<double> Pface(2*3*5, 0);
		double xi, eta, zeta, Pxi, Peta, Pzeta;

		ku4 = k*ku4_const;
		zeta=(double)(k)* dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
			ku3 = j*ku3_const;
			eta=(double)(j)*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
				ku2 = i*5;
				xi=(double)(i)*dnxm1;
				for(int ix=0; ix<2; ix++){
					exact_solution((double)ix, eta, zeta, {Pface.begin() + (ix*15), Pface.end()});
				}
				for(int iy=0; iy<2; iy++){
					exact_solution(xi, (double)iy , zeta, {Pface.begin() + (iy*15+5), Pface.end()});
				}
				for(int iz=0; iz<2; iz++){
					exact_solution(xi, eta, (double)iz, {Pface.begin() + (iz*15+10), Pface.end()});
				}
				for(int m=0; m<5; m++){
					Pxi=xi*Pface[15+m]+(1.0-xi)*Pface[m];
					Peta=eta*Pface[20+m]+(1.0-eta)*Pface[5+m];
					Pzeta=zeta*Pface[25+m]+(1.0-zeta)*Pface[10+m];
					u[ku4+ku3+ku2+m]=Pxi+Peta+Pzeta- 
						Pxi*Peta-Pxi*Pzeta-Peta*Pzeta+ 
						Pxi*Peta*Pzeta;
				}
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * now store the exact values on the boundaries        
	 * ---------------------------------------------------------------------
	 * west face                                                  
	 * ---------------------------------------------------------------------
	 */
	i=0;
	xi=0.0;
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku3, ku2;
  		std::vector<double> temp(5, 0);
		double zeta, eta;
		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
			ku3 = j*ku3_const;
			eta=(double)(j)*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+m]=temp[m];
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * east face                                                      
	 * ---------------------------------------------------------------------
	 */
	i=grid_points[0]-1;
	ku2 = i*5;
	xi=1.0;
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku3;
  		std::vector<double> temp(5, 0);
		double zeta, eta;
		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
			ku3 = j*ku3_const;
			eta=(double)(j)*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * south face                                                 
	 * ---------------------------------------------------------------------
	 */
	j=0;
	eta=0.0;
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku2;
		std::vector<double> temp(5, 0);						
		double zeta, xi;

		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int i=0; i<=grid_points[0]-1; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku2+m]=temp[m];
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * north face                                    
	 * ---------------------------------------------------------------------
	 */
	j=grid_points[1]-1;
	eta=1.0;
	ku3 = j*ku3_const;
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku2;
		std::vector<double> temp(5, 0);
		double zeta, xi;

		ku4 = k*ku4_const;
		zeta=(double)(k)*dnzm1;
		for(int i=0; i<=grid_points[0]-1; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * bottom face                                       
	 * ---------------------------------------------------------------------
	 */
	k=0;
	zeta=0.0;
	std::for_each_n(policy, iter.front(), grid_points[1], [&](int j){
		int ku3, ku2;
		std::vector<double> temp(5, 0);
		double eta, xi;

		ku3 = j*ku3_const;
		eta=(double)(j)*dnym1;
		for(int i=0; i<=grid_points[0]-1; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku3+ku2+m]=temp[m];
			}
		}
	});
	/* 
	 * ---------------------------------------------------------------------
	 * top face     
	 * ---------------------------------------------------------------------
	 */
	k=grid_points[2]-1;
	zeta=1.0;
	ku4 = k*ku4_const;
	std::for_each_n(policy, iter.front(), grid_points[1], [&](int j){
		int ku3, ku2;
		std::vector<double> temp(5, 0);
		double eta, xi;

		ku3 = j*ku3_const;
		eta=(double)(j)*dnym1;
		for(int i=0; i<=grid_points[0]-1; i++){
			ku2 = i*5;
			xi=(double)(i)*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	});
}

void lhsinit(std::vector<double>& lhs, int size){
	int i=size;
	int lh4=i*75, lh3T1=25, lh3T2=50, lh2;
	/* 
	 * ---------------------------------------------------------------------
	 * zero the whole left hand side for starters
	 * ---------------------------------------------------------------------
	 */
	for(int m=0; m<5; m++){
		for(int n=0; n<5; n++){		
			lh2 = n*5;		
			lhs[lh2+m]=0.0;
			lhs[lh3T1+lh2+m]=0.0;
			lhs[lh3T2+lh2+m]=0.0;
			lhs[lh4+lh2+m]=0.0;
			lhs[lh4+lh3T1+lh2+m]=0.0;
			lhs[lh4+lh3T2+lh2+m]=0.0;
		}
	}
	/* 
	 * ---------------------------------------------------------------------
	 * next, set all diagonal values to 1. This is overkill, but convenient
	 * ---------------------------------------------------------------------
	 */
	for(int m=0; m<5; m++){
		lh2 = m*5;
		lhs[lh3T1+lh2+m]=1.0;
		lhs[lh4+lh3T1+lh2+m]=1.0;
	}
}

/*
 * ---------------------------------------------------------------------
 * subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
 * ---------------------------------------------------------------------
 */
void matmul_sub(double ablock[], double bblock[], double cblock[]){
	 for(int i=0; i<5; i++){
		for(int j=0; j<5; j++){
			for(int k=0; k<5; k++){
				cblock[i*5+j] -= ablock[k*5+j]*bblock[i*5+k];
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * subtracts bvec=bvec - ablock*avec
 * ---------------------------------------------------------------------
 */
void matvec_sub(double ablock[], double avec[5], double bvec[5]){
	/*
	 * ---------------------------------------------------------------------
	 * rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] - lhs[ia][ablock][0][i]*
	 * ---------------------------------------------------------------------
	 */

	for(int i=0; i<5; i++){
		for(int j=0; j<5; j++){
			bvec[i] -= ablock[j*5+i]*avec[j];
		}
	}
}

void rhs_norm(std::vector<double> &rms){
	int i, j, k, d, m;
	double add;
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	std::fill(rms.begin(), rms.end(), 0.0);

	for(int k=1; k<=grid_points[2]-2; k++){
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++) {
					add=rhs[ku4+ku3+ku2+m];
					rms[m]=rms[m]+add*add;
				} 
			} 
		} 
	}
	for(int m=0; m<5; m++){
		for(int d=0; d<3; d++){
			rms[m]=rms[m]/(double)(grid_points[d]-2);
		}
		rms[m]=sqrt(rms[m]);
	}
}

void set_constants(){
	ce[(0*5)]=2.0;
	ce[(1*5)]=0.0;
	ce[(2*5)]=0.0;
	ce[(3*5)]=4.0;
	ce[(4*5)]=5.0;
	ce[(5*5)]=3.0;
	ce[(6*5)]=0.5;
	ce[(7*5)]=0.02;
	ce[(8*5)]=0.01;
	ce[(9*5)]=0.03;
	ce[(10*5)]=0.5;
	ce[(11*5)]=0.4;
	ce[(12*5)]=0.3;
	/* */
	ce[(0*5)+1]=1.0;
	ce[(1*5)+1]=0.0;
	ce[(2*5)+1]=0.0;
	ce[(3*5)+1]=0.0;
	ce[(4*5)+1]=1.0;
	ce[(5*5)+1]=2.0;
	ce[(6*5)+1]=3.0;
	ce[(7*5)+1]=0.01;
	ce[(8*5)+1]=0.03;
	ce[(9*5)+1]=0.02;
	ce[(10*5)+1]=0.4;
	ce[(11*5)+1]=0.3;
	ce[(12*5)+1]=0.5;
	/* */
	ce[(0*5)+2]=2.0;
	ce[(1*5)+2]=2.0;
	ce[(2*5)+2]=0.0;
	ce[(3*5)+2]=0.0;
	ce[(4*5)+2]=0.0;
	ce[(5*5)+2]=2.0;
	ce[(6*5)+2]=3.0;
	ce[(7*5)+2]=0.04;
	ce[(8*5)+2]=0.03;
	ce[(9*5)+2]=0.05;
	ce[(10*5)+2]=0.3;
	ce[(11*5)+2]=0.5;
	ce[(12*5)+2]=0.4;
	/* */
	ce[(0*5)+3]=2.0;
	ce[(1*5)+3]=2.0;
	ce[(2*5)+3]=0.0;
	ce[(3*5)+3]=0.0;
	ce[(4*5)+3]=0.0;
	ce[(5*5)+3]=2.0;
	ce[(6*5)+3]=3.0;
	ce[(7*5)+3]=0.03;
	ce[(8*5)+3]=0.05;
	ce[(9*5)+3]=0.04;
	ce[(10*5)+3]=0.2;
	ce[(11*5)+3]=0.1;
	ce[(12*5)+3]=0.3;
	/* */
	ce[(0*5)+4]=5.0;
	ce[(1*5)+4]=4.0;
	ce[(2*5)+4]=3.0;
	ce[(3*5)+4]=2.0;
	ce[(4*5)+4]=0.1;
	ce[(5*5)+4]=0.4;
	ce[(6*5)+4]=0.3;
	ce[(7*5)+4]=0.05;
	ce[(8*5)+4]=0.04;
	ce[(9*5)+4]=0.03;
	ce[(10*5)+4]=0.1;
	ce[(11*5)+4]=0.3;
	ce[(12*5)+4]=0.2;
	/* */
	c1=1.4;
	c2=0.4;
	c3=0.1;
	c4=1.0;
	c5=1.4;
	dnxm1=1.0/(double)(grid_points[0]-1);
	dnym1=1.0/(double)(grid_points[1]-1);
	dnzm1=1.0/(double)(grid_points[2]-1);
	c1c2=c1*c2;
	c1c5=c1*c5;
	c3c4=c3*c4;
	c1345=c1c5*c3c4;
	conz1=(1.0-c1c5);
	tx1=1.0/(dnxm1*dnxm1);
	tx2=1.0/(2.0*dnxm1);
	tx3=1.0/dnxm1;
	ty1=1.0/(dnym1*dnym1);
	ty2=1.0/(2.0*dnym1);
	ty3=1.0/dnym1;
	tz1=1.0/(dnzm1*dnzm1);
	tz2=1.0/(2.0*dnzm1);
	tz3=1.0/dnzm1;
	dx1=0.75;
	dx2=0.75;
	dx3=0.75;
	dx4=0.75;
	dx5=0.75;
	dy1=0.75;
	dy2=0.75;
	dy3=0.75;
	dy4=0.75;
	dy5=0.75;
	dz1=1.0;
	dz2=1.0;
	dz3=1.0;
	dz4=1.0;
	dz5=1.0;
	dxmax=std::max(dx3, dx4);
	dymax=std::max(dy2, dy4);
	dzmax=std::max(dz2, dz3);
	dssp=0.25*std::max(dx1, std::max(dy1, dz1));
	c4dssp=4.0*dssp;
	c5dssp=5.0*dssp;
	dttx1=dt*tx1;
	dttx2=dt*tx2;
	dtty1=dt*ty1;
	dtty2=dt*ty2;
	dttz1=dt*tz1;
	dttz2=dt*tz2;
	c2dttx1=2.0*dttx1;
	c2dtty1=2.0*dtty1;
	c2dttz1=2.0*dttz1;
	dtdssp=dt*dssp;
	comz1=dtdssp;
	comz4=4.0*dtdssp;
	comz5=5.0*dtdssp;
	comz6=6.0*dtdssp;
	c3c4tx3=c3c4*tx3;
	c3c4ty3=c3c4*ty3;
	c3c4tz3=c3c4*tz3;
	dx1tx1=dx1*tx1;
	dx2tx1=dx2*tx1;
	dx3tx1=dx3*tx1;
	dx4tx1=dx4*tx1;
	dx5tx1=dx5*tx1;
	dy1ty1=dy1*ty1;
	dy2ty1=dy2*ty1;
	dy3ty1=dy3*ty1;
	dy4ty1=dy4*ty1;
	dy5ty1=dy5*ty1;
	dz1tz1=dz1*tz1;
	dz2tz1=dz2*tz1;
	dz3tz1=dz3*tz1;
	dz4tz1=dz4*tz1;
	dz5tz1=dz5*tz1;
	c2iv=2.5;
	con43=4.0/3.0;
	con16=1.0/6.0;
	xxcon1=c3c4tx3*con43*tx3;
	xxcon2=c3c4tx3*tx3;
	xxcon3=c3c4tx3*conz1*tx3;
	xxcon4=c3c4tx3*con16*tx3;
	xxcon5=c3c4tx3*c1c5*tx3;
	yycon1=c3c4ty3*con43*ty3;
	yycon2=c3c4ty3*ty3;
	yycon3=c3c4ty3*conz1*ty3;
	yycon4=c3c4ty3*con16*ty3;
	yycon5=c3c4ty3*c1c5*ty3;
	zzcon1=c3c4tz3*con43*tz3;
	zzcon2=c3c4tz3*tz3;
	zzcon3=c3c4tz3*conz1*tz3;
	zzcon4=c3c4tz3*con16*tz3;
	zzcon5=c3c4tz3*c1c5*tz3;
}

/*
 * ---------------------------------------------------------------------
 * verification routine                         
 * ---------------------------------------------------------------------
 */
void verify(int no_time_steps, char* class_npb, boolean* verified){
	std::vector<double> xce(5), xcr(5), xcrref(5), xceref(5), xcrdif(5), xcedif(5); 
	double epsilon, dtref=0.0;
	/*
	 * ---------------------------------------------------------------------
	 * tolerance level
	 * ---------------------------------------------------------------------
	 */  
	epsilon=1.0e-08;
	/*
	 * ---------------------------------------------------------------------
	 * compute the error norm and the residual norm, and exit if not printing
	 * ---------------------------------------------------------------------
	 */  
	error_norm(xce);
	compute_rhs();
	rhs_norm(xcr);
	std::transform(xcr.begin(), xcr.end(), xcr.begin(), [&](double x){return x/dt;});
	
	*class_npb='U';
	*verified=TRUE;
	std::fill(xcrref.begin(), xcrref.end(), 1.0);
	std::fill(xceref.begin(), xceref.end(), 1.0);
	/*
	 * ---------------------------------------------------------------------
	 * reference data for 12X12X12 grids after 60 time steps, with DT = 1.0e-02
	 * ---------------------------------------------------------------------
	 */  
	if((grid_points[0]==12)&&
			(grid_points[1]==12)&&
			(grid_points[2]==12)&&
			(no_time_steps==60)){
		*class_npb='S';
		dtref=1.0e-2;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */       
		xcrref[0]=1.7034283709541311e-01;
		xcrref[1]=1.2975252070034097e-02;
		xcrref[2]=3.2527926989486055e-02;
		xcrref[3]=2.6436421275166801e-02;
		xcrref[4]=1.9211784131744430e-01;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=4.9976913345811579e-04;
		xceref[1]=4.5195666782961927e-05;
		xceref[2]=7.3973765172921357e-05;
		xceref[3]=7.3821238632439731e-05;
		xceref[4]=8.9269630987491446e-04;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
		 * ---------------------------------------------------------------------
		 */  
	}else if((grid_points[0]==24)&&
			(grid_points[1]==24)&&
			(grid_points[2]==24)&&
			(no_time_steps==200)){
		*class_npb='W';
		dtref = 0.8e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=0.1125590409344e+03;
		xcrref[1]=0.1180007595731e+02;
		xcrref[2]=0.2710329767846e+02;
		xcrref[3]=0.2469174937669e+02;
		xcrref[4]=0.2638427874317e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=0.4419655736008e+01;
		xceref[1]=0.4638531260002e+00;
		xceref[2]=0.1011551749967e+01;
		xceref[3]=0.9235878729944e+00;
		xceref[4]=0.1018045837718e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==64)&&
			(grid_points[1]==64)&&
			(grid_points[2]==64)&&
			(no_time_steps==200)){
		*class_npb='A';
		dtref=0.8e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=1.0806346714637264e+02;
		xcrref[1]=1.1319730901220813e+01;
		xcrref[2]=2.5974354511582465e+01;
		xcrref[3]=2.3665622544678910e+01;
		xcrref[4]=2.5278963211748344e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=4.2348416040525025e+00;
		xceref[1]=4.4390282496995698e-01;
		xceref[2]=9.6692480136345650e-01;
		xceref[3]=8.8302063039765474e-01;
		xceref[4]=9.7379901770829278e+00;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 102X102X102 grids after 200 time steps,
		 * with DT = 3.0e-04
		 * ---------------------------------------------------------------------
		 */  
	}else if((grid_points[0]==102)&&
			(grid_points[1]==102)&&
			(grid_points[2]==102)&&
			(no_time_steps==200)){
		*class_npb='B';
		dtref=3.0e-4;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=1.4233597229287254e+03;
		xcrref[1]=9.9330522590150238e+01;
		xcrref[2]=3.5646025644535285e+02;
		xcrref[3]=3.2485447959084092e+02;
		xcrref[4]=3.2707541254659363e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=5.2969847140936856e+01;
		xceref[1]=4.4632896115670668e+00;
		xceref[2]=1.3122573342210174e+01;
		xceref[3]=1.2006925323559144e+01;
		xceref[4]=1.2459576151035986e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 162X162X162 grids after 200 time steps,
		 * with DT = 1.0e-04
		 * ---------------------------------------------------------------------
		 */  
	}else if((grid_points[0]==162)&&
			(grid_points[1]==162)&&
			(grid_points[2]==162)&&
			(no_time_steps==200)){
		*class_npb='C';
		dtref=1.0e-4;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=0.62398116551764615e+04;
		xcrref[1]=0.50793239190423964e+03;
		xcrref[2]=0.15423530093013596e+04;
		xcrref[3]=0.13302387929291190e+04;
		xcrref[4]=0.11604087428436455e+05;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=0.16462008369091265e+03;
		xceref[1]=0.11497107903824313e+02;
		xceref[2]=0.41207446207461508e+02;
		xceref[3]=0.37087651059694167e+02;
		xceref[4]=0.36211053051841265e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 408x408x408 grids after 250 time steps,
		 * with DT = 0.2e-04
		 * ---------------------------------------------------------------------
		 */ 
	}else if((grid_points[0]==408)&&
			(grid_points[1]==408)&&
			(grid_points[2]==408)&&
			(no_time_steps==250)){
		*class_npb='D';
		dtref=0.2e-4;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=0.2533188551738e+05;
		xcrref[1]=0.2346393716980e+04;
		xcrref[2]=0.6294554366904e+04;
		xcrref[3]=0.5352565376030e+04;
		xcrref[4]=0.3905864038618e+05;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=0.3100009377557e+03;
		xceref[1]=0.2424086324913e+02;
		xceref[2]=0.7782212022645e+02;
		xceref[3]=0.6835623860116e+02;
		xceref[4]=0.6065737200368e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 1020x1020x1020 grids after 250 time steps,
		 * with DT = 0.4e-05
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==1020)&&
			(grid_points[1]==1020)&&
			(grid_points[2]==1020)&&
			(no_time_steps==250)){
		*class_npb='E';
		dtref=0.4e-5;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */    
		xcrref[0]=0.9795372484517e+05;
		xcrref[1]=0.9739814511521e+04;
		xcrref[2]=0.2467606342965e+05;
		xcrref[3]=0.2092419572860e+05;
		xcrref[4]=0.1392138856939e+06;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */    
		xceref[0]=0.4327562208414e+03;
		xceref[1]=0.3699051964887e+02;
		xceref[2]=0.1089845040954e+03;
		xceref[3]=0.9462517622043e+02;
		xceref[4]=0.7765512765309e+03;
	}else{
		*verified=FALSE;
	}
	/*
	 * ---------------------------------------------------------------------
	 * verification test for residuals if gridsize is one of 
	 * the defined grid sizes above (*class_npb != 'U')
	 * ---------------------------------------------------------------------
	 * compute the difference of solution values and the known reference values.
	 * ---------------------------------------------------------------------
	 */  
	for(int m=0; m<5; m++){
		xcrdif[m]=fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m]=fabs((xce[m]-xceref[m])/xceref[m]);
	}
	/*
	 * ---------------------------------------------------------------------
	 * output the comparison of computed results to known cases.
	 * ---------------------------------------------------------------------
	 */ 
	if(*class_npb!='U'){
    std::cout
      << " Verification being performed for class "
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
      << " Unknown class"
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

/*
 * ---------------------------------------------------------------------
 * performs line solves in X direction by first factoring
 * the block-tridiagonal matrix into an upper triangular matrix, 
 * and then performing back substitution to solve for the unknow
 * vectors of each line.  
 * 
 * make sure we treat elements zero to cell_size in the direction
 * of the sweep. 
 * ---------------------------------------------------------------------
 */
void x_solve(){
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	int ks3, ks2;
	int ks3_const=(JMAXP+1)*(IMAXP+1), ks2_const=(IMAXP+1);

	int lh4, lh3, lh2;
	int lh4_const=75, lh3_const=25;
	int lh3AA=lh3_const*AA, lh3BB=lh3_const*BB, lh3CC=lh3_const*CC;

	int kj3, kj2;
	int kj3_const=25;

	if(timeron){timer_start(T_XSOLVE);}
	/*
	 * ---------------------------------------------------------------------
	 * this function computes the left hand side in the xi-direction
	 * ---------------------------------------------------------------------
	 */
	int isize=grid_points[0]-1;
	/*
	 * ---------------------------------------------------------------------
	 * determine a (labeled f) and n jacobians
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2, kj3, kj2;
		int lh4, lh3, lh2;
		std::vector<double> fjac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double> njac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double>  lhs((PROBLEM_SIZE+1)*3*5*5, 0.0);
    	double tmp1, tmp2, tmp3;

		ks3=k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ks2=j*ks2_const;
			ku3 = j*ku3_const;
			for(int i=0; i<=isize; i++){
				kj3 = i*kj3_const;
				ku2 = i*5;
				tmp1=rho_i[ks3+ks2+i];
				tmp2=tmp1*tmp1;
				tmp3=tmp1*tmp2;
				fjac[kj3+(0*5)+0]=0.0;
				fjac[kj3+(1*5)+0]=1.0;
				fjac[kj3+(2*5)+0]=0.0;
				fjac[kj3+(3*5)+0]=0.0;
				fjac[kj3+(4*5)+0]=0.0;
				fjac[kj3+(0*5)+1]=-(u[ku4+ku3+ku2+1]*tmp2*u[ku4+ku3+ku2+1])+c2*qs[ks3+ks2+i];
				fjac[kj3+(1*5)+1]=(2.0-c2)*(u[ku4+ku3+ku2+1]/u[ku4+ku3+ku2]);
				fjac[kj3+(2*5)+1]=-c2*(u[ku4+ku3+ku2+2]*tmp1);
				fjac[kj3+(3*5)+1]=-c2*(u[ku4+ku3+ku2+3]*tmp1);
				fjac[kj3+(4*5)+1]=c2;
				fjac[kj3+(0*5)+2]=-(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+2])*tmp2;
				fjac[kj3+(1*5)+2]=u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+(2*5)+2]=u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+(3*5)+2]=0.0;
				fjac[kj3+(4*5)+2]=0.0;
				fjac[kj3+(0*5)+3]=-(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+(1*5)+3]=u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+(2*5)+3]=0.0;
				fjac[kj3+(3*5)+3]=u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+(4*5)+3]=0.0;
				fjac[kj3+(0*5)+4]=(c2*2.0*square[ks3+ks2+i]-c1*u[ku4+ku3+ku2+4])*(u[ku4+ku3+ku2+1]*tmp2);
				fjac[kj3+(1*5)+4]=c1*u[ku4+ku3+ku2+4]*tmp1-c2*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]*tmp2+qs[ks3+ks2+i]);
				fjac[kj3+(2*5)+4]=-c2*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+1])*tmp2;
				fjac[kj3+(3*5)+4]=-c2*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+1])*tmp2;
				fjac[kj3+(4*5)+4]=c1*(u[ku4+ku3+ku2+1]*tmp1);
				njac[kj3+(0*5)+0]=0.0;
				njac[kj3+(1*5)+0]=0.0;
				njac[kj3+(2*5)+0]=0.0;
				njac[kj3+(3*5)+0]=0.0;
				njac[kj3+(4*5)+0]=0.0;
				njac[kj3+(0*5)+1]=-con43*c3c4*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+(1*5)+1]=con43*c3c4*tmp1;
				njac[kj3+(2*5)+1]=0.0;
				njac[kj3+(3*5)+1]=0.0;
				njac[kj3+(4*5)+1]=0.0;
				njac[kj3+(0*5)+2]=-c3c4*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+(1*5)+2]=0.0;
				njac[kj3+(2*5)+2]=c3c4*tmp1;
				njac[kj3+(3*5)+2]=0.0;
				njac[kj3+(4*5)+2]=0.0;
				njac[kj3+(0*5)+3]=-c3c4*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+(1*5)+3]=0.0;
				njac[kj3+(2*5)+3]=0.0;
				njac[kj3+(3*5)+3]=c3c4*tmp1;
				njac[kj3+(4*5)+3]=0.0;
				njac[kj3+(0*5)+4]=-(con43*c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1])
					-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2])
					-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					-c1345*tmp2*u[ku4+ku3+ku2+4];
				njac[kj3+(1*5)+4]=(con43*c3c4-c1345)*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+(2*5)+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+(3*5)+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+(4*5)+4]=(c1345)*tmp1;
			}
			/*
			 * ---------------------------------------------------------------------
			 * now jacobians set, so form left hand side in x direction
			 * ---------------------------------------------------------------------
			 */
			lhsinit(lhs, isize);
			int lh2T2=10, lh2T3=15, lh2T4=20;
			for(int i=1; i<=isize-1; i++){
				lh4 = i*lh4_const;
				lh3 = AA*lh3_const;
				kj3 = (i-1)*kj3_const;
				tmp1=dt*tx1;
				tmp2=dt*tx2;
				lhs[lh4+lh3]=-tmp2*fjac[kj3+(0*5)+0]
					-tmp1*njac[kj3+(0*5)+0]
					-tmp1*dx1; 
				lhs[lh4+lh3+(1*5)+0]=-tmp2*fjac[kj3+(1*5)+0]
					-tmp1*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=-tmp2*fjac[kj3+(2*5)+0]
					-tmp1*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=-tmp2*fjac[kj3+(3*5)+0]
					-tmp1*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=-tmp2*fjac[kj3+(4*5)+0]
					-tmp1*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=-tmp2*fjac[kj3+(0*5)+1]
					-tmp1*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=-tmp2*fjac[kj3+(1*5)+1]
					-tmp1*njac[kj3+(1*5)+1]
					-tmp1*dx2;
				lhs[lh4+lh3+(2*5)+1]=-tmp2*fjac[kj3+(2*5)+1]
					-tmp1*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=-tmp2*fjac[kj3+(3*5)+1]
					-tmp1*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=-tmp2*fjac[kj3+(4*5)+1]
					-tmp1*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=-tmp2*fjac[kj3+(0*5)+2]
					-tmp1*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=-tmp2*fjac[kj3+(1*5)+2]
					-tmp1*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=-tmp2*fjac[kj3+(2*5)+2]
					-tmp1*njac[kj3+(2*5)+2]
					-tmp1*dx3;
				lhs[lh4+lh3+(3*5)+2]=-tmp2*fjac[kj3+(3*5)+2]
					-tmp1*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=-tmp2*fjac[kj3+(4*5)+2]
					-tmp1*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=-tmp2*fjac[kj3+(0*5)+3]
					-tmp1*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=-tmp2*fjac[kj3+(1*5)+3]
					-tmp1*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=-tmp2*fjac[kj3+(2*5)+3]
					-tmp1*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=-tmp2*fjac[kj3+(3*5)+3]
					-tmp1*njac[kj3+(3*5)+3]
					-tmp1*dx4;
				lhs[lh4+lh3+(4*5)+3]=-tmp2*fjac[kj3+(4*5)+3]
					-tmp1*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=-tmp2*fjac[kj3+(0*5)+4]
					-tmp1*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=-tmp2*fjac[kj3+(1*5)+4]
					-tmp1*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=-tmp2*fjac[kj3+(2*5)+4]
					-tmp1*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=-tmp2*fjac[kj3+(3*5)+4]
					-tmp1*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=-tmp2*fjac[kj3+(4*5)+4]
					-tmp1*njac[kj3+(4*5)+4]
					-tmp1*dx5;
				lh3 = BB*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3+(0*5)+0]=1.0
					+tmp1*2.0*njac[kj3+(0*5)+0]
					+tmp1*2.0*dx1;
				lhs[lh4+lh3+(1*5)+0]=tmp1*2.0*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=tmp1*2.0*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=tmp1*2.0*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=tmp1*2.0*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=tmp1*2.0*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=1.0
					+tmp1*2.0*njac[kj3+(1*5)+1]
					+tmp1*2.0*dx2;
				lhs[lh4+lh3+(2*5)+1]=tmp1*2.0*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=tmp1*2.0*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=tmp1*2.0*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=tmp1*2.0*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=tmp1*2.0*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=1.0
					+tmp1*2.0*njac[kj3+(2*5)+2]
					+tmp1*2.0*dx3;
				lhs[lh4+lh3+(3*5)+2]=tmp1*2.0*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=tmp1*2.0*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=tmp1*2.0*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=tmp1*2.0*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=tmp1*2.0*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=1.0
					+tmp1*2.0*njac[kj3+(3*5)+3]
					+tmp1*2.0*dx4;
				lhs[lh4+lh3+(4*5)+3]=tmp1*2.0*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=tmp1*2.0*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=tmp1*2.0*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=tmp1*2.0*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=tmp1*2.0*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=1.0
					+tmp1*2.0*njac[kj3+(4*5)+4]
					+tmp1*2.0*dx5;
				lh3 = CC*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3+(0*5)+0]=tmp2*fjac[kj3+(0*5)+0]
					-tmp1*njac[kj3+(0*5)+0]
					-tmp1*dx1;
				lhs[lh4+lh3+(1*5)+0]=tmp2*fjac[kj3+(1*5)+0]
					-tmp1*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=tmp2*fjac[kj3+(2*5)+0]
					-tmp1*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=tmp2*fjac[kj3+(3*5)+0]
					-tmp1*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=tmp2*fjac[kj3+(4*5)+0]
					-tmp1*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=tmp2*fjac[kj3+(0*5)+1]
					-tmp1*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=tmp2*fjac[kj3+(1*5)+1]
					-tmp1*njac[kj3+(1*5)+1]
					-tmp1*dx2;
				lhs[lh4+lh3+(2*5)+1]=tmp2*fjac[kj3+(2*5)+1]
					-tmp1*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=tmp2*fjac[kj3+(3*5)+1]
					-tmp1*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=tmp2*fjac[kj3+(4*5)+1]
					-tmp1*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=tmp2*fjac[kj3+(0*5)+2]
					-tmp1*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=tmp2*fjac[kj3+(1*5)+2]
					-tmp1*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=tmp2*fjac[kj3+(2*5)+2]
					-tmp1*njac[kj3+(2*5)+2]
					-tmp1*dx3;
				lhs[lh4+lh3+(3*5)+2]=tmp2*fjac[kj3+(3*5)+2]
					-tmp1*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=tmp2*fjac[kj3+(4*5)+2]
					-tmp1*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=tmp2*fjac[kj3+(0*5)+3]
					-tmp1*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=tmp2*fjac[kj3+(1*5)+3]
					-tmp1*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=tmp2*fjac[kj3+(2*5)+3]
					-tmp1*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=tmp2*fjac[kj3+(3*5)+3]
					-tmp1*njac[kj3+(3*5)+3]
					-tmp1*dx4;
				lhs[lh4+lh3+(4*5)+3]=tmp2*fjac[kj3+(4*5)+3]
					-tmp1*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=tmp2*fjac[kj3+(0*5)+4]
					-tmp1*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=tmp2*fjac[kj3+(1*5)+4]
					-tmp1*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=tmp2*fjac[kj3+(2*5)+4]
					-tmp1*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=tmp2 * fjac[kj3+(3*5)+4]
					-tmp1*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=tmp2*fjac[kj3+(4*5)+4]
					-tmp1*njac[kj3+(4*5)+4]
					-tmp1*dx5;
			}
			/*
			 * ---------------------------------------------------------------------
			 * performs guaussian elimination on this cell.
			 * 
			 * assumes that unpacking routines for non-first cells 
			 * preload C' and rhs' from previous cell.
			 * 
			 * assumed send happens outside this routine, but that
			 * c'(IMAX) and rhs'(IMAX) will be sent to next cell
			 * ---------------------------------------------------------------------
			 * outer most do loops - sweeping in i direction
			 * ---------------------------------------------------------------------
			 * multiply c(0,j,k) by b_inverse and copy back to c
			 * multiply rhs(0) by b_inverse(0) and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh3 = BB*lh3_const;
			binvcrhs(&lhs[lh3], &lhs[CC*lh3_const], &rhs[ku4+ku3]);
			/*
			 * ---------------------------------------------------------------------
			 * begin inner most do loop
			 * do all the elements of the cell unless last 
			 * ---------------------------------------------------------------------
			 */
			for(int i=1; i<=isize-1; i++){
				ku2 = i*5;
				lh4 = i*lh4_const;
				/*
				 * -------------------------------------------------------------------
				 * rhs(i) = rhs(i) - A*rhs(i-1)
				 * -------------------------------------------------------------------
				 */

				matvec_sub(&lhs[lh4+(lh3AA)], &rhs[ku4+ku3+ku2-5], &rhs[ku4+ku3+ku2]);
				/*
				 * -------------------------------------------------------------------
				 * B(i) = B(i) - C(i-1)*A(i)
				 * -------------------------------------------------------------------
				 */
				matmul_sub(&lhs[lh4+(lh3AA)], &lhs[(lh4-lh4_const)+CC*lh3_const], &lhs[lh4+lh3BB]);
				/*
				 * -------------------------------------------------------------------
				 * multiply c(i,j,k) by b_inverse and copy back to c
				 * multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
				 * -------------------------------------------------------------------
				 */
				binvcrhs(&lhs[lh4+lh3], &lhs[lh4+CC*lh3_const], &rhs[ku4+ku3+ku2]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * rhs(isize) = rhs(isize) - A*rhs(isize-1)
			 * ---------------------------------------------------------------------
			 */
			ku2 = isize*5;
			lh4 = isize*lh4_const;
			matvec_sub(&lhs[lh4+lh3AA], &rhs[ku4+ku3+ku2-5], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * B(isize) = B(isize) - C(isize-1)*A(isize)
			 * ---------------------------------------------------------------------
			 */
			matmul_sub(&lhs[lh4+lh3AA], &lhs[lh4-lh4_const+lh3CC], &lhs[lh4+lh3BB]);
			/*
			 * ---------------------------------------------------------------------
			 * multiply rhs() by b_inverse() and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh4 = isize*lh4_const;
			lh3 = BB*lh3_const;
			binvrhs(&lhs[lh4+lh3], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * back solve: if last cell, then generate U(isize)=rhs(isize)
			 * else assume U(isize) is loaded in un pack backsub_info
			 * so just use it
			 * after u(istart) will be sent to next cell
			 * ---------------------------------------------------------------------
			 */
			lh3 = CC*lh3_const;
			for(int i=isize-1; i>=0; i--){
				ku2 = i*5;
				lh4 = i*lh4_const;
				for(int m=0; m<BLOCK_SIZE; m++){
					for(int n=0; n<BLOCK_SIZE; n++){
						lh2 = n*5;
						rhs[ku4+ku3+ku2+m]-=lhs[lh4+lh3+lh2+m]*rhs[ku4+ku3+ku2+5+n];
					}
				}
			}
		}
	});
	if(timeron){timer_stop(T_XSOLVE);}
}

/*
 * ---------------------------------------------------------------------
 * performs line solves in y direction by first factoring
 * the block-tridiagonal matrix into an upper triangular matrix, 
 * and then performing back substitution to solve for the unknow
 * vectors of each line.  
 *  
 * make sure we treat elements zero to cell_size in the direction
 * of the sweep.
 * ---------------------------------------------------------------------
 */
void y_solve(){
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	int ks3, ks2;
	int ks3_const=(JMAXP+1)*(IMAXP+1), ks2_const=(IMAXP+1);

	int lh4, lh3, lh2;
	int lh4_const=75, lh3_const=25;
	int lh3AA=lh3_const*AA, lh3BB=lh3_const*BB, lh3CC=lh3_const*CC;

	int kj3, kj3_const=25;	

	if(timeron){timer_start(T_YSOLVE);}
	/*
	 * ---------------------------------------------------------------------
	 * this function computes the left hand side for the three y-factors   
	 * ---------------------------------------------------------------------
	 */
	int jsize=grid_points[1]-1;
	/*
	 * ---------------------------------------------------------------------
	 * compute the indices for storing the tri-diagonal matrix;
	 * determine a (labeled f) and n jacobians for cell c
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[2]-2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2, kj3, kj2;
		int lh4, lh3, lh2;
		std::vector<double> fjac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double> njac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double>  lhs((PROBLEM_SIZE+1)*3*5*5, 0.0);
    	double tmp1, tmp2, tmp3;

		ks3=k*ks3_const;
		ku4 = k*ku4_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int j=0; j<=jsize; j++){
				kj3 = j*kj3_const;
				ks2=j*ks2_const;
				ku3 = j*ku3_const;
				tmp1=rho_i[ks3+ks2+i];
				tmp2=tmp1*tmp1;
				tmp3=tmp1*tmp2;
				fjac[kj3+(0*5)+0]=0.0;
				fjac[kj3+(1*5)+0]=0.0;
				fjac[kj3+(2*5)+0]=1.0;
				fjac[kj3+(3*5)+0]=0.0;
				fjac[kj3+(4*5)+0]=0.0;
				fjac[kj3+(0*5)+1]=-(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+2])*tmp2;
				fjac[kj3+(1*5)+1]=u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+(2*5)+1]=u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+(3*5)+1]=0.0;
				fjac[kj3+(4*5)+1]=0.0;
				fjac[kj3+(0*5)+2]=-(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]*tmp2)+c2*qs[ks3+ks2+i];
				fjac[kj3+(1*5)+2]=-c2*u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+(2*5)+2]=(2.0-c2)*u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+(3*5)+2]=-c2*u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+(4*5)+2]=c2;
				fjac[kj3+(0*5)+3]=-(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+(1*5)+3]=0.0;
				fjac[kj3+(2*5)+3]=u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+(3*5)+3]=u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+(4*5)+3]=0.0;
				fjac[kj3+(0*5)+4]=(c2*2.0*square[ks3+ks2+i]-c1*u[ku4+ku3+ku2+4])*u[ku4+ku3+ku2+2]*tmp2;
				fjac[kj3+(1*5)+4]=-c2*u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+2]*tmp2;
				fjac[kj3+(2*5)+4]=c1*u[ku4+ku3+ku2+4]*tmp1-c2*(qs[ks3+ks2+i]+u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]*tmp2);
				fjac[kj3+(3*5)+4]=-c2*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+(4*5)+4]=c1*u[ku4+ku3+ku2+2]*tmp1;
				njac[kj3+(0*5)+0]=0.0;
				njac[kj3+(1*5)+0]=0.0;
				njac[kj3+(2*5)+0]=0.0;
				njac[kj3+(3*5)+0]=0.0;
				njac[kj3+(4*5)+0]=0.0;
				njac[kj3+(0*5)+1]=-c3c4*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+(1*5)+1]=c3c4*tmp1;
				njac[kj3+(2*5)+1]=0.0;
				njac[kj3+(3*5)+1]=0.0;
				njac[kj3+(4*5)+1]=0.0;
				njac[kj3+(0*5)+2]=-con43*c3c4*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+(1*5)+2]=0.0;
				njac[kj3+(2*5)+2]=con43*c3c4*tmp1;
				njac[kj3+(3*5)+2]=0.0;
				njac[kj3+(4*5)+2]=0.0;
				njac[kj3+(0*5)+3]=-c3c4*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+(1*5)+3]=0.0;
				njac[kj3+(2*5)+3]=0.0;
				njac[kj3+(3*5)+3]=c3c4*tmp1;
				njac[kj3+(4*5)+3]=0.0;
				njac[kj3+(0*5)+4]=-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1])
					-(con43*c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2])
					-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					-c1345*tmp2*u[ku4+ku3+ku2+4];
				njac[kj3+(1*5)+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+(2*5)+4]=(con43*c3c4-c1345)*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+(3*5)+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+(4*5)+4]=(c1345)*tmp1;
			}
			/*
			 * ---------------------------------------------------------------------
			 * now joacobians set, so form left hand side in y direction
			 * ---------------------------------------------------------------------
			 */
			lhsinit(lhs, jsize);
			for(int j=1; j<=jsize-1; j++){
				lh4 = j*lh4_const;
				lh3 = AA*lh3_const;
				kj3 = (j-1)*kj3_const;
				tmp1=dt*ty1;
				tmp2=dt*ty2;
				lhs[lh4+lh3+(0*5)+0]=-tmp2*fjac[kj3+(0*5)+0]
					-tmp1*njac[kj3+(0*5)+0]
					-tmp1*dy1; 
				lhs[lh4+lh3+(1*5)+0]=-tmp2*fjac[kj3+(1*5)+0]
					-tmp1*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=-tmp2*fjac[kj3+(2*5)+0]
					-tmp1*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=-tmp2*fjac[kj3+(3*5)+0]
					-tmp1*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=-tmp2*fjac[kj3+(4*5)+0]
					-tmp1*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=-tmp2*fjac[kj3+(0*5)+1]
					-tmp1*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=-tmp2*fjac[kj3+(1*5)+1]
					-tmp1*njac[kj3+(1*5)+1]
					-tmp1*dy2;
				lhs[lh4+lh3+(2*5)+1]=-tmp2*fjac[kj3+(2*5)+1]
					-tmp1*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=-tmp2*fjac[kj3+(3*5)+1]
					-tmp1*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=-tmp2*fjac[kj3+(4*5)+1]
					-tmp1*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=-tmp2*fjac[kj3+(0*5)+2]
					-tmp1*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=-tmp2*fjac[kj3+(1*5)+2]
					-tmp1*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=-tmp2*fjac[kj3+(2*5)+2]
					-tmp1*njac[kj3+(2*5)+2]
					-tmp1*dy3;
				lhs[lh4+lh3+(3*5)+2]=-tmp2*fjac[kj3+(3*5)+2]
					-tmp1*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=-tmp2*fjac[kj3+(4*5)+2]
					-tmp1*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=-tmp2*fjac[kj3+(0*5)+3]
					-tmp1*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=-tmp2*fjac[kj3+(1*5)+3]
					-tmp1*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=-tmp2*fjac[kj3+(2*5)+3]
					-tmp1*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=-tmp2*fjac[kj3+(3*5)+3]
					-tmp1*njac[kj3+(3*5)+3]
					-tmp1*dy4;
				lhs[lh4+lh3+(4*5)+3]=-tmp2*fjac[kj3+(4*5)+3]
					-tmp1*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=-tmp2*fjac[kj3+(0*5)+4]
					-tmp1*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=-tmp2*fjac[kj3+(1*5)+4]
					-tmp1*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=-tmp2*fjac[kj3+(2*5)+4]
					-tmp1*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=-tmp2*fjac[kj3+(3*5)+4]
					-tmp1*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=-tmp2*fjac[kj3+(4*5)+4]
					-tmp1*njac[kj3+(4*5)+4]
					-tmp1*dy5;
				lh3 = BB*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3+(0*5)+0]=1.0
					+tmp1*2.0*njac[kj3+(0*5)+0]
					+tmp1*2.0*dy1;
				lhs[lh4+lh3+(1*5)+0]=tmp1*2.0*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=tmp1*2.0*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=tmp1*2.0*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=tmp1*2.0*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=tmp1*2.0*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=1.0
					+tmp1*2.0*njac[kj3+(1*5)+1]
					+tmp1*2.0*dy2;
				lhs[lh4+lh3+(2*5)+1]=tmp1*2.0*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=tmp1*2.0*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=tmp1*2.0*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=tmp1*2.0*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=tmp1*2.0*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=1.0
					+tmp1*2.0*njac[kj3+(2*5)+2]
					+tmp1*2.0*dy3;
				lhs[lh4+lh3+(3*5)+2]=tmp1*2.0*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=tmp1*2.0*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=tmp1*2.0*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=tmp1*2.0*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=tmp1*2.0*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=1.0
					+tmp1*2.0*njac[kj3+(3*5)+3]
					+tmp1*2.0*dy4;
				lhs[lh4+lh3+(4*5)+3]=tmp1*2.0*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=tmp1*2.0*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=tmp1*2.0*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=tmp1*2.0*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=tmp1*2.0*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=1.0
					+tmp1*2.0*njac[kj3+(4*5)+4] 
					+tmp1*2.0*dy5;
				lh3 = CC*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3+(0*5)+0]=tmp2*fjac[kj3+(0*5)+0]
					-tmp1*njac[kj3+(0*5)+0]
					-tmp1*dy1;
				lhs[lh4+lh3+(1*5)+0]=tmp2*fjac[kj3+(1*5)+0]
					-tmp1*njac[kj3+(1*5)+0];
				lhs[lh4+lh3+(2*5)+0]=tmp2*fjac[kj3+(2*5)+0]
					-tmp1*njac[kj3+(2*5)+0];
				lhs[lh4+lh3+(3*5)+0]=tmp2*fjac[kj3+(3*5)+0]
					-tmp1*njac[kj3+(3*5)+0];
				lhs[lh4+lh3+(4*5)+0]=tmp2*fjac[kj3+(4*5)+0]
					-tmp1*njac[kj3+(4*5)+0];
				lhs[lh4+lh3+(0*5)+1]=tmp2*fjac[kj3+(0*5)+1]
					-tmp1*njac[kj3+(0*5)+1];
				lhs[lh4+lh3+(1*5)+1]=tmp2*fjac[kj3+(1*5)+1]
					-tmp1*njac[kj3+(1*5)+1]
					-tmp1*dy2;
				lhs[lh4+lh3+(2*5)+1]=tmp2*fjac[kj3+(2*5)+1]
					-tmp1*njac[kj3+(2*5)+1];
				lhs[lh4+lh3+(3*5)+1]=tmp2*fjac[kj3+(3*5)+1]
					-tmp1*njac[kj3+(3*5)+1];
				lhs[lh4+lh3+(4*5)+1]=tmp2*fjac[kj3+(4*5)+1]
					-tmp1*njac[kj3+(4*5)+1];
				lhs[lh4+lh3+(0*5)+2]=tmp2*fjac[kj3+(0*5)+2]
					-tmp1*njac[kj3+(0*5)+2];
				lhs[lh4+lh3+(1*5)+2]=tmp2*fjac[kj3+(1*5)+2]
					-tmp1*njac[kj3+(1*5)+2];
				lhs[lh4+lh3+(2*5)+2]=tmp2*fjac[kj3+(2*5)+2]
					-tmp1*njac[kj3+(2*5)+2]
					-tmp1*dy3;
				lhs[lh4+lh3+(3*5)+2]=tmp2*fjac[kj3+(3*5)+2]
					-tmp1*njac[kj3+(3*5)+2];
				lhs[lh4+lh3+(4*5)+2]=tmp2*fjac[kj3+(4*5)+2]
					-tmp1*njac[kj3+(4*5)+2];
				lhs[lh4+lh3+(0*5)+3]=tmp2*fjac[kj3+(0*5)+3]
					-tmp1*njac[kj3+(0*5)+3];
				lhs[lh4+lh3+(1*5)+3]=tmp2*fjac[kj3+(1*5)+3]
					-tmp1*njac[kj3+(1*5)+3];
				lhs[lh4+lh3+(2*5)+3]=tmp2*fjac[kj3+(2*5)+3]
					-tmp1*njac[kj3+(2*5)+3];
				lhs[lh4+lh3+(3*5)+3]=tmp2*fjac[kj3+(3*5)+3]
					-tmp1*njac[kj3+(3*5)+3]
					-tmp1*dy4;
				lhs[lh4+lh3+(4*5)+3]=tmp2*fjac[kj3+(4*5)+3]
					-tmp1*njac[kj3+(4*5)+3];
				lhs[lh4+lh3+(0*5)+4]=tmp2*fjac[kj3+(0*5)+4]
					-tmp1*njac[kj3+(0*5)+4];
				lhs[lh4+lh3+(1*5)+4]=tmp2*fjac[kj3+(1*5)+4]
					-tmp1*njac[kj3+(1*5)+4];
				lhs[lh4+lh3+(2*5)+4]=tmp2*fjac[kj3+(2*5)+4]
					-tmp1*njac[kj3+(2*5)+4];
				lhs[lh4+lh3+(3*5)+4]=tmp2*fjac[kj3+(3*5)+4]
					-tmp1*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+(4*5)+4]=tmp2*fjac[kj3+(4*5)+4]
					-tmp1*njac[kj3+(4*5)+4]
					-tmp1*dy5;
			}
			/*
			 * ---------------------------------------------------------------------
			 * performs guaussian elimination on this cell.
			 *
			 * assumes that unpacking routines for non-first cells 
			 * preload c' and rhs' from previous cell.
			 * 
			 * assumed send happens outside this routine, but that
			 * c'(JMAX) and rhs'(JMAX) will be sent to next cell
			 * ---------------------------------------------------------------------
			 * multiply c(i,0,k) by b_inverse and copy back to c
			 * multiply rhs(0) by b_inverse(0) and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh3 = BB*lh3_const;
			binvcrhs(&lhs[lh3], &lhs[CC*lh3_const], &rhs[ku4+ku3]);
			/*
			 * ---------------------------------------------------------------------
			 * begin inner most do loop
			 * do all the elements of the cell unless last 
			 * ---------------------------------------------------------------------
			 */
			for(int j=1; j<=jsize-1; j++){
				ku3 = j*ku3_const;
				lh4 = j*lh4_const;
				/*
				 * -------------------------------------------------------------------
				 * subtract A*lhs_vector(j-1) from lhs_vector(j)
				 *  
				 * rhs(j) = rhs(j) - A*rhs(j-1)
				 * -------------------------------------------------------------------
				 */
				matvec_sub(&lhs[lh4 + AA*lh3_const], &rhs[ku4+ku3-ku3_const+ku2], &rhs[ku4+ku3+ku2]);
				/*
				 * -------------------------------------------------------------------
				 * B(j) = B(j) - C(j-1)*A(j)
				 * -------------------------------------------------------------------
				 */
				matmul_sub(&lhs[lh4+lh3AA], &lhs[lh4-lh4_const+lh3CC], &lhs[lh4+lh3BB]);
				/*
				 * -------------------------------------------------------------------
				 * multiply c(i,j,k) by b_inverse and copy back to c
				 * multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
				 * -------------------------------------------------------------------
				 */
				binvcrhs(&lhs[lh4+lh3], &lhs[lh4+CC*lh3_const], &rhs[ku4+ku3+ku2]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
			 * ---------------------------------------------------------------------
			 */
			lh4 = jsize*lh4_const;
			ku3 = (jsize)*ku3_const;
			matvec_sub(&lhs[lh4+lh3AA], &rhs[ku4+ku3-ku3_const+ku2], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
			 * matmul_sub(aa,i,jsize,k,c,
			 * $ cc,i,jsize-1,k,c,bb,i,jsize,k)
			 * ---------------------------------------------------------------------
			 */
			matmul_sub(&lhs[lh4+lh3AA], &lhs[lh4-lh4_const+lh3CC], &lhs[lh4+lh3BB]);
			/*
			 * ---------------------------------------------------------------------
			 * multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh4 = jsize*lh4_const;
			// lh3 = BB*lh3_const;
			binvrhs(&lhs[lh4+lh3BB], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * back solve: if last cell, then generate U(jsize)=rhs(jsize)
			 * else assume U(jsize) is loaded in un pack backsub_info
			 * so just use it
			 * after u(jstart) will be sent to next cell
			 * ---------------------------------------------------------------------
			 */
			for(int j=jsize-1; j>=0; j--){
				lh4 = j*lh4_const;
				ku3 = j*ku3_const;
				for(int m=0; m<BLOCK_SIZE; m++){
					for(int n=0; n<BLOCK_SIZE; n++){
						lh2 = n*5;
						rhs[ku4+ku3+ku2+m]-=lhs[lh4+lh3CC+lh2+m]*rhs[ku4+ku3+ku3_const+ku2+n];
					}
				}
			}
		}
	});
	if(timeron){timer_stop(T_YSOLVE);}
}

/*
 * ---------------------------------------------------------------------
 * performs line solves in Z direction by first factoring
 * the block-tridiagonal matrix into an upper triangular matrix, 
 * and then performing back substitution to solve for the unknow
 * vectors of each line.  
 *  
 * make sure we treat elements zero to cell_size in the direction
 * of the sweep.
 * ---------------------------------------------------------------------
 */
void z_solve(){
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	int ks3, ks2;
	int ks3_const=(JMAXP+1)*(IMAXP+1), ks2_const=(IMAXP+1);

	int lh4, lh3, lh2;
	int lh4_const=75, lh3_const=25;
	int lh3AA=lh3_const*AA, lh3BB=lh3_const*BB, lh3CC=lh3_const*CC;

	int kj3, kj3_const=25;

	if(timeron){timer_start(T_ZSOLVE);}
	/*
	 * ---------------------------------------------------------------------
	 * this function computes the left hand side for the three z-factors   
	 * ---------------------------------------------------------------------
	 */
	int ksize = grid_points[2]-1;
	/*
	 * ---------------------------------------------------------------------
	 * compute the indices for storing the block-diagonal matrix;
	 * determine c (labeled f) and s jacobians
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front()+1, grid_points[1]-2, [&](int j){
		int ku4, ku3, ku2, ks3, ks2, kj3, kj2;
		int lh4, lh3, lh2;
		std::vector<double> fjac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double> njac((PROBLEM_SIZE+1)*5*5, 0.0);
		std::vector<double>  lhs((PROBLEM_SIZE+1)*3*5*5, 0.0);
    	double tmp1, tmp2, tmp3;
		
		ku3=j*ku3_const;
		ks2=j*ks2_const;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			for(int k=0; k<=ksize; k++){
				kj3 = k*kj3_const;
				ks3=k*ks3_const;
				ku4 = k*ku4_const;
				tmp1=1.0/u[ku4+ku3+ku2];
				tmp2=tmp1*tmp1;
				tmp3=tmp1*tmp2;
				fjac[kj3]=0.0;
				fjac[kj3+5]=0.0;
				fjac[kj3+10]=0.0;
				fjac[kj3+15]=1.0;
				fjac[kj3+20]=0.0;
				fjac[kj3+1]=-(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+5+1]=u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+10+1]=0.0;
				fjac[kj3+15+1]=u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+20+1]=0.0;
				fjac[kj3+2]=-(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+5+2]=0.0;
				fjac[kj3+10+2]=u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+15+2]=u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+20+2]=0.0;
				fjac[kj3+3]=-(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3]*tmp2)+c2*qs[ks3+ks2+i];
				fjac[kj3+5+3]=-c2*u[ku4+ku3+ku2+1]*tmp1;
				fjac[kj3+10+3]=-c2*u[ku4+ku3+ku2+2]*tmp1;
				fjac[kj3+15+3]=(2.0-c2)*u[ku4+ku3+ku2+3]*tmp1;
				fjac[kj3+20+3]=c2;
				fjac[kj3+4]=(c2*2.0*square[ks3+ks2+i]-c1*u[ku4+ku3+ku2+4])*u[ku4+ku3+ku2+3]*tmp2;
				fjac[kj3+5+4]=-c2*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+10+4]=-c2*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+3])*tmp2;
				fjac[kj3+15+4]=c1*(u[ku4+ku3+ku2+4]*tmp1)-c2*(qs[ks3+ks2+i]+u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3]*tmp2);
				fjac[kj3+20+4]=c1*u[ku4+ku3+ku2+3]*tmp1;
				njac[kj3]=0.0;
				njac[kj3+5]=0.0;
				njac[kj3+10]=0.0;
				njac[kj3+15]=0.0;
				njac[kj3+20]=0.0;
				njac[kj3+1]=-c3c4*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+5+1]=c3c4*tmp1;
				njac[kj3+10+1]=0.0;
				njac[kj3+15+1]=0.0;
				njac[kj3+20+1]=0.0;
				njac[kj3+2]=-c3c4*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+5+2]=0.0;
				njac[kj3+10+2]=c3c4*tmp1;
				njac[kj3+15+2]=0.0;
				njac[kj3+20+2]=0.0;
				njac[kj3+3]=-con43*c3c4*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+5+3]=0.0;
				njac[kj3+10+3]=0.0;
				njac[kj3+15+3]=con43*c3*c4*tmp1;
				njac[kj3+20+3]=0.0;
				njac[kj3+4]=-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1])
					-(c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2])
					-(con43*c3c4-c1345)*tmp3*(u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])
					-c1345*tmp2*u[ku4+ku3+ku2+4];
				njac[kj3+5+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+1];
				njac[kj3+10+4]=(c3c4-c1345)*tmp2*u[ku4+ku3+ku2+2];
				njac[kj3+15+4]=(con43*c3c4-c1345)*tmp2*u[ku4+ku3+ku2+3];
				njac[kj3+20+4]=(c1345)*tmp1;
			}
			/*
			 * ---------------------------------------------------------------------
			 * now jacobians set, so form left hand side in z direction
			 * ---------------------------------------------------------------------
			 */
			lhsinit(lhs, ksize);
			for(int k=1; k<=ksize-1; k++){
				lh4 = k*lh4_const;
				lh3 = AA*lh3_const;
				kj3 = (k-1)*kj3_const;
				tmp1=dt*tz1;
				tmp2=dt*tz2;
				lhs[lh4+lh3]=-tmp2*fjac[kj3]
					-tmp1*njac[kj3]
					-tmp1*dz1; 
				lhs[lh4+lh3+5]=-tmp2*fjac[kj3+5]
					-tmp1*njac[kj3+5];
				lhs[lh4+lh3+10]=-tmp2*fjac[kj3+10]
					-tmp1*njac[kj3+10];
				lhs[lh4+lh3+15]=-tmp2*fjac[kj3+15]
					-tmp1*njac[kj3+15];
				lhs[lh4+lh3+20]=-tmp2*fjac[kj3+20]
					-tmp1*njac[kj3+20];
				lhs[lh4+lh3+1]=-tmp2*fjac[kj3+1]
					-tmp1*njac[kj3+1];
				lhs[lh4+lh3+5+1]=-tmp2*fjac[kj3+5+1]
					-tmp1*njac[kj3+5+1]
					-tmp1*dz2;
				lhs[lh4+lh3+10+1]=-tmp2*fjac[kj3+10+1]
					-tmp1*njac[kj3+10+1];
				lhs[lh4+lh3+15+1]=-tmp2*fjac[kj3+15+1]
					-tmp1*njac[kj3+15+1];
				lhs[lh4+lh3+20+1]=-tmp2*fjac[kj3+20+1]
					-tmp1*njac[kj3+20+1];
				lhs[lh4+lh3+2]=-tmp2*fjac[kj3+2]
					-tmp1*njac[kj3+2];
				lhs[lh4+lh3+5+2]=-tmp2*fjac[kj3+5+2]
					-tmp1*njac[kj3+5+2];
				lhs[lh4+lh3+10+2]=-tmp2*fjac[kj3+10+2]
					-tmp1*njac[kj3+10+2]
					-tmp1*dz3;
				lhs[lh4+lh3+15+2]=-tmp2*fjac[kj3+15+2]
					-tmp1*njac[kj3+15+2];
				lhs[lh4+lh3+20+2]=-tmp2*fjac[kj3+20+2]
					-tmp1*njac[kj3+20+2];
				lhs[lh4+lh3+3]=-tmp2*fjac[kj3+3]
					-tmp1*njac[kj3+3];
				lhs[lh4+lh3+5+3]=-tmp2*fjac[kj3+5+3]
					-tmp1*njac[kj3+5+3];
				lhs[lh4+lh3+10+3]=-tmp2*fjac[kj3+10+3]
					-tmp1*njac[kj3+10+3];
				lhs[lh4+lh3+15+3]=-tmp2*fjac[kj3+15+3]
					-tmp1*njac[kj3+15+3]
					-tmp1*dz4;
				lhs[lh4+lh3+20+3]=-tmp2*fjac[kj3+20+3]
					-tmp1*njac[kj3+20+3];
				lhs[lh4+lh3+4]=-tmp2*fjac[kj3+4]
					-tmp1*njac[kj3+4];
				lhs[lh4+lh3+5+4]=-tmp2*fjac[kj3+5+4]
					-tmp1*njac[kj3+5+4];
				lhs[lh4+lh3+10+4]=-tmp2*fjac[kj3+10+4]
					-tmp1*njac[kj3+10+4];
				lhs[lh4+lh3+15+4]=-tmp2*fjac[kj3+15+4]
					-tmp1*njac[kj3+15+4];
				lhs[lh4+lh3+20+4]=-tmp2*fjac[kj3+20+4]
					-tmp1*njac[kj3+20+4]
					-tmp1*dz5;
				lh3 = BB*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3]=1.0
					+tmp1*2.0*njac[kj3]
					+tmp1*2.0*dz1;
				lhs[lh4+lh3+5]=tmp1*2.0*njac[kj3+5];
				lhs[lh4+lh3+10]=tmp1*2.0*njac[kj3+10];
				lhs[lh4+lh3+15]=tmp1*2.0*njac[kj3+15];
				lhs[lh4+lh3+20]=tmp1*2.0*njac[kj3+20];
				lhs[lh4+lh3+1]=tmp1*2.0*njac[kj3+1];
				lhs[lh4+lh3+5+1]=1.0
					+tmp1*2.0*njac[kj3+5+1]
					+tmp1*2.0*dz2;
				lhs[lh4+lh3+10+1]=tmp1*2.0*njac[kj3+10+1];
				lhs[lh4+lh3+15+1]=tmp1*2.0*njac[kj3+15+1];
				lhs[lh4+lh3+20+1]=tmp1*2.0*njac[kj3+20+1];
				lhs[lh4+lh3+2]=tmp1*2.0*njac[kj3+2];
				lhs[lh4+lh3+5+2]=tmp1*2.0*njac[kj3+5+2];
				lhs[lh4+lh3+10+2]=1.0
					+tmp1*2.0*njac[kj3+10+2]
					+tmp1*2.0*dz3;
				lhs[lh4+lh3+15+2]=tmp1*2.0*njac[kj3+15+2];
				lhs[lh4+lh3+20+2]=tmp1*2.0*njac[kj3+20+2];
				lhs[lh4+lh3+3]=tmp1*2.0*njac[kj3+3];
				lhs[lh4+lh3+5+3]=tmp1*2.0*njac[kj3+5+3];
				lhs[lh4+lh3+10+3]=tmp1*2.0*njac[kj3+10+3];
				lhs[lh4+lh3+15+3]=1.0
					+tmp1*2.0*njac[kj3+15+3]
					+tmp1*2.0*dz4;
				lhs[lh4+lh3+20+3]=tmp1*2.0*njac[kj3+20+3];
				lhs[lh4+lh3+4]=tmp1*2.0*njac[kj3+4];
				lhs[lh4+lh3+5+4]=tmp1*2.0*njac[kj3+5+4];
				lhs[lh4+lh3+10+4]=tmp1*2.0*njac[kj3+10+4];
				lhs[lh4+lh3+15+4]=tmp1*2.0*njac[kj3+15+4];
				lhs[lh4+lh3+20+4]=1.0
					+tmp1*2.0*njac[kj3+20+4] 
					+tmp1*2.0*dz5;
				lh3 = CC*lh3_const;
				kj3 += kj3_const;
				lhs[lh4+lh3]=tmp2*fjac[kj3]
					-tmp1*njac[kj3]
					-tmp1*dz1;
				lhs[lh4+lh3+5]=tmp2*fjac[kj3+5]
					-tmp1*njac[kj3+5];
				lhs[lh4+lh3+10]=tmp2*fjac[kj3+10]
					-tmp1*njac[kj3+10];
				lhs[lh4+lh3+15]=tmp2*fjac[kj3+15]
					-tmp1*njac[kj3+15];
				lhs[lh4+lh3+20]=tmp2*fjac[kj3+20]
					-tmp1*njac[kj3+20];
				lhs[lh4+lh3+1]=tmp2*fjac[kj3+1]
					-tmp1*njac[kj3+1];
				lhs[lh4+lh3+5+1]=tmp2*fjac[kj3+5+1]
					-tmp1*njac[kj3+5+1]
					-tmp1*dz2;
				lhs[lh4+lh3+10+1]=tmp2*fjac[kj3+10+1]
					-tmp1*njac[kj3+10+1];
				lhs[lh4+lh3+15+1]=tmp2*fjac[kj3+15+1]
					-tmp1*njac[kj3+15+1];
				lhs[lh4+lh3+20+1]=tmp2*fjac[kj3+20+1]
					-tmp1*njac[kj3+20+1];
				lhs[lh4+lh3+2]=tmp2*fjac[kj3+2]
					-tmp1*njac[kj3+2];
				lhs[lh4+lh3+5+2]= tmp2*fjac[kj3+5+2]
					-tmp1*njac[kj3+5+2];
				lhs[lh4+lh3+10+2]=tmp2*fjac[kj3+10+2]
					-tmp1*njac[kj3+10+2]
					-tmp1*dz3;
				lhs[lh4+lh3+15+2]=tmp2*fjac[kj3+15+2]
					-tmp1*njac[kj3+15+2];
				lhs[lh4+lh3+20+2]=tmp2*fjac[kj3+20+2]
					-tmp1*njac[kj3+20+2];					
				lhs[lh4+lh3+3]=tmp2*fjac[kj3+3]
					-tmp1*njac[kj3+3];
				lhs[lh4+lh3+5+3]=tmp2*fjac[kj3+5+3]
					-tmp1*njac[kj3+5+3];
				lhs[lh4+lh3+10+3]=tmp2*fjac[kj3+10+3]
					-tmp1*njac[kj3+10+3];
				lhs[lh4+lh3+15+3]=tmp2*fjac[kj3+15+3]
					-tmp1*njac[kj3+15+3]
					-tmp1*dz4;
				lhs[lh4+lh3+20+3]=tmp2*fjac[kj3+20+3]
					-tmp1*njac[kj3+20+3];
				lhs[lh4+lh3+4]=tmp2*fjac[kj3+4]
					-tmp1*njac[kj3+4];
				lhs[lh4+lh3+5+4]=tmp2*fjac[kj3+5+4]
					-tmp1*njac[kj3+5+4];
				lhs[lh4+lh3+10+4]=tmp2*fjac[kj3+10+4]
					-tmp1*njac[kj3+10+4];
				lhs[lh4+lh3+(3*5)+4]=tmp2*fjac[kj3+(3*5)+4]
					-tmp1*njac[kj3+(3*5)+4];
				lhs[lh4+lh3+20+4]=tmp2*fjac[kj3+20+4]
					-tmp1*njac[kj3+20+4]
					-tmp1*dz5;
			}
			/*
			 * ---------------------------------------------------------------------
			 * performs guaussian elimination on this cell.
			 *  
			 * assumes that unpacking routines for non-first cells 
			 * preload c' and rhs' from previous cell.
			 *  
			 * assumed send happens outside this routine, but that
			 * c'(KMAX) and rhs'(KMAX) will be sent to next cell.
			 * ---------------------------------------------------------------------
			 * outer most do loops - sweeping in i direction
			 * ---------------------------------------------------------------------
			 * multiply c(i,j,0) by b_inverse and copy back to c
			 * multiply rhs(0) by b_inverse(0) and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh3 = BB*lh3_const;
			binvcrhs(&lhs[lh3], &lhs[CC*lh3_const], &rhs[ku4+ku3]);
			/*
			 * ---------------------------------------------------------------------
			 * begin inner most do loop
			 * do all the elements of the cell unless last 
			 * ---------------------------------------------------------------------
			 */
			for(int k=1; k<=ksize-1; k++){
				ku4 = k*ku4_const;
				lh4 = k*lh4_const;
				/*
				 * -------------------------------------------------------------------
				 * subtract A*lhs_vector(k-1) from lhs_vector(k)
				 *  
				 * rhs(k) = rhs(k) - A*rhs(k-1)
				 * -------------------------------------------------------------------
				 */
				matvec_sub(&lhs[lh4+lh3AA], &rhs[ku4-ku4_const+ku3+ku2], &rhs[ku4+ku3+ku2]);
				/*
				 * -------------------------------------------------------------------
				 * B(k) = B(k) - C(k-1)*A(k)
				 * matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
				 * --------------------------------------------------------------------
				 */
				matmul_sub(&lhs[lh4+lh3AA], &lhs[lh4-lh4_const+lh3CC], &lhs[lh4+lh3BB]);
				/*
				 * -------------------------------------------------------------------
				 * multiply c(i,j,k) by b_inverse and copy back to c
				 * multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
				 * -------------------------------------------------------------------
				 */
				binvcrhs(&lhs[lh4+lh3], &lhs[lh4+lh3CC], &rhs[ku4+ku3+ku2]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * now finish up special cases for last cell
			 * ---------------------------------------------------------------------
			 * rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
			 * ---------------------------------------------------------------------
			 */
			ku4 = ksize*ku4_const;
			lh4 = ksize*lh4_const;
			matvec_sub(&lhs[lh4+lh3AA], &rhs[ku4-ku4_const+ku3+ku2], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
			 * matmul_sub(aa,i,j,ksize,c,
			 * $ cc,i,j,ksize-1,c,bb,i,j,ksize)
			 * ---------------------------------------------------------------------
			 */
			matmul_sub(&lhs[lh4+lh3AA], &lhs[lh4-lh4_const+lh3CC], &lhs[lh4+lh3BB]);
			/*
			 * ---------------------------------------------------------------------
			 * multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
			 * ---------------------------------------------------------------------
			 */
			lh4 = ksize*lh4_const;
			lh3 = BB*lh3_const;
			binvrhs(&lhs[lh4+lh3], &rhs[ku4+ku3+ku2]);
			/*
			 * ---------------------------------------------------------------------
			 * back solve: if last cell, then generate U(ksize)=rhs(ksize)
			 * else assume U(ksize) is loaded in un pack backsub_info
			 * so just use it
			 * after u(kstart) will be sent to next cell
			 * ---------------------------------------------------------------------
			 */
			for(int k=ksize-1; k>=0; k--){
				lh4 = k*lh4_const;
				ku4 = k*ku4_const;
				for(int m=0; m<BLOCK_SIZE; m++){
					for(int n=0; n<BLOCK_SIZE; n++){
						lh2 = n*5;
						rhs[ku4+ku3+ku2+m]-=lhs[lh4+lh3CC+lh2+m]*rhs[ku4+ku4_const+ku3+ku2+n];
					}
				}
			}
		}
	});
	if(timeron){timer_stop(T_ZSOLVE);}
}
