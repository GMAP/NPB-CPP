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
	W. Saphir 
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
#define T_TXINVR 11
#define T_PINVR 12
#define T_NINVR 13
#define T_TZETAR 14
#define T_ADD 15
#define T_LAST 15

/* global variables */
std::vector<double> u((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> us((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> vs((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> ws((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> qs((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> rho_i((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> speed((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> square((JMAXP+1)*(IMAXP+1)*KMAX);
std::vector<double> rhs((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> forcing((KMAX)*(JMAXP+1)*(IMAXP+1)*(5));
std::vector<double> cv(PROBLEM_SIZE);
std::vector<double> rhon(PROBLEM_SIZE);
std::vector<double> rhos(PROBLEM_SIZE);
std::vector<double> rhoq(PROBLEM_SIZE);
std::vector<double> cuf(PROBLEM_SIZE);
std::vector<double> q(PROBLEM_SIZE);
std::vector<double> ue(PROBLEM_SIZE*5);
std::vector<double> buf(PROBLEM_SIZE*5);
std::vector<double> lhs((IMAXP+1)*(IMAXP+1)*5, 0.0);
std::vector<double> lhsp((IMAXP+1)*(IMAXP+1)*5, 0.0);
std::vector<double> lhsm((IMAXP+1)*(IMAXP+1)*5, 0.0);
std::vector<double> ce(13*5);
static double tx1, tx2, tx3,ty1, ty2, ty3, tz1, tz2, tz3, 
	      dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, 
	      dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, 
	      dxmax, dymax, dzmax, xxcon1, xxcon2, 
	      xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
	      dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
	      yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
	      zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, 
	      dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, 
	      dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, 
	      c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1, bt,
	      dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, 
	      c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, 
	      c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;
static int nx2, ny2, nz2;
std::vector<int> grid_points(3, 0);
static boolean timeron;

CountIterator iter(0);
static auto policy = std::execution::par;

/* function prototypes */
void add();
void adi();
void compute_rhs();
void error_norm(std::vector<double> &rms);
void exact_rhs();
void exact_solution(double xi, double eta, double zeta, std::span<double> dtemp);
void initialize();
void lhsinit(int ni, int nj);
void lhsinitj(int nj, int ni);
void ninvr();
void pinvr();
void rhs_norm(std::vector<double> &rms);
void set_constants();
void txinvr();
void tzetar();
void verify(int no_time_steps, char* class_npb, boolean* verified);
void x_solve();
void y_solve();
void z_solve();

/* sp */
int main(int argc, char* argv[]){
	int niter, n3;
	double mflops, t, tmax;
  std::vector<double> trecs(T_LAST+1, 0);
	boolean verified;
	char class_npb;
	std::vector<std::string> t_names(T_LAST+1);
	/*
	 * ---------------------------------------------------------------------
	 * read input file (if it exists), else take
	 * defaults from parameters
	 * ---------------------------------------------------------------------
	 */
	std::fstream fp;
  fp.open("inputsp.data", std::ios::in);
	if(fp.is_open()){
		std::cout
      << " Reading from input file inputsp.data"
      << std::endl;
    fp >> niter;
		while(fp.get() != '\n');
		fp >> dt;
		while(fp.get() != '\n');
    fp >> grid_points[0] >> grid_points[1] >> grid_points[2];
    fp.close();
	}else{
		std::cout
      << " No input file inputsp.data. Using compiled defaults"
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
		t_names[T_TZETAR]="tzetar";
		t_names[T_NINVR]="ninvr";
		t_names[T_PINVR]="pinvr";
		t_names[T_TXINVR]="txinvr";
		t_names[T_ADD]="add";
		fp.close();
	}else{
		timeron = FALSE;
	}

	iter = CountIterator(std::max(std::max(niter+1, grid_points[0]), std::max(grid_points[1], grid_points[2])));

  std::cout
    << std::endl
    << std::endl;
  
  std::cout
    << " NAS Parallel Benchmarks 4.1 Serial C++ version - SP Benchmark";
  
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
	nx2=grid_points[0]-2;
	ny2=grid_points[1]-2;
	nz2=grid_points[2]-2;
	set_constants();
	
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	
	exact_rhs();
	initialize();
	/*
	 * ---------------------------------------------------------------------
	 * do one time step to touch all code, and reinitialize
	 * ---------------------------------------------------------------------
	 */
	adi();
	initialize();
	for(int i=1;i<=T_LAST;i++){timer_clear(i);}
	timer_start(1);
	
	for(int step=1;step<=niter;step++){
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
	if(tmax!=0.0){
		n3=grid_points[0]*grid_points[1]*grid_points[2];
		t=(grid_points[0]+grid_points[1]+grid_points[2])/3.0;
		mflops=(881.174*(double)n3-
				4683.91*(t*t)+
				11484.5*t-
				19272.4)*(double)niter/(tmax*1000000.0);
	}else{
		mflops=0.0;
	}
	c_print_results("SP",
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
				<< "    --> "
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
				<< "    --> "
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
				<< "    --> "
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
	if(timeron){timer_start(T_ADD);}
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		int ku4, ku3, ku2;
    	ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					u[ku4+ku3+ku2+m]+=rhs[ku4+ku3+ku2+m];
				}
			}
		}
	});
	if(timeron){timer_stop(T_ADD);}
}

void adi(){
	compute_rhs();
	txinvr();
	x_solve();
	y_solve();
	z_solve();
	add();
}

void compute_rhs(){
	int i, j, k;
	double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;
	int ku4P, ku4M, ku4, ku3, ku2, ks3, ks2;
	int ku4_const = (JMAXP+1)*(IMAXP+1)*5;
	int ku3_const = (IMAXP+1)*5;
	int ks3_const = (JMAXP+1)*(IMAXP+1);
	int ks2_const = (IMAXP+1);
	if(timeron){timer_start(T_RHS);}
	/*
	 * ---------------------------------------------------------------------
	 * compute the reciprocal of density, and the kinetic energy, 
	 * and the speed of sound. 
	 * ---------------------------------------------------------------------
	 */
	std::for_each_n(policy, iter.front(), grid_points[2], [&](int k){
		int ku4, ku3, ku2, ks3, ks2;
		double rho_inv, aux;
		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=0; j<=grid_points[1]-1; j++){
			ks2 = j*ks2_const;
			ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-1; i++){
        		ku2 = i*5;
				rho_inv=1.0/u[ku4+ku3+ku2];
				rho_i[ks3+ks2+i]=rho_inv;
				us[ks3+ks2+i]=u[ku4+ku3+ku2+1]*rho_inv;
				vs[ks3+ks2+i]=u[ku4+ku3+ku2+2]*rho_inv;
				ws[ks3+ks2+i]=u[ku4+ku3+ku2+3]*rho_inv;
				square[ks3+ks2+i]=0.5*(
						u[ku4+ku3+ku2+1]*u[ku4+ku3+ku2+1]+ 
						u[ku4+ku3+ku2+2]*u[ku4+ku3+ku2+2]+
						u[ku4+ku3+ku2+3]*u[ku4+ku3+ku2+3])*rho_inv;
				qs[ks3+ks2+i]=square[ks3+ks2+i]*rho_inv;
				/*
				 * ---------------------------------------------------------------------
				 * (don't need speed and ainx until the lhs computation)
				 * ---------------------------------------------------------------------
				 */
				aux=c1c2*rho_inv*(u[ku4+ku3+ku2+4]-square[ks3+ks2+i]);
				speed[ks3+ks2+i]=sqrt(aux);
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * copy the exact forcing term to the right hand side;  because 
	 * this forcing term is known, we can store it on the whole grid
	 * including the boundary                   
	 * ---------------------------------------------------------------------
	 */
	std::copy(policy, forcing.begin(), forcing.end(), rhs.begin());

	/*
	 * ---------------------------------------------------------------------
	 * compute xi-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron){timer_start(T_RHSX);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		int ku4, ku3, ku2, ks3, ks2;
		double up1, um1, uijk;
		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
			ks2 = j*ks2_const;
			ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				uijk=us[ks3+ks2+i];
				up1=us[ks3+ks2+i+1];
				um1=us[ks3+ks2+i-1];
				rhs[ku4+ku3+ku2]+=dx1tx1* 
					(u[ku4+ku3+ku2+5]-2.0*u[ku4+ku3+ku2]+u[ku4+ku3+ku2-5])-
					tx2*(u[ku4+ku3+ku2+5+1]-u[ku4+ku3+ku2-5+1]);
				rhs[ku4+ku3+ku2+1]+=dx2tx1* 
					(u[ku4+ku3+ku2+5+1]-2.0*u[ku4+ku3+ku2+1]+u[ku4+ku3+ku2-5+1])+
					xxcon2*con43*(up1-2.0*uijk+um1)-
					tx2*(u[ku4+ku3+ku2+5+1]*up1-u[ku4+ku3+ku2-5+1]*um1+
							(u[ku4+ku3+ku2+5+4]-square[ks3+ks2+i+1]-
							 u[ku4+ku3+ku2-5+4]+square[ks3+ks2+i-1])*c2);
				rhs[ku4+ku3+ku2+2]+=dx3tx1* 
					(u[ku4+ku3+ku2+5+2]-2.0*u[ku4+ku3+ku2+2]+u[ku4+ku3+ku2-5+2])+
					xxcon2*(vs[ks3+ks2+i+1]-2.0*vs[ks3+ks2+i]+vs[ks3+ks2+i-1])-
					tx2*(u[ku4+ku3+ku2+5+2]*up1-u[ku4+ku3+ku2-5+2]*um1);
				rhs[ku4+ku3+ku2+3]+=dx4tx1* 
					(u[ku4+ku3+ku2+5+3]-2.0*u[ku4+ku3+ku2+3]+u[ku4+ku3+ku2-5+3])+
					xxcon2*(ws[ks3+ks2+i+1]-2.0*ws[ks3+ks2+i]+ws[ks3+ks2+i-1])-
					tx2*(u[ku4+ku3+ku2+5+3]*up1-u[ku4+ku3+ku2-5+3]*um1);
				rhs[ku4+ku3+ku2+4]+=dx5tx1* 
					(u[ku4+ku3+ku2+5+4]-2.0*u[ku4+ku3+ku2+4]+u[ku4+ku3+ku2-5+4])+
					xxcon3*(qs[ks3+ks2+i+1]-2.0*qs[ks3+ks2+i]+qs[ks3+ks2+i-1])+
					xxcon4*(up1*up1-2.0*uijk*uijk+um1*um1)+
					xxcon5*(u[ku4+ku3+ku2+5+4]*rho_i[ks3+ks2+i+1]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4+ku3+ku2-5+4]*rho_i[ks3+ks2+i-1])-
					tx2*((c1*u[ku4+ku3+ku2+5+4]-c2*square[ks3+ks2+i+1])*up1-
							(c1*u[ku4+ku3+ku2-5+4]-c2*square[ks3+ks2+i-1])*um1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order xi-direction dissipation               
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			i=1;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+5+m]-=dssp* 
					(5.0*u[ku4+ku3+5+m]-4.0*u[ku4+ku3+10+m]+u[ku4+ku3+15+m]);
			}
			i=2;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+10+m]-=dssp* 
					(-4.0*u[ku4+ku3+5+m]+6.0*u[ku4+ku3+10+m]-
					 4.0*u[ku4+ku3+15+m]+u[ku4+ku3+20+m]);
			}
		}

		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=3; i<=nx2-2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]-=dssp* 
						(u[ku4+ku3+ku2-10+m]-4.0*u[ku4+ku3+ku2-5+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3+ku2+5+m]+ 
						 u[ku4+ku3+ku2+10+m]);
				}
			}
		}

		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
      		ku2 = (nx2-1)*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m] -= dssp*
					(u[ku4+ku3+ku2-10+m]-4.0*u[ku4+ku3+ku2-5+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3+ku2+5+m]);
			}
      		ku2 = nx2*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m] -= dssp*
					(u[ku4+ku3+ku2-10+m]-4.0*u[ku4+ku3+ku2-5+m]+5.0*u[ku4+ku3+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSX);}
	/*
	 * ---------------------------------------------------------------------
	 * compute eta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron){timer_start(T_RHSY);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		int ku3P, ku3M;
		int ku4, ku3, ku2, ks3, ks2;
		double vp1, vm1, vijk;

		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
			ks2 = j*ks2_const;
			ku3 = j*ku3_const;
			ku3P = ku3+ku3_const;
			ku3M = ku3-ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				vijk=vs[ks3+ks2+i];
				vp1=vs[ks3+ks2+ks2_const+i];
				vm1=vs[ks3+ks2-ks2_const+i];
				rhs[ku4+ku3+ku2]+=dy1ty1* 
					(u[ku4+ku3P+ku2]-2.0*u[ku4+ku3+ku2]+u[ku4+ku3M+ku2])-
					ty2*(u[ku4+ku3P+ku2+2]-u[ku4+ku3M+ku2+2]);
				rhs[ku4+ku3+ku2+1]+=dy2ty1* 
					(u[ku4+ku3P+ku2+1]-2.0*u[ku4+ku3+ku2+1]+u[ku4+ku3M+ku2+1])+
					yycon2*(us[ks3+ks2+ks2_const+i]-2.0*us[ks3+ks2+i]+us[ks3+ks2-ks2_const+i])-
					ty2*(u[ku4+ku3P+ku2+1]*vp1-u[ku4+ku3M+ku2+1]*vm1);
				rhs[ku4+ku3+ku2+2]+=dy3ty1* 
					(u[ku4+ku3P+ku2+2]-2.0*u[ku4+ku3+ku2+2]+u[ku4+ku3M+ku2+2])+
					yycon2*con43*(vp1-2.0*vijk+vm1)-
					ty2*(u[ku4+ku3P+ku2+2]*vp1-u[ku4+ku3M+ku2+2]*vm1+
							(u[ku4+ku3P+ku2+4]-square[ks3+ks2+ks2_const+i]- 
							 u[ku4+ku3M+ku2+4]+square[ks3+ks2-ks2_const+i])* c2);
				rhs[ku4+ku3+ku2+3]+=dy4ty1* 
					(u[ku4+ku3P+ku2+3]-2.0*u[ku4+ku3+ku2+3]+u[ku4+ku3M+ku2+3])+
					yycon2*(ws[ks3+ks2+ks2_const+i]-2.0*ws[ks3+ks2+i]+ws[ks3+ks2-ks2_const+i])-
					ty2*(u[ku4+ku3P+ku2+3]*vp1-u[ku4+ku3M+ku2+3]*vm1);
				rhs[ku4+ku3+ku2+4]+=dy5ty1* 
					(u[ku4+ku3P+ku2+4]-2.0*u[ku4+ku3+ku2+4]+u[ku4+ku3M+ku2+4])+
					yycon3*(qs[ks3+ks2+ks2_const+i]-2.0*qs[ks3+ks2+i]+qs[ks3+ks2-ks2_const+i])+
					yycon4*(vp1*vp1- 2.0*vijk*vijk + vm1*vm1)+
					yycon5*(u[ku4+ku3P+ku2+4]*rho_i[ks3+ks2+ks2_const+i]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4+ku3M+ku2+4]*rho_i[ks3+ks2-ks2_const+i])-
					ty2*((c1*u[ku4+ku3P+ku2+4]-c2*square[ks3+ks2+ks2_const+i])*vp1 -
							(c1*u[ku4+ku3M+ku2+4]-c2*square[ks3+ks2-ks2_const+i])*vm1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order eta-direction dissipation         
		 * ---------------------------------------------------------------------
		 */
		j = 1;
    	ku3 = 2*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3-ku3_const+ku2+m]-=dssp*
					(5.0*u[ku4+ku3-ku3_const+ku2+m]-4.0*u[ku4+ku3+ku2+m]+u[ku4+ku3+ku3_const+ku2+m]);
			}
		}
		j = 2;
    	ku3 = 3*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m = 0; m < 5; m++) {
				rhs[ku4+ku3-ku3_const+ku2+m]-=dssp*
					(-4.0*u[ku4+ku3_const+ku2+m]+6.0*u[ku4+ku3-ku3_const+ku2+m]-
					 4.0*u[ku4+ku3+ku2+m]+u[ku4+ku3+ku3_const+ku2+m]);
			}
		}
		for(int j=3; j<=ny2-2; j++){
			ku3 = j*ku3_const;
			ku3P = ku3+ku3_const;
			ku3M = ku3-ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]-=dssp*
						(u[ku4+ku3M-ku3_const+ku2+m]-4.0*u[ku4+ku3M+ku2+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3P+ku2+m]+ 
						 u[ku4+ku3P+ku3_const+ku2+m]);
				}
			}
		}
		j=ny2-1;
		ku3 = j*ku3_const;
		ku3P = ku3+ku3_const;
		ku3M = ku3-ku3_const; 
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]-=dssp*
					(u[ku4+ku3M-ku3_const+ku2+m]-4.0*u[ku4+ku3M+ku2+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4+ku3P+ku2+m]);
			}
		}
		j=ny2;
		ku3 = j*ku3_const;
		ku3M = ku3-ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]-=dssp*
					(u[ku4+ku3M-ku3_const+ku2+m]
           -4.0*u[ku4+ku3M+ku2+m] +5.0*u[ku4+ku3+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSY);}
	/*
	 * ---------------------------------------------------------------------
	 * compute zeta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron){timer_start(T_RHSZ);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		int ku4P, ku4M;
		int ku4, ku3, ku2, ks3, ks2;
		double wp1, wm1, wijk;

		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		ku4P = ku4+ku4_const;
		ku4M = ku4-ku4_const;
		for(int j=1; j<=ny2; j++){
			ks2 = j*ks2_const;
			ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				wijk=ws[ks3+ks2+i];
				wp1=ws[ks3+ks3_const+ks2+i];
				wm1=ws[ks3-ks3_const+ks2+i];
				rhs[ku4+ku3+ku2]+=dz1tz1* 
					(u[ku4P+ku3+ku2]-2.0*u[ku4+ku3+ku2]+u[ku4M+ku3+ku2])-
					tz2*(u[ku4P+ku3+ku2+3]-u[ku4M+ku3+ku2+3]);
				rhs[ku4+ku3+ku2+1]+=dz2tz1* 
					(u[ku4P+ku3+ku2+1]-2.0*u[ku4+ku3+ku2+1]+u[ku4M+ku3+ku2+1])+
					zzcon2*(us[ks3+ks3_const+ks2+i]-2.0*us[ks3+ks2+i]+us[ks3-ks3_const+ks2+i])-
					tz2*(u[ku4P+ku3+ku2+1]*wp1-u[ku4M+ku3+ku2+1]*wm1);
				rhs[ku4+ku3+ku2+2]+=dz3tz1* 
					(u[ku4P+ku3+ku2+2]-2.0*u[ku4+ku3+ku2+2]+u[ku4M+ku3+ku2+2])+
					zzcon2*(vs[ks3+ks3_const+ks2+i]-2.0*vs[ks3+ks2+i]+vs[ks3-ks3_const+ks2+i])-
					tz2*(u[ku4P+ku3+ku2+2]*wp1-u[ku4M+ku3+ku2+2]*wm1);
				rhs[ku4+ku3+ku2+3]+=dz4tz1* 
					(u[ku4P+ku3+ku2+3]-2.0*u[ku4+ku3+ku2+3]+u[ku4M+ku3+ku2+3])+
					zzcon2*con43*(wp1-2.0*wijk+wm1)-
					tz2*(u[ku4P+ku3+ku2+3]*wp1-u[ku4M+ku3+ku2+3]*wm1+
							(u[ku4P+ku3+ku2+4]-square[ks3+ks3_const+ks2+i]- 
							 u[ku4M+ku3+ku2+4]+square[ks3-ks3_const+ks2+i])*c2);
				rhs[ku4+ku3+ku2+4]+=dz5tz1* 
					(u[ku4P+ku3+ku2+4]-2.0*u[ku4+ku3+ku2+4]+u[ku4M+ku3+ku2+4])+
					zzcon3*(qs[ks3+ks3_const+ks2+i]-2.0*qs[ks3+ks2+i]+qs[ks3-ks3_const+ks2+i])+
					zzcon4*(wp1*wp1-2.0*wijk*wijk+wm1*wm1)+
					zzcon5*(u[ku4P+ku3+ku2+4]*rho_i[ks3+ks3_const+ks2+i]- 
							2.0*u[ku4+ku3+ku2+4]*rho_i[ks3+ks2+i]+
							u[ku4M+ku3+ku2+4]*rho_i[ks3-ks3_const+ks2+i])-
					tz2*((c1*u[ku4P+ku3+ku2+4]-c2*square[ks3+ks3_const+ks2+i])*wp1-
							(c1*u[ku4M+ku3+ku2+4]-c2*square[ks3-ks3_const+ks2+i])*wm1);
			}
		}
	});
	/*
	 * ---------------------------------------------------------------------
	 * add fourth order zeta-direction dissipation                
	 * ---------------------------------------------------------------------
	 */
	k=1;
	ku4 = 2*ku4_const;
	std::for_each_n(policy, iter.front()+1, ny2, [&](int j){
		int ku3, ku2;
    	ku3 = j*ku3_const;
		for(int i=1; i<=nx2; i++){
			ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4_const+ku3+ku2+m]-=dssp* 
					(5.0*u[ku4_const+ku3+ku2+m]
					-4.0*u[ku4+ku3+ku2+m] + u[ku4+ku4_const+ku3+ku2+m]);
			}
		}
	});
	k=2;
  	ku4 = 3*ku4_const;
	std::for_each_n(policy, iter.front()+1, ny2, [&](int j){
		int ku3, ku2;
    	ku3 = j*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4-ku4_const+ku3+ku2+m]-=dssp* 
					(-4.0*u[ku4_const+ku3+ku2+m]+6.0*u[ku4-ku4_const+ku3+ku2+m]-
					 4.0*u[ku4+ku3+ku2+m]+u[ku4+ku4_const+ku3+ku2+m]);
			}
		}
	});
	std::for_each_n(policy, iter.front()+3, nz2-4, [&](int k){
		int ku4P, ku4M;
		int ku4, ku3, ku2;
		ku4 = k*ku4_const;
		ku4P = ku4+ku4_const;
		ku4M = ku4-ku4_const;
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]-=dssp* 
						(u[ku4M-ku4_const+ku3+ku2+m]-4.0*u[ku4M+ku3+ku2+m]+ 
						 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4P+ku3+ku2+m]+ 
						 u[ku4P+ku4_const+ku3+ku2+m]);
				}
			}
		}
	});
	k=nz2-1;
	ku4 = k*ku4_const;
	ku4P = ku4+ku4_const;
	ku4M = ku4-ku4_const;
	std::for_each_n(policy, iter.front()+1, ny2, [&](int j){
		int ku3, ku2;
    	ku3 = j*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]-=dssp*
					(u[ku4M-ku4_const+ku3+ku2+m]-4.0*u[ku4M+ku3+ku2+m]+ 
					 6.0*u[ku4+ku3+ku2+m]-4.0*u[ku4P+ku3+ku2+m]);
			}
		}
	});
	k=nz2;
	ku4 = k*ku4_const;
	ku4M = ku4-ku4_const;
	std::for_each_n(policy, iter.front()+1, ny2, [&](int j){
		int ku3, ku2;
    	ku3 = j*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<5; m++){
				rhs[ku4+ku3+ku2+m]-=dssp*
					(u[ku4M-ku4_const+ku3+ku2+m]
           -4.0*u[ku4M+ku3+ku2+m]+5.0*u[ku4+ku3+ku2+m]);
			}
		}
	});
	if(timeron){timer_stop(T_RHSZ);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		int ku4, ku3, ku2;
    	ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
    		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					rhs[ku4+ku3+ku2+m]*=dt;
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
	std::vector<double> u_exact(5, 0.0);
	double xi, eta, zeta, add;
	int ku4, ku3, ku2;
  	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	std::fill(rms.begin(), rms.begin()+5, 0.0);
	for(int k=0; k<=grid_points[2]-1; k++){
    	ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
      		ku3 = j*ku3_const;
			eta=(double)j*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
        		ku2 = i*5;
				xi=(double)i*dnxm1;
				exact_solution(xi, eta, zeta, u_exact);
				for(int m=0; m<5; m++){
					add=u[ku4+ku3+ku2+m]-u_exact[m];
					rms[m]+=add*add;
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
	int i, j, k, ip1, im1, jp1, jm1, km1, kp1;
  int ku4, ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3, ku3_const=(IMAXP+1)*5, ku2;
  int km2;
  int PST4=(PROBLEM_SIZE*4), PST3=(PROBLEM_SIZE*3), PST2=(PROBLEM_SIZE*2);
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
		zeta=(double)k*dnzm1;
		for(int j=1; j<=grid_points[1]-2; j++){
      		ku3 = j*ku3_const;
			eta=(double)j*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
				xi=(double)i*dnxm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
          			km2 = m*PROBLEM_SIZE;
					ue[km2+i]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
          			km2 = PROBLEM_SIZE*m;
					buf[km2+i]=dtpp*dtemp[m];
				}
				cuf[i]=buf[PROBLEM_SIZE+i]*buf[PROBLEM_SIZE+i];
				buf[i]=cuf[i]+buf[PST2+i]*buf[PST2+i]+buf[PST3+i]*buf[PST3+i]; 
				q[i]=0.5*(buf[PROBLEM_SIZE+i]*ue[PROBLEM_SIZE+i]+buf[PST2+i]*ue[PST2+i]+
						buf[PST3+i]*ue[PST3+i]);
			}
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				im1=i-1;
				ip1=i+1;
				forcing[ku4+ku3+ku2]-=
					tx2*(ue[PROBLEM_SIZE+ip1]-ue[PROBLEM_SIZE+im1])-
					dx1tx1*(ue[ip1]-2.0*ue[i]+ue[im1]);
				forcing[ku4+ku3+ku2+1]-=tx2*(
						(ue[PROBLEM_SIZE+ip1]*buf[PROBLEM_SIZE+ip1]+c2*(ue[PST4+ip1]-q[ip1]))-
						(ue[PROBLEM_SIZE+im1]*buf[PROBLEM_SIZE+im1]+c2*(ue[PST4+im1]-q[im1])))-
					xxcon1*(buf[PROBLEM_SIZE+ip1]-2.0*buf[PROBLEM_SIZE+i]+buf[PROBLEM_SIZE+im1])-
					dx2tx1*(ue[PROBLEM_SIZE+ip1]-2.0*ue[PROBLEM_SIZE+i]+ue[PROBLEM_SIZE+im1]);
				forcing[ku4+ku3+ku2+2]-=tx2*(
						ue[PST2+ip1]*buf[PROBLEM_SIZE+ip1]-ue[PST2+im1]*buf[PROBLEM_SIZE+im1])-
					xxcon2*(buf[PST2+ip1]-2.0*buf[PST2+i]+buf[PST2+im1])-
					dx3tx1*(ue[PST2+ip1]-2.0*ue[PST2+i]+ue[PST2+im1]);
				forcing[ku4+ku3+ku2+3]-=tx2*(
						ue[PST3+ip1]*buf[PROBLEM_SIZE+ip1]-ue[PST3+im1]*buf[PROBLEM_SIZE+im1])-
					xxcon2*(buf[PST3+ip1]-2.0*buf[PST3+i]+buf[PST3+im1])-
					dx4tx1*(ue[PST3+ip1]-2.0*ue[PST3+i]+ue[PST3+im1]);
				forcing[ku4+ku3+ku2+4]-=tx2*(
						buf[PROBLEM_SIZE+ip1]*(c1*ue[PST4+ip1]-c2*q[ip1])-
						buf[PROBLEM_SIZE+im1]*(c1*ue[PST4+im1]-c2*q[im1]))-
					0.5*xxcon3*(buf[ip1]-2.0*buf[i]+buf[im1])-
					xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])-
					xxcon5*(buf[PST4+ip1]-2.0*buf[PST4+i]+buf[PST4+im1])-
					dx5tx1*(ue[PST4+ip1]-2.0*ue[PST4+i]+ue[PST4+im1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                         
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				i=1;
				forcing[ku4+ku3+5+m]-=dssp*
					(5.0*ue[km2+i]-4.0*ue[km2+i+1]+ue[km2+i+2]);
				i=2;
				forcing[ku4+ku3+10+m]-=dssp*
					(-4.0*ue[km2+i-1]+6.0*ue[km2+i]-
					 4.0*ue[km2+i+1]+ue[km2+i+2]);
			}
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				for(int i=3; i<=grid_points[0]-4; i++){	
          			ku2 = i*5;			
					forcing[ku4+ku3+ku2+m]-=dssp*
						(ue[km2+i-2]-4.0*ue[km2+i-1]+
						 6.0*ue[km2+i]-4.0*ue[km2+i+1]+ue[km2+i+2]);
				}
			}
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				i=grid_points[0]-3;
        		ku2 = i*5;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+i-2]-4.0*ue[km2+i-1]+
					 6.0*ue[km2+i]-4.0*ue[km2+i+1]);
				i=grid_points[0]-2;
        		ku2 = i*5;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+i-2]-4.0*ue[km2+i-1]+5.0*ue[km2+i]);
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
		zeta=(double)k*dnzm1;
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			xi=(double)i*dnxm1;
			for(int j=0;j<=grid_points[1]-1;j++){
				eta=(double)j*dnym1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
					km2 = m*PROBLEM_SIZE;
					ue[km2+j]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
					km2 = m*PROBLEM_SIZE;
					buf[km2+j]=dtpp*dtemp[m];
				}
				cuf[j]=buf[PST2+j]*buf[PST2+j];
				buf[j]=cuf[j]+buf[PROBLEM_SIZE+j]*buf[PROBLEM_SIZE+j]+buf[PST3+j]*buf[PST3+j];
				q[j]=0.5*(buf[PROBLEM_SIZE+j]*ue[PROBLEM_SIZE+j]+buf[PST2+j]*ue[PST2+j]+
						buf[PST3+j]*ue[PST3+j]);
			}
			for(int j=1; j<=grid_points[1]-2; j++){
        		ku3 = j*ku3_const;
				jm1=j-1;
				jp1=j+1;
				forcing[ku4+ku3+ku2]-=
					ty2*(ue[PST2+jp1]-ue[PST2+jm1])-
					dy1ty1*(ue[jp1]-2.0*ue[j]+ue[jm1]);
				forcing[ku4+ku3+ku2+1]-=ty2*(
						ue[PROBLEM_SIZE+jp1]*buf[PST2+jp1]-ue[PROBLEM_SIZE+jm1]*buf[PST2+jm1])-
					yycon2*(buf[PROBLEM_SIZE+jp1]-2.0*buf[PROBLEM_SIZE+j]+buf[PROBLEM_SIZE+jm1])-
					dy2ty1*(ue[PROBLEM_SIZE+jp1]-2.0*ue[PROBLEM_SIZE+j]+ue[PROBLEM_SIZE+jm1]);
				forcing[ku4+ku3+ku2+2]-=ty2*(
						(ue[PST2+jp1]*buf[PST2+jp1]+c2*(ue[PST4+jp1]-q[jp1]))-
						(ue[PST2+jm1]*buf[PST2+jm1]+c2*(ue[PST4+jm1]-q[jm1])))-
					yycon1*(buf[PST2+jp1]-2.0*buf[PST2+j]+buf[PST2+jm1])-
					dy3ty1*(ue[PST2+jp1]-2.0*ue[PST2+j]+ue[PST2+jm1]);
				forcing[ku4+ku3+ku2+3]-=ty2*(
						ue[PST3+jp1]*buf[PST2+jp1]-ue[PST3+jm1]*buf[PST2+jm1])-
					yycon2*(buf[PST3+jp1]-2.0*buf[PST3+j]+buf[PST3+jm1])-
					dy4ty1*(ue[PST3+jp1]-2.0*ue[PST3+j]+ue[PST3+jm1]);
				forcing[ku4+ku3+ku2+4]-=ty2*(
						buf[PST2+jp1]*(c1*ue[PST4+jp1]-c2*q[jp1])-
						buf[PST2+jm1]*(c1*ue[PST4+jm1]-c2*q[jm1]))-
					0.5*yycon3*(buf[jp1]-2.0*buf[j]+
							buf[jm1])-
					yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])-
					yycon5*(buf[PST4+jp1]-2.0*buf[PST4+j]+buf[PST4+jm1])-
					dy5ty1*(ue[PST4+jp1]-2.0*ue[PST4+j]+ue[PST4+jm1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                      
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
				km2 = m*PROBLEM_SIZE;
				j=1;
				ku3 = ku3_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(5.0*ue[km2+j]-4.0*ue[km2+j+1]+ue[km2+j+2]);
				j=2;
				forcing[ku4+(ku3*2)+ku2+m]-=dssp*
					(-4.0*ue[km2+j-1]+6.0*ue[km2+j]-
					 4.0*ue[km2+j+1]+ue[km2+j+2]);
			}
			for(int m=0; m<5; m++){
				km2 = m*PROBLEM_SIZE;
				for(int j=3; j<=grid_points[1]-4; j++){	
					ku3 = j*ku3_const;			
					forcing[ku4+ku3+ku2+m]-=dssp*
						(ue[km2+j-2]-4.0*ue[km2+j-1]+
						 6.0*ue[km2+j]-4.0*ue[km2+j+1]+ue[km2+j+2]);
				}
			}
			for(int m=0; m<5; m++){
				km2 = m*PROBLEM_SIZE;
				j=grid_points[1]-3;
        		ku3 = j*ku3_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+j-2]-4.0*ue[km2+j-1]+
					 6.0*ue[km2+j]-4.0*ue[km2+j+1]);
				j=grid_points[1]-2;
        		ku3 = j*ku3_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+j-2]-4.0*ue[km2+j-1]+5.0*ue[km2+j]);
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
		eta=(double)j*dnym1;
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			xi=(double)i*dnxm1;
			for(int k=0; k<=grid_points[2]-1; k++){
				zeta=(double)k*dnzm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(int m=0; m<5; m++){
          			km2 = m*PROBLEM_SIZE;
					ue[km2+k]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(int m=1; m<5; m++){
          			km2 = m*PROBLEM_SIZE;
					buf[km2+k]=dtpp*dtemp[m];
				}
				cuf[k]=buf[PST3+k]*buf[PST3+k];
				buf[k]=cuf[k]+buf[PROBLEM_SIZE+k]*buf[PROBLEM_SIZE+k]+buf[PST2+k]*buf[PST2+k];
				q[k]=0.5*(buf[PROBLEM_SIZE+k]*ue[PROBLEM_SIZE+k]+buf[PST2+k]*ue[PST2+k]+
						buf[PST3+k]*ue[PST3+k]);
			}
			for(int k=1; k<=grid_points[2]-2; k++){
        		ku4 = k*ku4_const;
				km1=k-1;
				kp1=k+1;
				forcing[ku4+ku3+ku2]-=
					tz2*(ue[PST3+kp1]-ue[PST3+km1])-
					dz1tz1*(ue[kp1]-2.0*ue[k]+ue[km1]);
				forcing[ku4+ku3+ku2+1]-=tz2*(
						ue[PROBLEM_SIZE+kp1]*buf[PST3+kp1]-ue[PROBLEM_SIZE+km1]*buf[PST3+km1])-
					zzcon2*(buf[PROBLEM_SIZE+kp1]-2.0*buf[PROBLEM_SIZE+k]+buf[PROBLEM_SIZE+km1])-
					dz2tz1*(ue[PROBLEM_SIZE+kp1]-2.0*ue[PROBLEM_SIZE+k]+ue[PROBLEM_SIZE+km1]);
				forcing[ku4+ku3+ku2+2]-=tz2*(
						ue[PST2+kp1]*buf[PST3+kp1]-ue[PST2+km1]*buf[PST3+km1])-
					zzcon2*(buf[PST2+kp1]-2.0*buf[PST2+k]+buf[PST2+km1])-
					dz3tz1*(ue[PST2+kp1]-2.0*ue[PST2+k]+ue[PST2+km1]);
				forcing[ku4+ku3+ku2+3]-=tz2*(
						(ue[PST3+kp1]*buf[PST3+kp1]+c2*(ue[PST4+kp1]-q[kp1]))-
						(ue[PST3+km1]*buf[PST3+km1]+c2*(ue[PST4+km1]-q[km1])))-
					zzcon1*(buf[PST3+kp1]-2.0*buf[PST3+k]+buf[PST3+km1])-
					dz4tz1*(ue[PST3+kp1]-2.0*ue[PST3+k]+ue[PST3+km1]);
				forcing[ku4+ku3+ku2+4]-=tz2*(
						buf[PST3+kp1]*(c1*ue[PST4+kp1]-c2*q[kp1])-
						buf[PST3+km1]*(c1*ue[PST4+km1]-c2*q[km1]))-
					0.5*zzcon3*(buf[kp1]-2.0*buf[k]+buf[km1])-
					zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])-
					zzcon5*(buf[PST4+kp1]-2.0*buf[PST4+k]+buf[PST4+km1])-
					dz5tz1*(ue[PST4+kp1]-2.0*ue[PST4+k]+ue[PST4+km1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				k=1;
        		ku4 = ku4_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(5.0*ue[km2+k]-4.0*ue[km2+k+1]+ue[km2+k+2]);
				k=2;
				forcing[(ku4*2)+ku3+ku2+m]-=dssp*
					(-4.0*ue[km2+k-1]+6.0*ue[km2+k]-
					 4.0*ue[km2+k+1]+ue[km2+k+2]);
			}
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				for(int k=3; k<=grid_points[2]-4; k++){	
          			ku4 = k*ku4_const;		
					forcing[ku4+ku3+ku2+m]-=dssp*
						(ue[km2+k-2]-4.0*ue[km2+k-1]+
						 6.0*ue[km2+k]-4.0*ue[km2+k+1]+ue[km2+k+2]);
				}
			}
			for(int m=0; m<5; m++){
        		km2 = m*PROBLEM_SIZE;
				k=grid_points[2]-3;
        		ku4 = k*ku4_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+k-2]-4.0*ue[km2+k-1]+
					 6.0*ue[km2+k]-4.0*ue[km2+k+1]);
				k=grid_points[2]-2;
        		ku4 = k*ku4_const;
				forcing[ku4+ku3+ku2+m]-=dssp*
					(ue[km2+k-2]-4.0*ue[km2+k-1]+5.0*ue[km2+k]);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * now change the sign of the forcing function
	 * ---------------------------------------------------------------------
	 */
	for(int k=1; k<=grid_points[2]-2; k++){
		ku4 = k*ku4_const;
		for(int j=1; j<=grid_points[1]-2; j++){
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
				ku2 = i*5;
				for(int m=0; m<5; m++){
					forcing[ku4+ku3+ku2+m]=-1.0*forcing[ku4+ku3+ku2+m];
				}
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * this function returns the exact solution at point xi, eta, zeta  
 * ---------------------------------------------------------------------
 */
void exact_solution(double xi, double eta, double zeta, std::span<double> dtemp){
	for(int m=0; m<5; m++){
		dtemp[m]=ce[m]+xi*
			(ce[5+m]+xi*
			 (ce[4*5+m]+xi*
			  (ce[7*5+m]+xi*
			   ce[10*5+m])))+eta*
			(ce[2*5+m]+eta*
			 (ce[5*5+m]+eta*
			  (ce[8*5+m]+eta*
			   ce[11*5+m])))+zeta*
			(ce[3*5+m]+zeta*
			 (ce[6*5+m]+zeta*
			  (ce[9*5+m]+zeta*
			   ce[12*5+m])));
	}
}

/*
 * ---------------------------------------------------------------------
 * this subroutine initializes the field variable u using 
 * tri-linear transfinite interpolation of the boundary values     
 * ---------------------------------------------------------------------
 */
void initialize(){
  
	int i, j, k, ix, iy, iz;
	double xi, eta, zeta, Pxi, Peta, Pzeta;
	std::vector<double> temp(5), Pface(2*3*5);
	int kp3;
	int ku4, ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3, ku3_const=(IMAXP+1)*5, ku2;
	/*
	 * ---------------------------------------------------------------------
	 * later (in compute_rhs) we compute 1/u for every element. a few of 
	 * the corner elements are not used, but it convenient (and faster) 
	 * to compute the whole thing with a simple loop. make sure those 
	 * values are nonzero by initializing the whole thing here. 
	 * ---------------------------------------------------------------------
	 */
	for(int k=0; k<=grid_points[2]-1; k++){
		ku4 = k*ku4_const;
		for(int j=0; j<=grid_points[1]-1; j++){
      		ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-1; i++){
        		ku2 = i*5;
				u[ku4+ku3+ku2]=1.0;
				u[ku4+ku3+ku2+1]=0.0;
				u[ku4+ku3+ku2+2]=0.0;
				u[ku4+ku3+ku2+3]=0.0;
				u[ku4+ku3+ku2+4]=1.0;
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * first store the "interpolated" values everywhere on the grid    
	 * ---------------------------------------------------------------------
	 */
	for(int k=0; k<=grid_points[2]-1; k++){
    	ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
      		ku3 = j*ku3_const;
			eta=(double)j*dnym1;
			for(int i=0; i<=grid_points[0]-1; i++){
        		ku2 = i*5;
				xi=(double)i*dnxm1;
				for(int ix=0; ix<2; ix++){
					Pxi=(double)ix;
					exact_solution(Pxi, eta, zeta, {Pface.begin() + (ix*15), Pface.end()});
				}
				for(int iy=0; iy<2; iy++){
					Peta=(double)iy;
					exact_solution(xi, Peta, zeta, {Pface.begin() + (iy*15+5), Pface.end()});
				}
				for(int iz=0; iz<2; iz++){
					Pzeta=(double)iz;
					exact_solution(xi, eta, Pzeta, {Pface.begin() + (iz*15+10), Pface.end()});
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
	}
	/*
	 * ---------------------------------------------------------------------
	 * now store the exact values on the boundaries        
	 * ---------------------------------------------------------------------
	 * west face                                                  
	 * ---------------------------------------------------------------------
	 */
	xi=0.0;
	i=0;
	for(int k=0; k<=grid_points[2]-1; k++){
    	ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
      		ku3 = j*ku3_const;
			eta=(double)j*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+m]=temp[m];
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * east face                                                      
	 * ---------------------------------------------------------------------
	 */
	xi=1.0;
	i=grid_points[0]-1;
  	ku2 = i*5;
	for(int k=0; k<=grid_points[2]-1; k++){
    	ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int j=0; j<=grid_points[1]-1; j++){
      		ku3 = j*ku3_const;
			eta=(double)j*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * south face                                                 
	 * ---------------------------------------------------------------------
	 */
	eta=0.0;
	j=0;
	for(int k=0; k<=grid_points[2]-1; k++){
		ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int i=0; i<=grid_points[0]-1; i++){
      		ku2 = i*5;
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku2+m]=temp[m];
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * north face                                    
	 * ---------------------------------------------------------------------
	 */
	eta=1.0;
	j=grid_points[1]-1;
  	ku3 = j*ku3_const;
	for(int k=0; k<=grid_points[2]-1; k++){
    	ku4 = k*ku4_const;
		zeta=(double)k*dnzm1;
		for(int i=0; i<=grid_points[0]-1; i++){
      		ku2 = i*5;
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * bottom face                                       
	 * ---------------------------------------------------------------------
	 */
	zeta=0.0;
	k=0;
	for(int j=0; j<=grid_points[1]-1; j++){
    	ku3 = j*ku3_const;
		eta=(double)j*dnym1;
		for(int i=0; i<=grid_points[0]-1; i++){
      		ku2 = i*5;
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku3+ku2+m]=temp[m];
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * top face     
	 * ---------------------------------------------------------------------
	 */
	zeta=1.0;
	k=grid_points[2]-1;
  	ku4 = k*ku4_const;
	for(int j=0; j<=grid_points[1]-1; j++){
    	ku3 = j*ku3_const;
		eta=(double)j*dnym1;
		for(int i=0; i<=grid_points[0]-1; i++){
      		ku2 = i*5;
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(int m=0; m<5; m++){
				u[ku4+ku3+ku2+m]=temp[m];
			}
		}
	}
}

void lhsinit(int ni, int nj){
  int kl3, kl3_const=(IMAXP+1)*5, kl2=ni*5;
	/*
	 * ---------------------------------------------------------------------
	 * zap the whole left hand side for starters
	 * set all diagonal values to 1. This is overkill, but convenient
	 * ---------------------------------------------------------------------
	 */
	for(int j=1; j<=nj; j++){
    	kl3 = j*kl3_const;
		for(int m=0; m<5; m++){
			lhs[kl3+m]=0.0;
			lhsp[kl3+m]=0.0;
			lhsm[kl3+m]=0.0;
			lhs[kl3+kl2+m]=0.0;
			lhsp[kl3+kl2+m]=0.0;
			lhsm[kl3+kl2+m]=0.0;
		}
		lhs[kl3+2]=1.0;
		lhsp[kl3+2]=1.0;
		lhsm[kl3+2]=1.0;
		lhs[kl3+kl2+2]=1.0;
		lhsp[kl3+kl2+2]=1.0;
		lhsm[kl3+kl2+2]=1.0;
	}
}

void lhsinitj(int nj, int ni){
  int kl3=(IMAXP+1)*5*nj, kl2;
	/*
	 * ---------------------------------------------------------------------
	 * zap the whole left hand side for starters
	 * set all diagonal values to 1. This is overkill, but convenient
	 * ---------------------------------------------------------------------
	 */
	for(int i=1; i<=ni; i++){
    	kl2 = i*5;
		for(int m=0; m<5; m++){
			lhs[kl2+m]=0.0;
			lhsp[kl2+m]=0.0;
			lhsm[kl2+m]=0.0;
			lhs[kl3+kl2+m]=0.0;
			lhsp[kl3+kl2+m]=0.0;
			lhsm[kl3+kl2+m]=0.0;
		}
		lhs[kl2+2]=1.0;
		lhsp[kl2+2]=1.0;
		lhsm[kl2+2]=1.0;
		lhs[kl3+kl2+2]=1.0;
		lhsp[kl3+kl2+2]=1.0;
		lhsm[kl3+kl2+2]=1.0;
	}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication              
 * ---------------------------------------------------------------------
 */
void ninvr(){
	if(timeron){timer_start(T_NINVR);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		double r1, r2, r3, r4, r5, t1, t2;
		int ku4, ku3, ku2;
		int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
 
		ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				r1=rhs[ku4+ku3+ku2];
				r2=rhs[ku4+ku3+ku2+1];
				r3=rhs[ku4+ku3+ku2+2];
				r4=rhs[ku4+ku3+ku2+3];
				r5=rhs[ku4+ku3+ku2+4];
				t1=bt*r3;
				t2=0.5*(r4+r5);
				rhs[ku4+ku3+ku2+0]=-r2;
				rhs[ku4+ku3+ku2+1]=r1;
				rhs[ku4+ku3+ku2+2]=bt*(r4-r5);
				rhs[ku4+ku3+ku2+3]=-t1+t2;
				rhs[ku4+ku3+ku2+4]=t1+t2;
			}
		}
	});
	if(timeron){timer_stop(T_NINVR);}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication                       
 * ---------------------------------------------------------------------
 */
void pinvr(){
	if(timeron){timer_start(T_PINVR);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		double r1, r2, r3, r4, r5, t1, t2;
		int ku4, ku3, ku2;
		int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
	
    	ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				r1=rhs[ku4+ku3+ku2];
				r2=rhs[ku4+ku3+ku2+1];
				r3=rhs[ku4+ku3+ku2+2];
				r4=rhs[ku4+ku3+ku2+3];
				r5=rhs[ku4+ku3+ku2+4];
				t1=bt*r1;
				t2=0.5*(r4+r5);
				rhs[ku4+ku3+ku2+0]=bt*(r4-r5);
				rhs[ku4+ku3+ku2+1]=-r3;
				rhs[ku4+ku3+ku2+2]=r2;
				rhs[ku4+ku3+ku2+3]=-t1+t2;
				rhs[ku4+ku3+ku2+4]=t1+t2;
			}
		}
	});
	if(timeron){timer_stop(T_PINVR);}
}

void rhs_norm(std::vector<double> &rms){
	double add;
	int ku4, ku3, ku2;
	int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;

	std::fill(rms.begin(), rms.end(), 0.0);
	for(int k=1; k<=nz2; k++){
    	ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m=0; m<5; m++){
					add=rhs[ku4+ku3+ku2+m];
					rms[m]+=add*add;
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

void set_constants(){
	ce[0]=2.0;
	ce[5]=0.0;
	ce[2*5]=0.0;
	ce[3*5]=4.0;
	ce[4*5]=5.0;
	ce[5*5]=3.0;
	ce[6*5]=0.5;
	ce[7*5]=0.02;
	ce[8*5]=0.01;
	ce[9*5]=0.03;
	ce[10*5]=0.5;
	ce[11*5]=0.4;
	ce[12*5]=0.3;
	/* */
	ce[0*5+1]=1.0;
	ce[1*5+1]=0.0;
	ce[2*5+1]=0.0;
	ce[3*5+1]=0.0;
	ce[4*5+1]=1.0;
	ce[5*5+1]=2.0;
	ce[6*5+1]=3.0;
	ce[7*5+1]=0.01;
	ce[8*5+1]=0.03;
	ce[9*5+1]=0.02;
	ce[10*5+1]=0.4;
	ce[11*5+1]=0.3;
	ce[12*5+1]=0.5;
	/* */
	ce[2]=2.0;
	ce[5+2]=2.0;
	ce[2*5+2]=0.0;
	ce[3*5+2]=0.0;
	ce[4*5+2]=0.0;
	ce[5*5+2]=2.0;
	ce[6*5+2]=3.0;
	ce[7*5+2]=0.04;
	ce[8*5+2]=0.03;
	ce[9*5+2]=0.05;
	ce[10*5+2]=0.3;
	ce[11*5+2]=0.5;
	ce[12*5+2]=0.4;
	/* */
	ce[3]=2.0;
	ce[5+3]=2.0;
	ce[2*5+3]=0.0;
	ce[3*5+3]=0.0;
	ce[4*5+3]=0.0;
	ce[5*5+3]=2.0;
	ce[6*5+3]=3.0;
	ce[7*5+3]=0.03;
	ce[8*5+3]=0.05;
	ce[9*5+3]=0.04;
	ce[10*5+3]=0.2;
	ce[11*5+3]=0.1;
	ce[12*5+3]=0.3;
	/* */
	ce[4]=5.0;
	ce[5+4]=4.0;
	ce[2*5+4]=3.0;
	ce[3*5+4]=2.0;
	ce[4*5+4]=0.1;
	ce[5*5+4]=0.4;
	ce[6*5+4]=0.3;
	ce[7*5+4]=0.05;
	ce[8*5+4]=0.04;
	ce[9*5+4]=0.03;
	ce[10*5+4]=0.1;
	ce[11*5+4]=0.3;
	ce[12*5+4]=0.2;
	/* */
	c1=1.4;
	c2=0.4;
	c3=0.1;
	c4=1.0;
	c5=1.4;
	/* */
	bt=sqrt(0.5);
	/* */
	dnxm1=1.0/(double)(grid_points[0]-1);
	dnym1=1.0/(double)(grid_points[1]-1);
	dnzm1=1.0/(double)(grid_points[2]-1);
	/* */
	c1c2=c1*c2;
	c1c5=c1*c5;
	c3c4=c3*c4;
	c1345=c1c5*c3c4;
	/* */
	conz1=(1.0-c1c5);
	/* */
	tx1=1.0/(dnxm1*dnxm1);
	tx2=1.0/(2.0*dnxm1);
	tx3=1.0/dnxm1;
	/* */
	ty1=1.0/(dnym1*dnym1);
	ty2=1.0/(2.0*dnym1);
	ty3=1.0/dnym1;
	/* */
	tz1=1.0/(dnzm1*dnzm1);
	tz2=1.0/(2.0*dnzm1);
	tz3=1.0/dnzm1;
	/* */
	dx1=0.75;
	dx2=0.75;
	dx3=0.75;
	dx4=0.75;
	dx5=0.75;
	/* */
	dy1=0.75;
	dy2=0.75;
	dy3=0.75;
	dy4=0.75;
	dy5=0.75;
	/* */
	dz1=1.0;
	dz2=1.0;
	dz3=1.0;
	dz4=1.0;
	dz5=1.0;
	/* */
	dxmax=std::max(dx3, dx4);
	dymax=std::max(dy2, dy4);
	dzmax=std::max(dz2, dz3);
	/* */
	dssp=0.25*std::max(dx1, std::max(dy1, dz1));
	/* */
	c4dssp=4.0*dssp;
	c5dssp=5.0*dssp;
	/* */
	dttx1=dt*tx1;
	dttx2=dt*tx2;
	dtty1=dt*ty1;
	dtty2=dt*ty2;
	dttz1=dt*tz1;
	dttz2=dt*tz2;
	/* */
	c2dttx1=2.0*dttx1;
	c2dtty1=2.0*dtty1;
	c2dttz1=2.0*dttz1;
	/* */
	dtdssp=dt*dssp;
	/* */
	comz1=dtdssp;
	comz4=4.0*dtdssp;
	comz5=5.0*dtdssp;
	comz6=6.0*dtdssp;
	/* */
	c3c4tx3=c3c4*tx3;
	c3c4ty3=c3c4*ty3;
	c3c4tz3=c3c4*tz3;
	/* */
	dx1tx1=dx1*tx1;
	dx2tx1=dx2*tx1;
	dx3tx1=dx3*tx1;
	dx4tx1=dx4*tx1;
	dx5tx1=dx5*tx1;
	/* */
	dy1ty1=dy1*ty1;
	dy2ty1=dy2*ty1;
	dy3ty1=dy3*ty1;
	dy4ty1=dy4*ty1;
	dy5ty1=dy5*ty1;
	/* */
	dz1tz1=dz1*tz1;
	dz2tz1=dz2*tz1;
	dz3tz1=dz3*tz1;
	dz4tz1=dz4*tz1;
	dz5tz1=dz5*tz1;
	/* */
	c2iv=2.5;
	con43=4.0/3.0;
	con16=1.0/6.0;
	/* */
	xxcon1=c3c4tx3*con43*tx3;
	xxcon2=c3c4tx3*tx3;
	xxcon3=c3c4tx3*conz1*tx3;
	xxcon4=c3c4tx3*con16*tx3;
	xxcon5=c3c4tx3*c1c5*tx3;
	/* */
	yycon1=c3c4ty3*con43*ty3;
	yycon2=c3c4ty3*ty3;
	yycon3=c3c4ty3*conz1*ty3;
	yycon4=c3c4ty3*con16*ty3;
	yycon5=c3c4ty3*c1c5*ty3;
	/* */
	zzcon1=c3c4tz3*con43*tz3;
	zzcon2=c3c4tz3*tz3;
	zzcon3=c3c4tz3*conz1*tz3;
	zzcon4=c3c4tz3*con16*tz3;
	zzcon5=c3c4tz3*c1c5*tz3;
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication                  
 * ---------------------------------------------------------------------
 */
void txinvr(){

	if(timeron){timer_start(T_TXINVR);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		double t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, r4, r5, ac2inv;
		int ku4, ku3, ku2, ks3, ks2;
		int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
		int ks3_const=(JMAXP+1)*(IMAXP+1);
		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
			ks2 = j*(IMAXP+1);
			ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				ru1=rho_i[ks3+ks2+i];
				uu=us[ks3+ks2+i];
				vv=vs[ks3+ks2+i];
				ww=ws[ks3+ks2+i];
				ac=speed[ks3+ks2+i];
				ac2inv=ac*ac;
				r1=rhs[ku4+ku3+ku2];
				r2=rhs[ku4+ku3+ku2+1];
				r3=rhs[ku4+ku3+ku2+2];
				r4=rhs[ku4+ku3+ku2+3];
				r5=rhs[ku4+ku3+ku2+4];
				t1=c2/ac2inv*(qs[ks3+ks2+i]*r1-uu*r2-vv*r3-ww*r4+r5);
				t2=bt*ru1*(uu*r1-r2);
				t3=(bt*ru1*ac)*t1;
				rhs[ku4+ku3+ku2]=r1-t1;
				rhs[ku4+ku3+ku2+1]=-ru1*(ww*r1-r4);
				rhs[ku4+ku3+ku2+2]=ru1*(vv*r1-r3);
				rhs[ku4+ku3+ku2+3]=-t2+t3;
				rhs[ku4+ku3+ku2+4]=t2+t3;
			}
		}
	});
	if(timeron){timer_stop(T_TXINVR);}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication                       
 * ---------------------------------------------------------------------
 */
void tzetar(){

  if(timeron){timer_start(T_TZETAR);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		double t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3, r4, r5, btuz, ac2u, uzik1;
		int ku4, ku3, ku2, ks3, ks2;
		int ku4_const=(JMAXP+1)*(IMAXP+1)*5, ku3_const=(IMAXP+1)*5;
		int ks3_const=(JMAXP+1)*(IMAXP+1);

		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		for(int j=1; j<=ny2; j++){
			ks2 = j*(IMAXP+1);
			ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				xvel=us[ks3+ks2+i];
				yvel=vs[ks3+ks2+i];
				zvel=ws[ks3+ks2+i];
				ac=speed[ks3+ks2+i];
				ac2u=ac*ac;
				r1=rhs[ku4+ku3+ku2];
				r2=rhs[ku4+ku3+ku2+1];
				r3=rhs[ku4+ku3+ku2+2];
				r4=rhs[ku4+ku3+ku2+3];
				r5=rhs[ku4+ku3+ku2+4];
				uzik1=u[ku4+ku3+ku2];
				btuz=bt*uzik1;
				t1=btuz/ac*(r4+r5);
				t2=r3+t1;
				t3=btuz*(r4-r5);
				rhs[ku4+ku3+ku2]=t2;
				rhs[ku4+ku3+ku2+1]=-uzik1*r2+xvel*t2;
				rhs[ku4+ku3+ku2+2]=uzik1*r1+yvel*t2;
				rhs[ku4+ku3+ku2+3]=zvel*t2+t3;
				rhs[ku4+ku3+ku2+4]=uzik1*(-xvel*r2+yvel*r1) + 
					qs[ks3+ks2+i]*t2+c2iv*ac2u*t1+zvel*t3;
			}
		}
	});
	if(timeron){timer_stop(T_TZETAR);}
}

/*
 * ---------------------------------------------------------------------
 * verification routine                         
 * ---------------------------------------------------------------------
 */
void verify(int no_time_steps, char* class_npb, boolean* verified){
	double epsilon, dtref;
	std::vector<double> xce(5), xcr(5), xcrref(5), xceref(5), xcrdif(5), xcedif(5);
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
	 * reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
	 * ---------------------------------------------------------------------
	 */
	if((grid_points[0]==12)&&(grid_points[1]==12)&&(grid_points[2]==12)&&(no_time_steps==100)){
		*class_npb='S';
		dtref=1.5e-2;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=2.7470315451339479e-02;
		xcrref[1]=1.0360746705285417e-02;
		xcrref[2]=1.6235745065095532e-02;
		xcrref[3]=1.5840557224455615e-02;
		xcrref[4]=3.4849040609362460e-02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=2.7289258557377227e-05;
		xceref[1]=1.0364446640837285e-05;
		xceref[2]=1.6154798287166471e-05;
		xceref[3]=1.5750704994480102e-05;
		xceref[4]=3.4177666183390531e-05;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 36X36X36 grids after 400 time steps, with DT = 1.5d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==36)&&(grid_points[1]==36)&&(grid_points[2]==36)&&(no_time_steps==400)){
		*class_npb='W';
		dtref=1.5e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.1893253733584e-02;
		xcrref[1]=0.1717075447775e-03;
		xcrref[2]=0.2778153350936e-03;
		xcrref[3]=0.2887475409984e-03;
		xcrref[4]=0.3143611161242e-02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.7542088599534e-04;
		xceref[1]=0.6512852253086e-05;
		xceref[2]=0.1049092285688e-04;
		xceref[3]=0.1128838671535e-04;
		xceref[4]=0.1212845639773e-03;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 64X64X64 grids after 400 time steps, with DT = 1.5d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==64)&&(grid_points[1]==64)&&(grid_points[2]==64)&&(no_time_steps==400)){
		*class_npb='A';
		dtref=1.5e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual.
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=2.4799822399300195;
		xcrref[1]=1.1276337964368832;
		xcrref[2]=1.5028977888770491;
		xcrref[3]=1.4217816211695179;
		xcrref[4]=2.1292113035138280;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error.
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=1.0900140297820550e-04;
		xceref[1]=3.7343951769282091e-05;
		xceref[2]=5.0092785406541633e-05;
		xceref[3]=4.7671093939528255e-05;
		xceref[4]=1.3621613399213001e-04;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 102X102X102 grids after 400 time steps,
		 * with DT = 1.0d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==102)&&(grid_points[1]==102)&&(grid_points[2]==102)&&(no_time_steps==400)){
		*class_npb='B';
		dtref=1.0e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.6903293579998e+02;
		xcrref[1]=0.3095134488084e+02;
		xcrref[2]=0.4103336647017e+02;
		xcrref[3]=0.3864769009604e+02;
		xcrref[4]=0.5643482272596e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.9810006190188e-02;
		xceref[1]=0.1022827905670e-02;
		xceref[2]=0.1720597911692e-02;
		xceref[3]=0.1694479428231e-02;
		xceref[4]=0.1847456263981e-01;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 162X162X162 grids after 400 time steps,
		 * with DT = 0.67d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==162)&&(grid_points[1]==162)&&(grid_points[2]==162)&&(no_time_steps==400)){
		*class_npb='C';
		dtref=0.67e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.5881691581829e+03;
		xcrref[1]=0.2454417603569e+03;
		xcrref[2]=0.3293829191851e+03;
		xcrref[3]=0.3081924971891e+03;
		xcrref[4]=0.4597223799176e+03;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.2598120500183e+00;
		xceref[1]=0.2590888922315e-01;
		xceref[2]=0.5132886416320e-01;
		xceref[3]=0.4806073419454e-01;
		xceref[4]=0.5483377491301e+00;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 408X408X408 grids after 500 time steps,
		 * with DT = 0.3d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==408)&&(grid_points[1]==408)&&(grid_points[2]==408)&&(no_time_steps==500)){
		*class_npb='D';
		dtref=0.30e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.1044696216887e+05;
		xcrref[1]=0.3204427762578e+04;
		xcrref[2]=0.4648680733032e+04;
		xcrref[3]=0.4238923283697e+04;
		xcrref[4]=0.7588412036136e+04;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.5089471423669e+01;
		xceref[1]=0.5323514855894e+00;
		xceref[2]=0.1187051008971e+01;
		xceref[3]=0.1083734951938e+01;
		xceref[4]=0.1164108338568e+02;
		/*
		 * ---------------------------------------------------------------------
		 * reference data for 1020X1020X1020 grids after 500 time steps,
		 * with DT = 0.1d-03
		 * ---------------------------------------------------------------------
		 */
	}else if((grid_points[0]==1020)&&(grid_points[1]==1020)&&(grid_points[2]==1020)&&(no_time_steps==500)){
		*class_npb='E';
		dtref=0.10e-3;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of residual
		 * ---------------------------------------------------------------------
		 */
		xcrref[0]=0.6255387422609e+05;
		xcrref[1]=0.1495317020012e+05;
		xcrref[2]=0.2347595750586e+05;
		xcrref[3]=0.2091099783534e+05;
		xcrref[4]=0.4770412841218e+05;
		/*
		 * ---------------------------------------------------------------------
		 * reference values of RMS-norms of solution error
		 * ---------------------------------------------------------------------
		 */
		xceref[0]=0.6742735164909e+02;
		xceref[1]=0.5390656036938e+01;
		xceref[2]=0.1680647196477e+02;
		xceref[3]=0.1536963126457e+02;
		xceref[4]=0.1575330146156e+03;
	}else{
		*verified=FALSE;
	}
	/*
	 * ---------------------------------------------------------------------
	 * verification test for residuals if gridsize is one of 
	 * the defined grid sizes above (class .ne. 'U')
	 * ---------------------------------------------------------------------
	 * compute the difference of solution values and the known reference values
	 * ---------------------------------------------------------------------
	 */
	for(int m=0; m<5; m++){
		xcrdif[m]=fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m]=fabs((xce[m]-xceref[m])/xceref[m]);
	}
	/*
	 * ---------------------------------------------------------------------
	 * output the comparison of computed results to known cases
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
	for(int m=0;m<5;m++){
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
 * this function performs the solution of the approximate factorization
 * step in the x-direction for all five matrix components
 * simultaneously. the thomas algorithm is employed to solve the
 * systems for the x-lines. boundary conditions are non-periodic
 * ---------------------------------------------------------------------
 */
void x_solve(){
	if(timeron){timer_start(T_XSOLVE);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		thread_local std::vector<double> cv(PROBLEM_SIZE);
		thread_local std::vector<double> rhon(PROBLEM_SIZE);
		thread_local std::vector<double> lhs((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsp((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsm((IMAXP+1)*(IMAXP+1)*5);
		int ku4, ku3, ku2, ks3, ks2;
		
		// int ku3_const=(IMAXP+1)*5;
		int ku3_const=(IMAXP+1)*5;
		int ku4_const=(IMAXP+1)*ku3_const;
		int ks3_const=(IMAXP+1)*(IMAXP+1);
		int i, j, i1, i2, m;
		double ru1, fac1, fac2;
		
		for(int j=1; j<=ny2; j++){
			int ku2 = (nx2+1)*5;
			int ku3 = j*ku3_const;
			for(int m=0; m<5; m++){
				lhs[ku3 + m]=0.0;
				lhsp[ku3 + m]=0.0;
				lhsm[ku3 + m]=0.0;
				lhs[ku3 + ku2 + m]=0.0;
				lhsp[ku3 + ku2 + m]=0.0;
				lhsm[ku3 + ku2 + m]=0.0;
			}
			lhs	[ku3 + 2]=1.0;
			lhsp[ku3 + 2]=1.0;
			lhsm[ku3 + 2]=1.0;
			lhs	[ku3 + ku2 + 2]=1.0;
			lhsp[ku3 + ku2 + 2]=1.0;
			lhsm[ku3 + ku2 + 2]=1.0;
		}

		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three x-factors  
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue                   
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ks2 = j*(IMAXP+1);
			ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-1; i++){
				ru1=c3c4*rho_i[ks3+ks2+i];
				cv[i]=us[ks3+ks2+i];
				rhon[i]=std::max(std::max(dx2+con43*ru1,dx5+c1c5*ru1), std::max(dxmax+ru1,dx1));
			}
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				lhs[ku3+ku2]=0.0;
				lhs[ku3+ku2+1]=-dttx2*cv[i-1]-dttx1*rhon[i-1];
				lhs[ku3+ku2+2]=1.0+c2dttx1*rhon[i];
				lhs[ku3+ku2+3]=dttx2*cv[i+1]-dttx1*rhon[i+1];
				lhs[ku3+ku2+4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                             
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			i=1;
			lhs[ku3+5+2]+=comz5;
			lhs[ku3+5+3]-=comz4;
			lhs[ku3+5+4]+=comz1;
			lhs[ku3+10+1]-=comz4;
			lhs[ku3+10+2]+=comz6;
			lhs[ku3+10+3]-=comz4;
			lhs[ku3+10+4]+=comz1;
		}
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			for(int i=3; i<=grid_points[0]-4; i++){
        		ku2 = i*5;
				lhs[ku3+ku2]+=comz1;
				lhs[ku3+ku2+1]-=comz4;
				lhs[ku3+ku2+2]+=comz6;
				lhs[ku3+ku2+3]-=comz4;
				lhs[ku3+ku2+4]+=comz1;
			}
		}
		for(int j=1; j<=ny2; j++){
      		ku3 = j*ku3_const;
			i=grid_points[0]-3;
      		ku2 = i*5;
			lhs[ku3+ku2]+=comz1;
			lhs[ku3+ku2+1]-=comz4;
			lhs[ku3+ku2+2]+=comz6;
			lhs[ku3+ku2+3]-=comz4;
			lhs[ku3+ku2+5]+=comz1;
			lhs[ku3+ku2+5+1]-=comz4;
			lhs[ku3+ku2+5+2]+=comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, fill the other factors (u+c), (u-c) by adding to 
		 * the first  
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ks2 = j*(IMAXP+1);
			ku3 = j*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				lhsp[ku3+ku2]=lhs[ku3+ku2];
				lhsp[ku3+ku2+1]=lhs[ku3+ku2+1]-dttx2*speed[ks3+ks2+i-1];
				lhsp[ku3+ku2+2]=lhs[ku3+ku2+2];
				lhsp[ku3+ku2+3]=lhs[ku3+ku2+3]+dttx2*speed[ks3+ks2+i+1];
				lhsp[ku3+ku2+4]=lhs[ku3+ku2+4];
				lhsm[ku3+ku2]=lhs[ku3+ku2];
				lhsm[ku3+ku2+1]=lhs[ku3+ku2+1]+dttx2*speed[ks3+ks2+i-1];
				lhsm[ku3+ku2+2]=lhs[ku3+ku2+2];
				lhsm[ku3+ku2+3]=lhs[ku3+ku2+3]-dttx2*speed[ks3+ks2+i+1];
				lhsm[ku3+ku2+4]=lhs[ku3+ku2+4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 * perform the thomas algorithm; first, FORWARD ELIMINATION     
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-3; i++){
        		ku2 = i*5;
				i1=i+1;
				i2=i+2;
				fac1=1.0/lhs[ku3+ku2+2];
				lhs[ku3+ku2+3]*=fac1;
				lhs[ku3+ku2+4]*=fac1;
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+m]*=fac1;
				}
				lhs[ku3+ku2+5+2]-=lhs[ku3+ku2+5+1]*lhs[ku3+ku2+3];
				lhs[ku3+ku2+5+3]-=lhs[ku3+ku2+5+1]*lhs[ku3+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+5+m]-=lhs[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
				}
				lhs[ku3+ku2+10+1]-=lhs[ku3+ku2+10]*lhs[ku3+ku2+3];
				lhs[ku3+ku2+10+2]-=lhs[ku3+ku2+10]*lhs[ku3+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+10+m]-=lhs[ku3+ku2+10]*rhs[ku4+ku3+ku2+m];
				}
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * the last two rows in this grid block are a bit different, 
		 * since they do not have two more rows available for the
		 * elimination of off-diagonal entries
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			i=grid_points[0]-2;
      		ku2 = i*5;
			fac1=1.0/lhs[ku3+ku2+2];
			lhs[ku3+ku2+3]*=fac1;
			lhs[ku3+ku2+4]*=fac1;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+m]*=fac1;
			}
			i1=grid_points[0]-1;
			lhs[ku3+ku2+5+2]-=lhs[ku3+ku2+5+1]*lhs[ku3+ku2+3];
			lhs[ku3+ku2+5+3]-=lhs[ku3+ku2+5+1]*lhs[ku3+ku2+4];
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+5+m]-=lhs[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			fac2 = 1.0/lhs[ku3+ku2+5+2];
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+5+m]*=fac2;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors                 
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			for(int i=0; i<=grid_points[0]-3; i++){
        		ku2 = i*5;
				i1=i+1;
				i2=i+2;
				m=3;
				fac1=1.0/lhsp[ku3+ku2+2];
				lhsp[ku3+ku2+3]*=fac1;
				lhsp[ku3+ku2+4]*=fac1;
				rhs[ku4+ku3+ku2+m]*=fac1;
				lhsp[ku3+ku2+5+2]-=lhsp[ku3+ku2+5+1]*lhsp[ku3+ku2+3];
				lhsp[ku3+ku2+5+3]-=lhsp[ku3+ku2+5+1]*lhsp[ku3+ku2+4];
				rhs[ku4+ku3+ku2+5+m]-=lhsp[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
				lhsp[ku3+ku2+10+1]-=lhsp[ku3+ku2+10]*lhsp[ku3+ku2+3];
				lhsp[ku3+ku2+10+2]-=lhsp[ku3+ku2+10]*lhsp[ku3+ku2+4];
				rhs[ku4+ku3+ku2+10+m]-=lhsp[ku3+ku2+10]*rhs[ku4+ku3+ku2+m];
				m=4;
				fac1=1.0/lhsm[ku3+ku2+2];
				lhsm[ku3+ku2+3]*=fac1;
				lhsm[ku3+ku2+4]*=fac1;
				rhs[ku4+ku3+ku2+m]*=fac1;
				lhsm[ku3+ku2+5+2]-=lhsm[ku3+ku2+5+1]*lhsm[ku3+ku2+3];
				lhsm[ku3+ku2+5+3]-=lhsm[ku3+ku2+5+1]*lhsm[ku3+ku2+4];
				rhs[ku4+ku3+ku2+5+m]-=lhsm[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
				lhsm[ku3+ku2+10+1]-=lhsm[ku3+ku2+10]*lhsm[ku3+ku2+3];
				lhsm[ku3+ku2+10+2]-=lhsm[ku3+ku2+10]*lhsm[ku3+ku2+4];
				rhs[ku4+ku3+ku2+10+m]-=lhsm[ku3+ku2+10]*rhs[ku4+ku3+ku2+m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			i=grid_points[0]-2;
			i1=grid_points[0]-1;
      		ku2 = i*5;
			m=3;
			fac1=1.0/lhsp[ku3+ku2+2];
			lhsp[ku3+ku2+3]*=fac1;
			lhsp[ku3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsp[ku3+ku2+5+2]-=lhsp[ku3+ku2+5+1]*lhsp[ku3+ku2+3];
			lhsp[ku3+ku2+5+3]-=lhsp[ku3+ku2+5+1]*lhsp[ku3+ku2+4];
			rhs[ku4+ku3+ku2+5+m]-=lhsp[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
			m=4;
			fac1=1.0/lhsm[ku3+ku2+2];
			lhsm[ku3+ku2+3]*=fac1;
			lhsm[ku3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsm[ku3+ku2+5+2]-=lhsm[ku3+ku2+5+1]*lhsm[ku3+ku2+3];
			lhsm[ku3+ku2+5+3]-=lhsm[ku3+ku2+5+1]*lhsm[ku3+ku2+4];
			rhs[ku4+ku3+ku2+5+m]-=lhsm[ku3+ku2+5+1]*rhs[ku4+ku3+ku2+m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately
			 * ---------------------------------------------------------------------
			 */
			rhs[ku4+ku3+ku2+5+3]/=lhsp[ku3+ku2+5+2];
			rhs[ku4+ku3+ku2+5+4]/=lhsm[ku3+ku2+5+2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			i=grid_points[0]-2;
			i1=grid_points[0]-1;
      		ku2 = i*5;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+5+m]-=lhs[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+m];
			}
			rhs[ku4+ku3+ku2+5+3]-=lhsp[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+3];
			rhs[ku4+ku3+ku2+5+4]-=lhsm[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=ny2; j++){
			ku3 = j*ku3_const;
			for(int i=grid_points[0]-3; i>=0; i--){
        		ku2 = i*5;
				i1=i+1;
				i2=i+2;
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+m]-=
						lhs[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+m]+
						lhs[ku3+ku2+4]*rhs[ku4+ku3+ku2+10+m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[ku4+ku3+ku2+3]-=
					lhsp[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+3]+
					lhsp[ku3+ku2+4]*rhs[ku4+ku3+ku2+10+3];
				rhs[ku4+ku3+ku2+4]-=
					lhsm[ku3+ku2+3]*rhs[ku4+ku3+ku2+5+4]+
					lhsm[ku3+ku2+4]*rhs[ku4+ku3+ku2+10+4];
			}
		}
	});
	if(timeron){timer_stop(T_XSOLVE);}
	/*
	 * ---------------------------------------------------------------------
	 * do the block-diagonal inversion          
	 * ---------------------------------------------------------------------
	 */
	ninvr();
}

/*
 * ---------------------------------------------------------------------
 * this function performs the solution of the approximate factorization
 * step in the y-direction for all five matrix components
 * simultaneously. the thomas algorithm is employed to solve the
 * systems for the y-lines. boundary conditions are non-periodic
 * ---------------------------------------------------------------------
 */
void y_solve(){
	if(timeron){timer_start(T_YSOLVE);}
	std::for_each_n(policy, iter.front()+1, nz2, [&](int k){
		thread_local std::vector<double> cv(PROBLEM_SIZE);
		thread_local std::vector<double> rhoq(PROBLEM_SIZE);
		thread_local std::vector<double> lhs((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsp((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsm((IMAXP+1)*(IMAXP+1)*5);
		int ku4, ku3, ku2, ks3, ks2;
		int ku3_const=(IMAXP+1)*5;
		int ku4_const=(IMAXP+1)*ku3_const;
		int ks3_const=(IMAXP+1)*(IMAXP+1);
		int i, j, j1, j2, m;
		double ru1, fac1, fac2;
		
		for(int i=1; i<=nx2; i++){
			int ku2 = i*5;
			int ku3 = (nz2+1)*ku3_const;
			for(m=0; m<5; m++){
				lhs[ku2 + m]=0.0;
				lhsp[ku2 + m]=0.0;
				lhsm[ku2 + m]=0.0;
				lhs[ku3 + ku2 + m]=0.0;
				lhsp[ku3 + ku2 + m]=0.0;
				lhsm[ku3 + ku2 + m]=0.0;
			}
			lhs	[ku2 + 2]=1.0;
			lhsp[ku2 + 2]=1.0;
			lhsm[ku2 + 2]=1.0;
			lhs	[ku3 + ku2 + 2]=1.0;
			lhsp[ku3 + ku2 + 2]=1.0;
			lhsm[ku3 + ku2 + 2]=1.0;
		}

		ks3 = k*ks3_const;
		ku4 = k*ku4_const;
		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three y-factors   
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue         
		 * ---------------------------------------------------------------------
		 */
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			for(int j=0; j<=grid_points[1]-1; j++){
        		ks2 = j*(IMAXP+1);
				ru1=c3c4*rho_i[ks3+ks2+i];
				cv[j]=vs[ks3+ks2+i];
				rhoq[j]=std::max(std::max(dy3+con43*ru1, dy5+c1c5*ru1), std::max(dymax+ru1, dy1));
			}
			for(int j=1; j<=grid_points[1]-2; j++){
        		ku3 = j*ku3_const;
				lhs[ku3+ku2]=0.0;
				lhs[ku3+ku2+1]=-dtty2*cv[j-1]-dtty1*rhoq[j-1];
				lhs[ku3+ku2+2]=1.0+c2dtty1*rhoq[j];
				lhs[ku3+ku2+3]=dtty2*cv[j+1]-dtty1*rhoq[j+1];
				lhs[ku3+ku2+4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                             
		 * ---------------------------------------------------------------------
		 */
		for(int i=1; i<=grid_points[0]-2; i++){
			ku2 = i*5;
			ku3 = 2*ku3_const;
			j=1;
			lhs[ku3_const+ku2+2]+=comz5;
			lhs[ku3_const+ku2+3]-=comz4;
			lhs[ku3_const+ku2+4]+=comz1;
			lhs[ku3+ku2+1]-=comz4;
			lhs[ku3+ku2+2]+=comz6;
			lhs[ku3+ku2+3]-=comz4;
			lhs[ku3+ku2+4]+=comz1;
		}
		for(int j=3; j<=grid_points[1]-4; j++){
      		ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				lhs[ku3+ku2+0]+=comz1;
				lhs[ku3+ku2+1]-=comz4;
				lhs[ku3+ku2+2]+=comz6;
				lhs[ku3+ku2+3]-=comz4;
				lhs[ku3+ku2+4]+=comz1;
			}
		}
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			j=grid_points[1]-3;
      		ku3 = j*ku3_const;
			lhs[ku3+ku2]+=comz1;
			lhs[ku3+ku2+1]-=comz4;
			lhs[ku3+ku2+2]+=comz6;
			lhs[ku3+ku2+3]-=comz4;
			lhs[ku3+ku3_const+ku2]+=comz1;
			lhs[ku3+ku3_const+ku2+1]-=comz4;
			lhs[ku3+ku3_const+ku2+2]+=comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, do the other two factors                    
		 * ---------------------------------------------------------------------
		 */
		for(int j=1; j<=grid_points[1]-2; j++){
			ks2 = j*(IMAXP+1);
			ku3 = j*ku3_const;
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				lhsp[ku3+ku2]=lhs[ku3+ku2];
				lhsp[ku3+ku2+1]=lhs[ku3+ku2+1]-dtty2*speed[ks3+ks2-(IMAXP+1)+i];
				lhsp[ku3+ku2+2]=lhs[ku3+ku2+2];
				lhsp[ku3+ku2+3]=lhs[ku3+ku2+3]+dtty2*speed[ks3+ks2+(IMAXP+1)+i];
				lhsp[ku3+ku2+4]=lhs[ku3+ku2+4];
				lhsm[ku3+ku2]=lhs[ku3+ku2];
				lhsm[ku3+ku2+1]=lhs[ku3+ku2+1]+dtty2*speed[ks3+ks2-(IMAXP+1)+i];
				lhsm[ku3+ku2+2]=lhs[ku3+ku2+2];
				lhsm[ku3+ku2+3]=lhs[ku3+ku2+3]-dtty2*speed[ks3+ks2+(IMAXP+1)+i];
				lhsm[ku3+ku2+4]=lhs[ku3+ku2+4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 */
		for(int j=0; j<=grid_points[1]-3; j++){
			ku3 = (j+1)*ku3_const;
			j1=j+1;
			j2=j+2;
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				fac1=1.0/lhs[ku3-ku3_const+ku2+2];
				lhs[ku3-ku3_const+ku2+3]*=fac1;
				lhs[ku3-ku3_const+ku2+4]*=fac1;
				for(int m=0; m<3; m++){
					rhs[ku4+ku3-ku3_const+ku2+m]*=fac1;
				}
				lhs[ku3+ku2+2]-=lhs[ku3+ku2+1]*lhs[ku3-ku3_const+ku2+3];
				lhs[ku3+ku2+3]-=lhs[ku3+ku2+1]*lhs[ku3-ku3_const+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+m]-=lhs[ku3+ku2+1]*rhs[ku4+ku3-ku3_const+ku2+m];
				}
				lhs[ku3+ku3_const+ku2+1]-=lhs[ku3+ku3_const+ku2]*lhs[ku3-ku3_const+ku2+3];
				lhs[ku3+ku3_const+ku2+2]-=lhs[ku3+ku3_const+ku2]*lhs[ku3-ku3_const+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku3_const+ku2+m]-=lhs[ku3+ku3_const+ku2]*rhs[ku4+ku3-ku3_const+ku2+m];
				}
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * the last two rows in this grid block are a bit different, 
		 * since they do not have two more rows available for the
		 * elimination of off-diagonal entries
		 * ---------------------------------------------------------------------
		 */
		j=grid_points[1]-2;
		j1=grid_points[1]-1;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			fac1=1.0/lhs[ku3+ku2+2];
			lhs[ku3+ku2+3]*=fac1;
			lhs[ku3+ku2+4]*=fac1;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+m]*=fac1;
			}
			lhs[ku3+ku3_const+ku2+2]-=lhs[ku3+ku3_const+ku2+1]*lhs[ku3+ku2+3];
			lhs[ku3+ku3_const+ku2+3]-=lhs[ku3+ku3_const+ku2+1]*lhs[ku3+ku2+4];
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku3_const+ku2+m]-=lhs[ku3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			fac2 = 1.0/lhs[ku3+ku3_const+ku2+2];
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku3_const+ku2+m]*=fac2;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors                 
		 * ---------------------------------------------------------------------
		 */
		for(int j=0; j<=grid_points[1]-3; j++){
			ku3 = (j+1)*ku3_const;
			j1=j+1;
			j2=j+2;
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				m=3;
				fac1=1.0/lhsp[ku3-ku3_const+ku2+2];
				lhsp[ku3-ku3_const+ku2+3]*=fac1;
				lhsp[ku3-ku3_const+ku2+4]*=fac1;
				rhs[ku4+ku3-ku3_const+ku2+m]*=fac1;
				lhsp[ku3+ku2+2]-=lhsp[ku3+ku2+1]*lhsp[ku3-ku3_const+ku2+3];
				lhsp[ku3+ku2+3]-=lhsp[ku3+ku2+1]*lhsp[ku3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku2+m]-=lhsp[ku3+ku2+1]*rhs[ku4+ku3-ku3_const+ku2+m];
				lhsp[ku3+ku3_const+ku2+1]-=lhsp[ku3+ku3_const+ku2]*lhsp[ku3-ku3_const+ku2+3];
				lhsp[ku3+ku3_const+ku2+2]-=lhsp[ku3+ku3_const+ku2]*lhsp[ku3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku3_const+ku2+m]-=lhsp[ku3+ku3_const+ku2]*rhs[ku4+ku3-ku3_const+ku2+m];
				m=4;
				fac1=1.0/lhsm[ku3-ku3_const+ku2+2];
				lhsm[ku3-ku3_const+ku2+3]*=fac1;
				lhsm[ku3-ku3_const+ku2+4]*=fac1;
				rhs[ku4+ku3-ku3_const+ku2+m]*=fac1;
				lhsm[ku3+ku2+2]-=lhsm[ku3+ku2+1]*lhsm[ku3-ku3_const+ku2+3];
				lhsm[ku3+ku2+3]-=lhsm[ku3+ku2+1]*lhsm[ku3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku2+m]-=lhsm[ku3+ku2+1]*rhs[ku4+ku3-ku3_const+ku2+m];
				lhsm[ku3+ku3_const+ku2+1]-=lhsm[ku3+ku3_const+ku2]*lhsm[ku3-ku3_const+ku2+3];
				lhsm[ku3+ku3_const+ku2+2]-=lhsm[ku3+ku3_const+ku2]*lhsm[ku3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku3_const+ku2+m]-=lhsm[ku3+ku3_const+ku2]*rhs[ku4+ku3-ku3_const+ku2+m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		j=grid_points[1]-2;
		j1=grid_points[1]-1;
		ku3 = j*ku3_const;
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			m=3;
			fac1=1.0/lhsp[ku3+ku2+2];
			lhsp[ku3+ku2+3]*=fac1;
			lhsp[ku3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsp[ku3+ku3_const+ku2+2]-=lhsp[ku3+ku3_const+ku2+1]*lhsp[ku3+ku2+3];
			lhsp[ku3+ku3_const+ku2+3]-=lhsp[ku3+ku3_const+ku2+1]*lhsp[ku3+ku2+4];
			rhs[ku4+ku3+ku3_const+ku2+m]-=lhsp[ku3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			m=4;
			fac1=1.0/lhsm[ku3+ku2+2];
			lhsm[ku3+ku2+3]*=fac1;
			lhsm[ku3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsm[ku3+ku3_const+ku2+2]-=lhsm[ku3+ku3_const+ku2+1]*lhsm[ku3+ku2+3];
			lhsm[ku3+ku3_const+ku2+3]-=lhsm[ku3+ku3_const+ku2+1]*lhsm[ku3+ku2+4];
			rhs[ku4+ku3+ku3_const+ku2+m]-=lhsm[ku3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			rhs[ku4+ku3+ku3_const+ku2+3]/=lhsp[ku3+ku3_const+ku2+2];
			rhs[ku4+ku3+ku3_const+ku2+4]/=lhsm[ku3+ku3_const+ku2+2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		j=grid_points[1]-2;
		j1=grid_points[1]-1;
		for(int i=1; i<=grid_points[0]-2; i++){
      		ku2 = i*5;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+m]-=lhs[ku3+ku2+3]*rhs[ku4+ku3+ku3_const+ku2+m];
			}
			rhs[ku4+ku3+ku2+3]-=lhsp[ku3+ku2+3]*rhs[ku4+ku3+ku3_const+ku2+3];
			rhs[ku4+ku3+ku2+4]-=lhsm[ku3+ku2+3]*rhs[ku4+ku3+ku3_const+ku2+4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(int j=grid_points[1]-3; j>=0; j--){
			ks3 = (j+1)*ku3_const;
			ku3 = j*ku3_const;
			j1=j+1;
			j2=j+2;
			for(int i=1; i<=grid_points[0]-2; i++){
        		ku2 = i*5;
				for(int m=0; m<3; m++){
					rhs[ku4+ks3-ku3_const+ku2+m]-=
						lhs[ku3+ku2+3]*rhs[ku4+ks3+ku2+m]+
						lhs[ku3+ku2+4]*rhs[ku4+ks3+ku3_const+ku2+m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[ku4+ks3-ku3_const+ku2+3]-=
					lhsp[ku3+ku2+3]*rhs[ku4+ks3+ku2+3]+
					lhsp[ku3+ku2+4]*rhs[ku4+ks3+ku3_const+ku2+3];
				rhs[ku4+ks3-ku3_const+ku2+4]-=
					lhsm[ku3+ku2+3]*rhs[ku4+ks3+ku2+4]+
					lhsm[ku3+ku2+4]*rhs[ku4+ks3+ku3_const+ku2+4];
			}
		}
	});
	if(timeron){timer_stop(T_YSOLVE);}
	pinvr();
}

/*
 * ---------------------------------------------------------------------
 * this function performs the solution of the approximate factorization
 * step in the z-direction for all five matrix components
 * simultaneously. The Thomas algorithm is employed to solve the
 * systems for the z-lines. Boundary conditions are non-periodic
 * ---------------------------------------------------------------------
 */
void z_solve(){
	if(timeron){timer_start(T_ZSOLVE);}
	std::for_each_n(policy, iter.front()+1, ny2, [&](int j){
		thread_local std::vector<double> cv(PROBLEM_SIZE);
		thread_local std::vector<double> rhos(PROBLEM_SIZE);
		thread_local std::vector<double> lhs((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsp((IMAXP+1)*(IMAXP+1)*5);
		thread_local std::vector<double> lhsm((IMAXP+1)*(IMAXP+1)*5);
		int ku4, ku3, ku2, kl3, kl2, ks3, ks2;
		int ku3_const=(IMAXP+1)*5;
		int ku4_const=(IMAXP+1)*ku3_const;
		int ks3_const=(IMAXP+1)*(IMAXP+1);
		int i, k, k1, k2, m;
		double ru1, fac1, fac2;
		
		for(int i=1; i<=nx2; i++){
			int ku2 = i*5;
			int ku3 = (ny2+1)*ku3_const;
			for(m=0; m<5; m++){
				lhs[ku2 + m]=0.0;
				lhsp[ku2 + m]=0.0;
				lhsm[ku2 + m]=0.0;
				lhs[ku3 + ku2 + m]=0.0;
				lhsp[ku3 + ku2 + m]=0.0;
				lhsm[ku3 + ku2 + m]=0.0;
			}
			lhs	[ku2 + 2]=1.0;
			lhsp[ku2 + 2]=1.0;
			lhsm[ku2 + 2]=1.0;
			lhs	[ku3 + ku2 + 2]=1.0;
			lhsp[ku3 + ku2 + 2]=1.0;
			lhsm[ku3 + ku2 + 2]=1.0;
		}

		ks2 = j*(IMAXP+1);
		ku3 = j*ku3_const;
		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three z-factors   
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue                          
		 * ---------------------------------------------------------------------
		 */
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int k=0; k<=nz2+1; k++){
        		ks3 = k*ks3_const;
				ru1=c3c4*rho_i[ks3+ks2+i];
				cv[k]=ws[ks3+ks2+i];
				rhos[k]=std::max(std::max(dz4+con43*ru1, dz5+c1c5*ru1), std::max(dzmax+ru1, dz1));
			}
			for(int k=1; k<=nz2; k++){
        		kl3 = k*ku3_const;
				lhs[kl3+ku2]=0.0;
				lhs[kl3+ku2+1]=-dttz2*cv[k-1]-dttz1*rhos[k-1];
				lhs[kl3+ku2+2]=1.0+c2dttz1*rhos[k];
				lhs[kl3+ku2+3]=dttz2*cv[k+1]-dttz1*rhos[k+1];
				lhs[kl3+ku2+4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                                  
		 * ---------------------------------------------------------------------
		 */
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			k=1;
			lhs[ku3_const+ku2+2]+=comz5;
			lhs[ku3_const+ku2+3]-=comz4;
			lhs[ku3_const+ku2+4]+=comz1;
			k=2;
			lhs[ku3_const+ku3_const+ku2+1]-=comz4;
			lhs[ku3_const+ku3_const+ku2+2]+=comz6;
			lhs[ku3_const+ku3_const+ku2+3]-=comz4;
			lhs[ku3_const+ku3_const+ku2+4]+=comz1;
		}
		for(int k=3; k<=nz2-2; k++){
      		kl3 = k*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				lhs[kl3+ku2]+=comz1;
				lhs[kl3+ku2+1]-=comz4;
				lhs[kl3+ku2+2]+=comz6;
				lhs[kl3+ku2+3]-=comz4;
				lhs[kl3+ku2+4]+=comz1;
			}
		}
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			k=nz2-1;
      		kl3 = k*ku3_const;
			lhs[kl3+ku2]+=comz1;
			lhs[kl3+ku2+1]-=comz4;
			lhs[kl3+ku2+2]+=comz6;
			lhs[kl3+ku2+3]-=comz4;
			k=nz2;
			lhs[kl3+ku3_const+ku2]+=comz1;
			lhs[kl3+ku3_const+ku2+1]-=comz4;
			lhs[kl3+ku3_const+ku2+2]+=comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, fill the other factors (u+c), (u-c) 
		 * ---------------------------------------------------------------------
		 */
		for(int k=1; k<=nz2; k++){
			ks3 = k*ks3_const;
			kl3 = k*ku3_const;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				lhsp[kl3+ku2]=lhs[kl3+ku2];
				lhsp[kl3+ku2+1]=lhs[kl3+ku2+1]-dttz2*speed[ks3-ks3_const+ks2+i];
				lhsp[kl3+ku2+2]=lhs[kl3+ku2+2];
				lhsp[kl3+ku2+3]=lhs[kl3+ku2+3]+dttz2*speed[ks3+ks3_const+ks2+i];
				lhsp[kl3+ku2+4]=lhs[kl3+ku2+4];
				lhsm[kl3+ku2]=lhs[kl3+ku2];
				lhsm[kl3+ku2+1]=lhs[kl3+ku2+1]+dttz2*speed[ks3-ks3_const+ks2+i];
				lhsm[kl3+ku2+2]=lhs[kl3+ku2+2];
				lhsm[kl3+ku2+3]=lhs[kl3+ku2+3]-dttz2*speed[ks3+ks3_const+ks2+i];
				lhsm[kl3+ku2+4]=lhs[kl3+ku2+4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 */
		for(int k=0; k<=grid_points[2]-3; k++){
			kl3 = (k+1)*ku3_const;
			ku4 = (k+1)*ku4_const;
			k1=k+1;
			k2=k+2;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				fac1=1.0/lhs[kl3-ku3_const+ku2+2];
				lhs[kl3-ku3_const+ku2+3]*=fac1;
				lhs[kl3-ku3_const+ku2+4]*=fac1;
				for(int m=0; m<3; m++){
					rhs[ku4-ku4_const+ku3+ku2+m]*=fac1;
				}
				lhs[kl3+ku2+2]-=lhs[kl3+ku2+1]*lhs[kl3-ku3_const+ku2+3];
				lhs[kl3+ku2+3]-=lhs[kl3+ku2+1]*lhs[kl3-ku3_const+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku3+ku2+m]-=lhs[kl3+ku2+1]*rhs[ku4-ku4_const+ku3+ku2+m];
				}
				lhs[kl3+ku3_const+ku2+1]-=lhs[kl3+ku3_const+ku2]*lhs[kl3-ku3_const+ku2+3];
				lhs[kl3+ku3_const+ku2+2]-=lhs[kl3+ku3_const+ku2]*lhs[kl3-ku3_const+ku2+4];
				for(int m=0; m<3; m++){
					rhs[ku4+ku4_const+ku3+ku2+m]-=lhs[kl3+ku3_const+ku2]*rhs[ku4-ku4_const+ku3+ku2+m];
				}
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * the last two rows in this grid block are a bit different, 
		 * since they do not have two more rows available for the
		 * elimination of off-diagonal entries
		 * ---------------------------------------------------------------------
		 */
		k=grid_points[2]-2;
		k1=grid_points[2]-1;
		ku4 = k*ku4_const;
		kl3 = k*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			fac1=1.0/lhs[kl3+ku2+2];
			lhs[kl3+ku2+3]*=fac1;
			lhs[kl3+ku2+4]*=fac1;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+m]*=fac1;
			}
			lhs[kl3+ku3_const+ku2+2]-=lhs[kl3+ku3_const+ku2+1]*lhs[kl3+ku2+3];
			lhs[kl3+ku3_const+ku2+3]-=lhs[kl3+ku3_const+ku2+1]*lhs[kl3+ku2+4];
			for(int m=0; m<3; m++){
				rhs[ku4+ku4_const+ku3+ku2+m]-=lhs[kl3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately
			 * ---------------------------------------------------------------------
			 */
			fac2=1.0/lhs[kl3+ku3_const+ku2+2];
			for(int m=0; m<3; m++){
				rhs[ku4+ku4_const+ku3+ku2+m]*=fac2;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors               
		 * ---------------------------------------------------------------------
		 */
		for(int k=0; k<=grid_points[2]-3; k++){
			ku4 = (k+1)*ku4_const;
			kl3 = (k+1)*ku3_const;
			k1=k+1;
			k2=k+2;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				m=3;
				fac1=1.0/lhsp[kl3-ku3_const+ku2+2];
				lhsp[kl3-ku3_const+ku2+3]*=fac1;
				lhsp[kl3-ku3_const+ku2+4]*=fac1;
				rhs[ku4-ku4_const+ku3+ku2+m]*=fac1;
				lhsp[kl3+ku2+2]-=lhsp[kl3+ku2+1]*lhsp[kl3-ku3_const+ku2+3];
				lhsp[kl3+ku2+3]-=lhsp[kl3+ku2+1]*lhsp[kl3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku2+m]-=lhsp[kl3+ku2+1]*rhs[ku4-ku4_const+ku3+ku2+m];
				lhsp[kl3+ku3_const+ku2+1]-=lhsp[kl3+ku3_const+ku2]*lhsp[kl3-ku3_const+ku2+3];
				lhsp[kl3+ku3_const+ku2+2]-=lhsp[kl3+ku3_const+ku2]*lhsp[kl3-ku3_const+ku2+4];
				rhs[ku4+ku4_const+ku3+ku2+m]-=lhsp[kl3+ku3_const+ku2]*rhs[ku4-ku4_const+ku3+ku2+m];
				m=4;
				fac1=1.0/lhsm[kl3-ku3_const+ku2+2];
				lhsm[kl3-ku3_const+ku2+3]*=fac1;
				lhsm[kl3-ku3_const+ku2+4]*=fac1;
				rhs[ku4-ku4_const+ku3+ku2+m]*=fac1;
				lhsm[kl3+ku2+2]-=lhsm[kl3+ku2+1]*lhsm[kl3-ku3_const+ku2+3];
				lhsm[kl3+ku2+3]-=lhsm[kl3+ku2+1]*lhsm[kl3-ku3_const+ku2+4];
				rhs[ku4+ku3+ku2+m]-=lhsm[kl3+ku2+1]*rhs[ku4-ku4_const+ku3+ku2+m];
				lhsm[kl3+ku3_const+ku2+1]-=lhsm[kl3+ku3_const+ku2]*lhsm[kl3-ku3_const+ku2+3];
				lhsm[kl3+ku3_const+ku2+2]-=lhsm[kl3+ku3_const+ku2]*lhsm[kl3-ku3_const+ku2+4];
				rhs[ku4+ku4_const+ku3+ku2+m]-=lhsm[kl3+ku3_const+ku2]*rhs[ku4-ku4_const+ku3+ku2+m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		k=grid_points[2]-2;
		k1=grid_points[2]-1;
		ku4 = k*ku4_const;
		kl3 = k*ku3_const;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			m=3;
			fac1=1.0/lhsp[kl3+ku2+2];
			lhsp[kl3+ku2+3]*=fac1;
			lhsp[kl3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsp[kl3+ku3_const+ku2+2]-=lhsp[kl3+ku3_const+ku2+1]*lhsp[kl3+ku2+3];
			lhsp[kl3+ku3_const+ku2+3]-=lhsp[kl3+ku3_const+ku2+1]*lhsp[kl3+ku2+4];
			rhs[ku4+ku4_const+ku3+ku2+m]-=lhsp[kl3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			m=4;
			fac1=1.0/lhsm[kl3+ku2+2];
			lhsm[kl3+ku2+3]*=fac1;
			lhsm[kl3+ku2+4]*=fac1;
			rhs[ku4+ku3+ku2+m]*=fac1;
			lhsm[kl3+ku3_const+ku2+2]-=lhsm[kl3+ku3_const+ku2+1]*lhsm[kl3+ku2+3];
			lhsm[kl3+ku3_const+ku2+3]-=lhsm[kl3+ku3_const+ku2+1]*lhsm[kl3+ku2+4];
			rhs[ku4+ku4_const+ku3+ku2+m]-=lhsm[kl3+ku3_const+ku2+1]*rhs[ku4+ku3+ku2+m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately (some of this is overkill
			 * if this is the last cell)
			 * ---------------------------------------------------------------------
			 */
			rhs[ku4+ku4_const+ku3+ku2+3]/=lhsp[kl3+ku3_const+ku2+2];
			rhs[ku4+ku4_const+ku3+ku2+4]/=lhsm[kl3+ku3_const+ku2+2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		k=grid_points[2]-2;
		k1=grid_points[2]-1;
		for(int i=1; i<=nx2; i++){
      		ku2 = i*5;
			for(int m=0; m<3; m++){
				rhs[ku4+ku3+ku2+m]-=lhs[kl3+ku2+3]*rhs[ku4+ku4_const+ku3+ku2+m];
			}
			rhs[ku4+ku3+ku2+3]-=lhsp[kl3+ku2+3]*rhs[ku4+ku4_const+ku3+ku2+3];
			rhs[ku4+ku3+ku2+4]-=lhsm[kl3+ku2+3]*rhs[ku4+ku4_const+ku3+ku2+4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * whether or not this is the last processor, we always have
		 * to complete the back-substitution 
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(int k=grid_points[2]-3; k>=0; k--){
			ku4 = (k+1)*ku4_const;
			kl3 = k*ku3_const;
			k1=k+1;
			k2=k+2;
			for(int i=1; i<=nx2; i++){
        		ku2 = i*5;
				for(int m = 0; m < 3; m++) {
					rhs[ku4-ku4_const+ku3+ku2+m]-=
						lhs[kl3+ku2+3]*rhs[ku4+ku3+ku2+m]+
						lhs[kl3+ku2+4]*rhs[ku4+ku4_const+ku3+ku2+m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[ku4-ku4_const+ku3+ku2+3]-=
					lhsp[kl3+ku2+3]*rhs[ku4+ku3+ku2+3]+
					lhsp[kl3+ku2+4]*rhs[ku4+ku4_const+ku3+ku2+3];
				rhs[ku4-ku4_const+ku3+ku2+4]-=
					lhsm[kl3+ku2+3]*rhs[ku4+ku3+ku2+4]+
					lhsm[kl3+ku2+4]*rhs[ku4+ku4_const+ku3+ku2+4];
			}
		}
	});
	if(timeron){timer_stop(T_ZSOLVE);}
	tzetar();
}
