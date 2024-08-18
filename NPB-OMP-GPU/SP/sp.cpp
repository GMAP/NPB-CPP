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

------------------------------------------------------------------------------

The OpenMP version is a parallel implementation of the serial C++ version
OpenMP version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-OMP

Authors of the OpenMP code:
	Júnior Löff <loffjh@gmail.com>
	
*/

#include "omp.h"
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
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
static double u[KMAX][JMAXP+1][IMAXP+1][5];
static double us[KMAX][JMAXP+1][IMAXP+1];
static double vs[KMAX][JMAXP+1][IMAXP+1];
static double ws[KMAX][JMAXP+1][IMAXP+1];
static double qs[KMAX][JMAXP+1][IMAXP+1];
static double rho_i[KMAX][JMAXP+1][IMAXP+1];
static double speed[KMAX][JMAXP+1][IMAXP+1];
static double square[KMAX][JMAXP+1][IMAXP+1];
static double rhs[KMAX][JMAXP+1][IMAXP+1][5];
static double forcing[KMAX][JMAXP+1][IMAXP+1][5];
static double cv[PROBLEM_SIZE];
static double rhon[PROBLEM_SIZE];
static double rhos[PROBLEM_SIZE];
static double rhoq[PROBLEM_SIZE];
static double cuf[PROBLEM_SIZE];
static double q[PROBLEM_SIZE];
static double ue[5][PROBLEM_SIZE];
static double buf[5][PROBLEM_SIZE];
static double lhs[IMAXP+1][IMAXP+1][5];
static double lhsp[IMAXP+1][IMAXP+1][5];
static double lhsm[IMAXP+1][IMAXP+1][5];
static double ce[13][5];
#else
static double (*u)[JMAXP+1][IMAXP+1][5]=(double(*)[JMAXP+1][IMAXP+1][5])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)*(5)));
static double (*us)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*vs)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*ws)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*qs)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*rho_i)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*speed)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*square)[JMAXP+1][IMAXP+1]=(double(*)[JMAXP+1][IMAXP+1])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)));
static double (*rhs)[JMAXP+1][IMAXP+1][5]=(double(*)[JMAXP+1][IMAXP+1][5])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)*(5)));
static double (*forcing)[JMAXP+1][IMAXP+1][5]=(double(*)[JMAXP+1][IMAXP+1][5])malloc(sizeof(double)*((KMAX)*(JMAXP+1)*(IMAXP+1)*(5)));
static double (*cv)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*rhon)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*rhos)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*rhoq)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*cuf)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*q)=(double*)malloc(sizeof(double)*(PROBLEM_SIZE));
static double (*ue)[PROBLEM_SIZE]=(double(*)[PROBLEM_SIZE])malloc(sizeof(double)*((PROBLEM_SIZE)*(5)));
static double (*buf)[PROBLEM_SIZE]=(double(*)[PROBLEM_SIZE])malloc(sizeof(double)*((PROBLEM_SIZE)*(5)));
static double (*lhs)[IMAXP+1][5]=(double(*)[IMAXP+1][5])malloc(sizeof(double)*((IMAXP+1)*(IMAXP+1)*(5)));
static double (*lhsp)[IMAXP+1][5]=(double(*)[IMAXP+1][5])malloc(sizeof(double)*((IMAXP+1)*(IMAXP+1)*(5)));
static double (*lhsm)[IMAXP+1][5]=(double(*)[IMAXP+1][5])malloc(sizeof(double)*((IMAXP+1)*(IMAXP+1)*(5)));
static double (*ce)[5]=(double(*)[5])malloc(sizeof(double)*((13)*(5)));
#endif
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
static int grid_points[3], nx2, ny2, nz2;
static boolean timeron;

/* function prototypes */
void add();
void adi();
void compute_rhs();
void error_norm(double rms[]);
void exact_rhs();
void exact_solution(double xi, double eta, double zeta, double dtemp[]);
void initialize();
void lhsinit(int ni, int nj);
void lhsinitj(int nj, int ni);
void ninvr();
void pinvr();
void rhs_norm(double rms[]);
void set_constants();
void txinvr();
void tzetar();
void verify(int no_time_steps, char* class_npb, boolean* verified);
void x_solve();
void y_solve();
void z_solve();

/* sp */
int main(int argc, char* argv[]){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	int i, niter, step, n3;
	double mflops, t, tmax, trecs[T_LAST+1];
	boolean verified;
	char class_npb;
	char* t_names[T_LAST+1];
	/*
	 * ---------------------------------------------------------------------
	 * read input file (if it exists), else take
	 * defaults from parameters
	 * ---------------------------------------------------------------------
	 */
	FILE* fp;
	if((fp=fopen("inputsp.data","r"))!=NULL){
		int result;
		printf(" Reading from input file inputsp.data\n");
		result=fscanf(fp,"%d", &niter);
		while(fgetc(fp)!='\n');
		result=fscanf(fp,"%lf",&dt);
		while(fgetc(fp)!='\n');
		result=fscanf(fp,"%d%d%d",&grid_points[0],&grid_points[1],&grid_points[2]);
		fclose(fp);
	}else{
		printf(" No input file inputsp.data. Using compiled defaults\n");
		niter=NITER_DEFAULT;
		dt=DT_DEFAULT;
		grid_points[0]=PROBLEM_SIZE;
		grid_points[1]=PROBLEM_SIZE;
		grid_points[2]=PROBLEM_SIZE;
	}
	if((fp=fopen("timer.flag","r"))!=NULL){
		timeron=TRUE;
		t_names[T_TOTAL]=(char*)"total";
		t_names[T_RHSX]=(char*)"rhsx";
		t_names[T_RHSY]=(char*)"rhsy";
		t_names[T_RHSZ]=(char*)"rhsz";
		t_names[T_RHS]=(char*)"rhs";
		t_names[T_XSOLVE]=(char*)"xsolve";
		t_names[T_YSOLVE]=(char*)"ysolve";
		t_names[T_ZSOLVE]=(char*)"zsolve";
		t_names[T_RDIS1]=(char*)"redist1";
		t_names[T_RDIS2]=(char*)"redist2";
		t_names[T_TZETAR]=(char*)"tzetar";
		t_names[T_NINVR]=(char*)"ninvr";
		t_names[T_PINVR]=(char*)"pinvr";
		t_names[T_TXINVR]=(char*)"txinvr";
		t_names[T_ADD]=(char*)"add";
		fclose(fp);
	}else{
		timeron = FALSE;
	}	
	printf("\n\n NAS Parallel Benchmarks 4.1 Parallel C++ version with OpenMP - SP Benchmark\n\n");
	printf(" Size: %4dx%4dx%4d\n",grid_points[0],grid_points[1],grid_points[2]);
	printf(" Iterations: %4d    dt: %10.6f\n",niter,dt);
	printf("\n");
	if((grid_points[0]>IMAX)||(grid_points[1]>JMAX)||(grid_points[2]>KMAX)){
		printf(" %d, %d, %d\n",grid_points[0],grid_points[1],grid_points[2]);
		printf(" Problem size too big for compiled array sizes\n");
		return 0;
	}
	nx2=grid_points[0]-2;
	ny2=grid_points[1]-2;
	nz2=grid_points[2]-2;
	set_constants();
	for(i=1;i<=T_LAST;i++){timer_clear(i);}
	exact_rhs();
	initialize();
	/*
	 * ---------------------------------------------------------------------
	 * do one time step to touch all code, and reinitialize
	 * ---------------------------------------------------------------------
	 */
	#pragma omp parallel
  	{
		adi();
	}
	initialize();
	for(i=1;i<=T_LAST;i++){timer_clear(i);}
	timer_start(1);
	#pragma omp parallel firstprivate(niter) private(step)
  	{
		for(step=1;step<=niter;step++){
			if((step%20)==0||step==1){
				#pragma omp master
					printf(" Time step %4d\n",step);
			}
			adi();
		}
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
	setenv("OMP_NUM_THREADS","1",0);
	c_print_results((char*)"SP",
			class_npb,
			grid_points[0],
			grid_points[1],
			grid_points[2],
			niter,
			tmax,
			mflops,
			(char*)"          floating point",
			verified,
			(char*)NPBVERSION,
			(char*)COMPILETIME,
			(char*)COMPILERVERSION,
			(char*)LIBVERSION,
			std::getenv("OMP_NUM_THREADS"),
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
		for(i=1; i<=T_LAST; i++){
			trecs[i]=timer_read(i);
		}
		if(tmax==0.0){tmax=1.0;}
		printf("  SECTION   Time (secs)\n");
		for(i=1; i<=T_LAST; i++){
			printf("  %-8s:%9.3f  (%6.2f%%)\n",t_names[i],trecs[i],trecs[i]*100./tmax);
			if(i==T_RHS){
				t=trecs[T_RHSX]+trecs[T_RHSY]+trecs[T_RHSZ];
				printf("    --> %8s:%9.3f  (%6.2f%%)\n","sub-rhs",t,t*100./tmax);
				t=trecs[T_RHS]-t;
				printf("    --> %8s:%9.3f  (%6.2f%%)\n","rest-rhs",t,t*100./tmax);
			}else if(i==T_ZSOLVE){
				t=trecs[T_ZSOLVE]-trecs[T_RDIS1]-trecs[T_RDIS2];
				printf("    --> %8s:%9.3f  (%6.2f%%)\n","sub-zsol",t,t*100./tmax);
			}else if(i==T_RDIS2){
				t=trecs[T_RDIS1]+trecs[T_RDIS2];
				printf("    --> %8s:%9.3f  (%6.2f%%)\n","redist",t,t*100./tmax);
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
	int i, j, k, m;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_ADD);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				for(m=0; m<5; m++){
					u[k][j][i][m]=u[k][j][i][m]+rhs[k][j][i][m];
				}
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_ADD);}
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
	int i, j, k, m;
	double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_RHS);}
	/*
	 * ---------------------------------------------------------------------
	 * compute the reciprocal of density, and the kinetic energy, 
	 * and the speed of sound. 
	 * ---------------------------------------------------------------------
	 */
	#pragma omp for
	for(k=0; k<=grid_points[2]-1; k++){
		for(j=0; j<=grid_points[1]-1; j++){
			for(i=0; i<=grid_points[0]-1; i++){
				rho_inv=1.0/u[k][j][i][0];
				rho_i[k][j][i]=rho_inv;
				us[k][j][i]=u[k][j][i][1]*rho_inv;
				vs[k][j][i]=u[k][j][i][2]*rho_inv;
				ws[k][j][i]=u[k][j][i][3]*rho_inv;
				square[k][j][i]=0.5*(
						u[k][j][i][1]*u[k][j][i][1]+ 
						u[k][j][i][2]*u[k][j][i][2]+
						u[k][j][i][3]*u[k][j][i][3])*rho_inv;
				qs[k][j][i]=square[k][j][i]*rho_inv;
				/*
				 * ---------------------------------------------------------------------
				 * (don't need speed and ainx until the lhs computation)
				 * ---------------------------------------------------------------------
				 */
				aux=c1c2*rho_inv*(u[k][j][i][4]-square[k][j][i]);
				speed[k][j][i]=sqrt(aux);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * copy the exact forcing term to the right hand side;  because 
	 * this forcing term is known, we can store it on the whole grid
	 * including the boundary                   
	 * ---------------------------------------------------------------------
	 */
	#pragma omp for
	for(k=0; k<=grid_points[2]-1; k++){
		for(j=0; j<=grid_points[1]-1; j++){
			for(i=0; i<=grid_points[0]-1; i++){
				for(m=0; m<5; m++){
					rhs[k][j][i][m]=forcing[k][j][i][m];
				}
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * compute xi-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron && thread_id==0){timer_start(T_RHSX);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				uijk=us[k][j][i];
				up1=us[k][j][i+1];
				um1=us[k][j][i-1];
				rhs[k][j][i][0]=rhs[k][j][i][0]+dx1tx1* 
					(u[k][j][i+1][0]-2.0*u[k][j][i][0]+u[k][j][i-1][0])-
					tx2*(u[k][j][i+1][1]-u[k][j][i-1][1]);
				rhs[k][j][i][1]=rhs[k][j][i][1]+dx2tx1* 
					(u[k][j][i+1][1]-2.0*u[k][j][i][1]+u[k][j][i-1][1])+
					xxcon2*con43*(up1-2.0*uijk+um1)-
					tx2*(u[k][j][i+1][1]*up1-u[k][j][i-1][1]*um1+
							(u[k][j][i+1][4]-square[k][j][i+1]-
							 u[k][j][i-1][4]+square[k][j][i-1])*c2);
				rhs[k][j][i][2]=rhs[k][j][i][2]+dx3tx1* 
					(u[k][j][i+1][2]-2.0*u[k][j][i][2]+u[k][j][i-1][2])+
					xxcon2*(vs[k][j][i+1]-2.0*vs[k][j][i]+vs[k][j][i-1])-
					tx2*(u[k][j][i+1][2]*up1-u[k][j][i-1][2]*um1);
				rhs[k][j][i][3]=rhs[k][j][i][3]+dx4tx1* 
					(u[k][j][i+1][3]-2.0*u[k][j][i][3]+u[k][j][i-1][3])+
					xxcon2*(ws[k][j][i+1]-2.0*ws[k][j][i]+ws[k][j][i-1])-
					tx2*(u[k][j][i+1][3]*up1-u[k][j][i-1][3]*um1);
				rhs[k][j][i][4]=rhs[k][j][i][4]+dx5tx1* 
					(u[k][j][i+1][4]-2.0*u[k][j][i][4]+u[k][j][i-1][4])+
					xxcon3*(qs[k][j][i+1]-2.0*qs[k][j][i]+qs[k][j][i-1])+
					xxcon4*(up1*up1-2.0*uijk*uijk+um1*um1)+
					xxcon5*(u[k][j][i+1][4]*rho_i[k][j][i+1]- 
							2.0*u[k][j][i][4]*rho_i[k][j][i]+
							u[k][j][i-1][4]*rho_i[k][j][i-1])-
					tx2*((c1*u[k][j][i+1][4]-c2*square[k][j][i+1])*up1-
							(c1*u[k][j][i-1][4]-c2*square[k][j][i-1])*um1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order xi-direction dissipation               
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			i=1;
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(5.0*u[k][j][i][m]-4.0*u[k][j][i+1][m]+u[k][j][i+2][m]);
			}
			i=2;
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(-4.0*u[k][j][i-1][m]+6.0*u[k][j][i][m]-
					 4.0*u[k][j][i+1][m]+u[k][j][i+2][m]);
			}
		}
		for(j=1; j<=ny2; j++){
			for(i=3; i<=nx2-2; i++){
				for(m=0; m<5; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
						(u[k][j][i-2][m]-4.0*u[k][j][i-1][m]+ 
						 6.0*u[k][j][i][m]-4.0*u[k][j][i+1][m]+ 
						 u[k][j][i+2][m]);
				}
			}
		}
		for(j=1; j<=ny2; j++){
			i=nx2-1;
			for(m=0; m<5; m++){
				rhs[k][j][i][m] = rhs[k][j][i][m]-dssp*
					(u[k][j][i-2][m]-4.0*u[k][j][i-1][m]+ 
					 6.0*u[k][j][i][m]-4.0*u[k][j][i+1][m]);
			}
			i=nx2;
			for(m=0; m<5; m++){
				rhs[k][j][i][m] = rhs[k][j][i][m]-dssp*
					(u[k][j][i-2][m]-4.0*u[k][j][i-1][m]+5.0*u[k][j][i][m]);
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_RHSX);}
	/*
	 * ---------------------------------------------------------------------
	 * compute eta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron && thread_id==0){timer_start(T_RHSY);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				vijk=vs[k][j][i];
				vp1=vs[k][j+1][i];
				vm1=vs[k][j-1][i];
				rhs[k][j][i][0]=rhs[k][j][i][0]+dy1ty1* 
					(u[k][j+1][i][0]-2.0*u[k][j][i][0]+u[k][j-1][i][0])-
					ty2*(u[k][j+1][i][2]-u[k][j-1][i][2]);
				rhs[k][j][i][1]=rhs[k][j][i][1]+dy2ty1* 
					(u[k][j+1][i][1]-2.0*u[k][j][i][1]+u[k][j-1][i][1])+
					yycon2*(us[k][j+1][i]-2.0*us[k][j][i]+us[k][j-1][i])-
					ty2*(u[k][j+1][i][1]*vp1-u[k][j-1][i][1]*vm1);
				rhs[k][j][i][2]=rhs[k][j][i][2]+dy3ty1* 
					(u[k][j+1][i][2]-2.0*u[k][j][i][2]+u[k][j-1][i][2])+
					yycon2*con43*(vp1-2.0*vijk+vm1)-
					ty2*(u[k][j+1][i][2]*vp1-u[k][j-1][i][2]*vm1+
							(u[k][j+1][i][4]-square[k][j+1][i]- 
							 u[k][j-1][i][4]+square[k][j-1][i])* c2);
				rhs[k][j][i][3]=rhs[k][j][i][3]+dy4ty1* 
					(u[k][j+1][i][3]-2.0*u[k][j][i][3]+u[k][j-1][i][3])+
					yycon2*(ws[k][j+1][i]-2.0*ws[k][j][i]+ws[k][j-1][i])-
					ty2*(u[k][j+1][i][3]*vp1-u[k][j-1][i][3]*vm1);
				rhs[k][j][i][4]=rhs[k][j][i][4]+dy5ty1* 
					(u[k][j+1][i][4]-2.0*u[k][j][i][4]+u[k][j-1][i][4])+
					yycon3*(qs[k][j+1][i]-2.0*qs[k][j][i]+qs[k][j-1][i])+
					yycon4*(vp1*vp1- 2.0*vijk*vijk + vm1*vm1)+
					yycon5*(u[k][j+1][i][4]*rho_i[k][j+1][i]- 
							2.0*u[k][j][i][4]*rho_i[k][j][i]+
							u[k][j-1][i][4]*rho_i[k][j-1][i])-
					ty2*((c1*u[k][j+1][i][4]-c2*square[k][j+1][i])*vp1 -
							(c1*u[k][j-1][i][4]-c2*square[k][j-1][i])*vm1);
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order eta-direction dissipation         
		 * ---------------------------------------------------------------------
		 */
		j = 1;
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(5.0*u[k][j][i][m]-4.0*u[k][j+1][i][m]+u[k][j+2][i][m]);
			}
		}
		j = 2;
		for(i=1; i<=nx2; i++){
			for(m = 0; m < 5; m++) {
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(-4.0*u[k][j-1][i][m]+6.0*u[k][j][i][m]-
					 4.0*u[k][j+1][i][m]+u[k][j+2][i][m]);
			}
		}
		for (j=3; j<=ny2-2; j++){
			for(i=1; i<=nx2; i++){
				for(m=0; m<5; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
						(u[k][j-2][i][m]-4.0*u[k][j-1][i][m]+ 
						 6.0*u[k][j][i][m]-4.0*u[k][j+1][i][m]+ 
						 u[k][j+2][i][m]);
				}
			}
		}
		j=ny2-1;
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp*
					(u[k][j-2][i][m]-4.0*u[k][j-1][i][m]+ 
					 6.0*u[k][j][i][m]-4.0*u[k][j+1][i][m]);
			}
		}
		j=ny2;
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp*
					(u[k][j-2][i][m]-4.0*u[k][j-1][i][m]+5.0*u[k][j][i][m]);
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_RHSY);}
	/*
	 * ---------------------------------------------------------------------
	 * compute zeta-direction fluxes 
	 * ---------------------------------------------------------------------
	 */
	if(timeron && thread_id==0){timer_start(T_RHSZ);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				wijk=ws[k][j][i];
				wp1=ws[k+1][j][i];
				wm1=ws[k-1][j][i];
				rhs[k][j][i][0]=rhs[k][j][i][0]+dz1tz1* 
					(u[k+1][j][i][0]-2.0*u[k][j][i][0]+u[k-1][j][i][0])-
					tz2*(u[k+1][j][i][3]-u[k-1][j][i][3]);
				rhs[k][j][i][1]=rhs[k][j][i][1]+dz2tz1* 
					(u[k+1][j][i][1]-2.0*u[k][j][i][1]+u[k-1][j][i][1])+
					zzcon2*(us[k+1][j][i]-2.0*us[k][j][i]+us[k-1][j][i])-
					tz2*(u[k+1][j][i][1]*wp1-u[k-1][j][i][1]*wm1);
				rhs[k][j][i][2]=rhs[k][j][i][2]+dz3tz1* 
					(u[k+1][j][i][2]-2.0*u[k][j][i][2]+u[k-1][j][i][2])+
					zzcon2*(vs[k+1][j][i]-2.0*vs[k][j][i]+vs[k-1][j][i])-
					tz2*(u[k+1][j][i][2]*wp1-u[k-1][j][i][2]*wm1);
				rhs[k][j][i][3]=rhs[k][j][i][3]+dz4tz1* 
					(u[k+1][j][i][3]-2.0*u[k][j][i][3]+u[k-1][j][i][3])+
					zzcon2*con43*(wp1-2.0*wijk+wm1)-
					tz2*(u[k+1][j][i][3]*wp1-u[k-1][j][i][3]*wm1+
							(u[k+1][j][i][4]-square[k+1][j][i]- 
							 u[k-1][j][i][4]+square[k-1][j][i])*c2);
				rhs[k][j][i][4]=rhs[k][j][i][4]+dz5tz1* 
					(u[k+1][j][i][4]-2.0*u[k][j][i][4]+u[k-1][j][i][4])+
					zzcon3*(qs[k+1][j][i]-2.0*qs[k][j][i]+qs[k-1][j][i])+
					zzcon4*(wp1*wp1-2.0*wijk*wijk+wm1*wm1)+
					zzcon5*(u[k+1][j][i][4]*rho_i[k+1][j][i]- 
							2.0*u[k][j][i][4]*rho_i[k][j][i]+
							u[k-1][j][i][4]*rho_i[k-1][j][i])-
					tz2*((c1*u[k+1][j][i][4]-c2*square[k+1][j][i])*wp1-
							(c1*u[k-1][j][i][4]-c2*square[k-1][j][i])*wm1);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * add fourth order zeta-direction dissipation                
	 * ---------------------------------------------------------------------
	 */
	k=1;
	#pragma omp for
	for(j=1; j<=ny2; j++){
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(5.0*u[k][j][i][m]-4.0*u[k+1][j][i][m]+u[k+2][j][i][m]);
			}
		}
	}
	k=2;
	#pragma omp for
	for(j=1; j<=ny2; j++){
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
					(-4.0*u[k-1][j][i][m]+6.0*u[k][j][i][m]-
					 4.0*u[k+1][j][i][m]+u[k+2][j][i][m]);
			}
		}
	}
	#pragma omp for
	for(k=3; k<=nz2-2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				for(m=0; m<5; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]-dssp* 
						(u[k-2][j][i][m]-4.0*u[k-1][j][i][m]+ 
						 6.0*u[k][j][i][m]-4.0*u[k+1][j][i][m]+ 
						 u[k+2][j][i][m]);
				}
			}
		}
	}
	k=nz2-1;
	#pragma omp for
	for(j=1; j<=ny2; j++){
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp*
					(u[k-2][j][i][m]-4.0*u[k-1][j][i][m]+ 
					 6.0*u[k][j][i][m]-4.0*u[k+1][j][i][m]);
			}
		}
	}
	k=nz2;
	#pragma omp for
	for(j=1; j<=ny2; j++){
		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-dssp*
					(u[k-2][j][i][m]-4.0*u[k-1][j][i][m]+5.0*u[k][j][i][m]);
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_RHSZ);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				for(m=0; m<5; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]*dt;
				}
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_RHS);}
}

/*
 * ---------------------------------------------------------------------
 * this function computes the norm of the difference between the
 * computed solution and the exact solution
 * ---------------------------------------------------------------------
 */
void error_norm(double rms[]){
	int i, j, k, m, d;
	double xi, eta, zeta, u_exact[5], add;
	for(m=0; m<5; m++){rms[m]=0.0;}
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(j=0; j<=grid_points[1]-1; j++){
			eta=(double)j*dnym1;
			for(i=0; i<=grid_points[0]-1; i++){
				xi=(double)i*dnxm1;
				exact_solution(xi, eta, zeta, u_exact);
				for(m=0; m<5; m++){
					add=u[k][j][i][m]-u_exact[m];
					rms[m]=rms[m]+add*add;
				}
			}
		}
	}
	for(m=0; m<5; m++){
		for(d=0; d<3; d++){
			rms[m]=rms[m]/(double)(grid_points[d]-2);
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
	double dtemp[5], xi, eta, zeta, dtpp;
	int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
	/*
	 * ---------------------------------------------------------------------
	 * initialize                                  
	 * ---------------------------------------------------------------------
	 */
	for(k=0; k<=grid_points[2]-1; k++){
		for(j=0; j<= grid_points[1]-1; j++){
			for(i=0; i<=grid_points[0]-1; i++){
				for(m=0; m<5; m++){
					forcing[k][j][i][m]=0.0;
				}
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * xi-direction flux differences                      
	 * ---------------------------------------------------------------------
	 */
	for(k=1; k<=grid_points[2]-2; k++){
		zeta=(double)k*dnzm1;
		for(j=1; j<=grid_points[1]-2; j++){
			eta=(double)j*dnym1;
			for(i=0; i<=grid_points[0]-1; i++){
				xi=(double)i*dnxm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(m=0; m<5; m++){
					ue[m][i]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(m=1; m<5; m++){
					buf[m][i]=dtpp*dtemp[m];
				}
				cuf[i]=buf[1][i]*buf[1][i];
				buf[0][i]=cuf[i]+buf[2][i]*buf[2][i]+buf[3][i]*buf[3][i]; 
				q[i]=0.5*(buf[1][i]*ue[1][i]+buf[2][i]*ue[2][i]+
						buf[3][i]*ue[3][i]);
			}
			for(i=1; i<=grid_points[0]-2; i++){
				im1=i-1;
				ip1=i+1;
				forcing[k][j][i][0]=forcing[k][j][i][0]-
					tx2*(ue[1][ip1]-ue[1][im1])+
					dx1tx1*(ue[0][ip1]-2.0*ue[0][i]+ue[0][im1]);
				forcing[k][j][i][1]=forcing[k][j][i][1]-tx2*(
						(ue[1][ip1]*buf[1][ip1]+c2*(ue[4][ip1]-q[ip1]))-
						(ue[1][im1]*buf[1][im1]+c2*(ue[4][im1]-q[im1])))+
					xxcon1*(buf[1][ip1]-2.0*buf[1][i]+buf[1][im1])+
					dx2tx1*(ue[1][ip1]-2.0*ue[1][i]+ue[1][im1]);
				forcing[k][j][i][2]=forcing[k][j][i][2]-tx2*(
						ue[2][ip1]*buf[1][ip1]-ue[2][im1]*buf[1][im1])+
					xxcon2*(buf[2][ip1]-2.0*buf[2][i]+buf[2][im1])+
					dx3tx1*(ue[2][ip1]-2.0*ue[2][i]+ue[2][im1]);
				forcing[k][j][i][3]=forcing[k][j][i][3]-tx2*(
						ue[3][ip1]*buf[1][ip1]-ue[3][im1]*buf[1][im1])+
					xxcon2*(buf[3][ip1]-2.0*buf[3][i]+buf[3][im1])+
					dx4tx1*(ue[3][ip1]-2.0*ue[3][i]+ue[3][im1]);
				forcing[k][j][i][4]=forcing[k][j][i][4]-tx2*(
						buf[1][ip1]*(c1*ue[4][ip1]-c2*q[ip1])-
						buf[1][im1]*(c1*ue[4][im1]-c2*q[im1]))+
					0.5*xxcon3*(buf[0][ip1]-2.0*buf[0][i]+buf[0][im1])+
					xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
					xxcon5*(buf[4][ip1]-2.0*buf[4][i]+buf[4][im1])+
					dx5tx1*(ue[4][ip1]-2.0*ue[4][i]+ue[4][im1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                         
			 * ---------------------------------------------------------------------
			 */
			for(m=0; m<5; m++){
				i=1;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(5.0*ue[m][i]-4.0*ue[m][i+1]+ue[m][i+2]);
				i=2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(-4.0*ue[m][i-1]+6.0*ue[m][i]-
					 4.0*ue[m][i+1]+ue[m][i+2]);
			}
			for(m=0; m<5; m++){
				for(i=3; i<=grid_points[0]-4; i++){				
					forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
						(ue[m][i-2]-4.0*ue[m][i-1]+
						 6.0*ue[m][i]-4.0*ue[m][i+1]+ue[m][i+2]);
				}
			}
			for(m=0; m<5; m++){
				i=grid_points[0]-3;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][i-2]-4.0*ue[m][i-1]+
					 6.0*ue[m][i]-4.0*ue[m][i+1]);
				i=grid_points[0]-2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][i-2]-4.0*ue[m][i-1]+5.0*ue[m][i]);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * eta-direction flux differences             
	 * ---------------------------------------------------------------------
	 */
	for(k=1; k<=grid_points[2]-2; k++){
		zeta=(double)k*dnzm1;
		for(i=1; i<=grid_points[0]-2; i++){
			xi=(double)i*dnxm1;
			for(j=0;j<=grid_points[1]-1;j++){
				eta=(double)j*dnym1;
				exact_solution(xi, eta, zeta, dtemp);
				for(m=0; m<5; m++){
					ue[m][j]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(m=1; m<5; m++){
					buf[m][j]=dtpp*dtemp[m];
				}
				cuf[j]=buf[2][j]*buf[2][j];
				buf[0][j]=cuf[j]+buf[1][j]*buf[1][j]+buf[3][j]*buf[3][j];
				q[j]=0.5*(buf[1][j]*ue[1][j]+buf[2][j]*ue[2][j]+
						buf[3][j]*ue[3][j]);
			}
			for(j=1; j<=grid_points[1]-2; j++){
				jm1=j-1;
				jp1=j+1;
				forcing[k][j][i][0]=forcing[k][j][i][0]-
					ty2*(ue[2][jp1]-ue[2][jm1])+
					dy1ty1*(ue[0][jp1]-2.0*ue[0][j]+ue[0][jm1]);
				forcing[k][j][i][1]=forcing[k][j][i][1]-ty2*(
						ue[1][jp1]*buf[2][jp1]-ue[1][jm1]*buf[2][jm1])+
					yycon2*(buf[1][jp1]-2.0*buf[1][j]+buf[1][jm1])+
					dy2ty1*(ue[1][jp1]-2.0*ue[1][j]+ue[1][jm1]);
				forcing[k][j][i][2]=forcing[k][j][i][2]-ty2*(
						(ue[2][jp1]*buf[2][jp1]+c2*(ue[4][jp1]-q[jp1]))-
						(ue[2][jm1]*buf[2][jm1]+c2*(ue[4][jm1]-q[jm1])))+
					yycon1*(buf[2][jp1]-2.0*buf[2][j]+buf[2][jm1])+
					dy3ty1*(ue[2][jp1]-2.0*ue[2][j]+ue[2][jm1]);
				forcing[k][j][i][3]=forcing[k][j][i][3]-ty2*(
						ue[3][jp1]*buf[2][jp1]-ue[3][jm1]*buf[2][jm1])+
					yycon2*(buf[3][jp1]-2.0*buf[3][j]+buf[3][jm1])+
					dy4ty1*(ue[3][jp1]-2.0*ue[3][j]+ue[3][jm1]);
				forcing[k][j][i][4]=forcing[k][j][i][4]-ty2*(
						buf[2][jp1]*(c1*ue[4][jp1]-c2*q[jp1])-
						buf[2][jm1]*(c1*ue[4][jm1]-c2*q[jm1]))+
					0.5*yycon3*(buf[0][jp1]-2.0*buf[0][j]+
							buf[0][jm1])+
					yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
					yycon5*(buf[4][jp1]-2.0*buf[4][j]+buf[4][jm1])+
					dy5ty1*(ue[4][jp1]-2.0*ue[4][j]+ue[4][jm1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation                      
			 * ---------------------------------------------------------------------
			 */
			for(m=0; m<5; m++){
				j=1;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(5.0*ue[m][j]-4.0*ue[m][j+1]+ue[m][j+2]);
				j=2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(-4.0*ue[m][j-1]+6.0*ue[m][j]-
					 4.0*ue[m][j+1]+ue[m][j+2]);
			}
			for(m=0; m<5; m++){
				for(j=3; j<=grid_points[1]-4; j++){				
					forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
						(ue[m][j-2]-4.0*ue[m][j-1]+
						 6.0*ue[m][j]-4.0*ue[m][j+1]+ue[m][j+2]);
				}
			}
			for(m=0; m<5; m++){
				j=grid_points[1]-3;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][j-2]-4.0*ue[m][j-1]+
					 6.0*ue[m][j]-4.0*ue[m][j+1]);
				j=grid_points[1]-2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][j-2]-4.0*ue[m][j-1]+5.0*ue[m][j]);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * zeta-direction flux differences                      
	 * ---------------------------------------------------------------------
	 */
	for(j=1; j<=grid_points[1]-2; j++){
		eta=(double)j*dnym1;
		for(i=1; i<=grid_points[0]-2; i++){
			xi=(double)i*dnxm1;
			for(k=0; k<=grid_points[2]-1; k++){
				zeta=(double)k*dnzm1;
				exact_solution(xi, eta, zeta, dtemp);
				for(m=0; m<5; m++){
					ue[m][k]=dtemp[m];
				}
				dtpp=1.0/dtemp[0];
				for(m=1; m<5; m++){
					buf[m][k]=dtpp*dtemp[m];
				}
				cuf[k]=buf[3][k]*buf[3][k];
				buf[0][k]=cuf[k]+buf[1][k]*buf[1][k]+buf[2][k]*buf[2][k];
				q[k]=0.5*(buf[1][k]*ue[1][k]+buf[2][k]*ue[2][k]+
						buf[3][k]*ue[3][k]);
			}
			for(k=1; k<=grid_points[2]-2; k++){
				km1=k-1;
				kp1=k+1;
				forcing[k][j][i][0]=forcing[k][j][i][0]-
					tz2*(ue[3][kp1]-ue[3][km1])+
					dz1tz1*(ue[0][kp1]-2.0*ue[0][k]+ue[0][km1]);
				forcing[k][j][i][1]=forcing[k][j][i][1]-tz2*(
						ue[1][kp1]*buf[3][kp1]-ue[1][km1]*buf[3][km1])+
					zzcon2*(buf[1][kp1]-2.0*buf[1][k]+buf[1][km1])+
					dz2tz1*(ue[1][kp1]-2.0*ue[1][k]+ue[1][km1]);
				forcing[k][j][i][2]=forcing[k][j][i][2]-tz2*(
						ue[2][kp1]*buf[3][kp1]-ue[2][km1]*buf[3][km1])+
					zzcon2*(buf[2][kp1]-2.0*buf[2][k]+buf[2][km1])+
					dz3tz1*(ue[2][kp1]-2.0*ue[2][k]+ue[2][km1]);
				forcing[k][j][i][3]=forcing[k][j][i][3]-tz2*(
						(ue[3][kp1]*buf[3][kp1]+c2*(ue[4][kp1]-q[kp1]))-
						(ue[3][km1]*buf[3][km1]+c2*(ue[4][km1]-q[km1])))+
					zzcon1*(buf[3][kp1]-2.0*buf[3][k]+buf[3][km1])+
					dz4tz1*(ue[3][kp1]-2.0*ue[3][k]+ue[3][km1]);
				forcing[k][j][i][4]=forcing[k][j][i][4]-tz2*(
						buf[3][kp1]*(c1*ue[4][kp1]-c2*q[kp1])-
						buf[3][km1]*(c1*ue[4][km1]-c2*q[km1]))+
					0.5*zzcon3*(buf[0][kp1]-2.0*buf[0][k]+buf[0][km1])+
					zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
					zzcon5*(buf[4][kp1]-2.0*buf[4][k]+buf[4][km1])+
					dz5tz1*(ue[4][kp1]-2.0*ue[4][k]+ue[4][km1]);
			}
			/*
			 * ---------------------------------------------------------------------
			 * fourth-order dissipation
			 * ---------------------------------------------------------------------
			 */
			for(m=0; m<5; m++){
				k=1;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(5.0*ue[m][k]-4.0*ue[m][k+1]+ue[m][k+2]);
				k=2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(-4.0*ue[m][k-1]+6.0*ue[m][k]-
					 4.0*ue[m][k+1]+ue[m][k+2]);
			}
			for(m=0; m<5; m++){
				for(k=3; k<=grid_points[2]-4; k++){				
					forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
						(ue[m][k-2]-4.0*ue[m][k-1]+
						 6.0*ue[m][k]-4.0*ue[m][k+1]+ue[m][k+2]);
				}
			}
			for(m=0; m<5; m++){
				k=grid_points[2]-3;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][k-2]-4.0*ue[m][k-1]+
					 6.0*ue[m][k]-4.0*ue[m][k+1]);
				k=grid_points[2]-2;
				forcing[k][j][i][m]=forcing[k][j][i][m]-dssp*
					(ue[m][k-2]-4.0*ue[m][k-1]+5.0*ue[m][k]);
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * now change the sign of the forcing function
	 * ---------------------------------------------------------------------
	 */
	for(k=1; k<=grid_points[2]-2; k++){
		for(j=1; j<=grid_points[1]-2; j++){
			for(i=1; i<=grid_points[0]-2; i++){
				for(m=0; m<5; m++){
					forcing[k][j][i][m]=-1.0*forcing[k][j][i][m];
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
void exact_solution(double xi, double eta, double zeta, double dtemp[]){
	int m;
	for(m=0; m<5; m++){
		dtemp[m]=ce[0][m]+xi*
			(ce[1][m]+xi*
			 (ce[4][m]+xi*
			  (ce[7][m]+xi*
			   ce[10][m])))+eta*
			(ce[2][m]+eta*
			 (ce[5][m]+eta*
			  (ce[8][m]+eta*
			   ce[11][m])))+zeta*
			(ce[3][m]+zeta*
			 (ce[6][m]+zeta*
			  (ce[9][m]+zeta*
			   ce[12][m])));
	}
}

/*
 * ---------------------------------------------------------------------
 * this subroutine initializes the field variable u using 
 * tri-linear transfinite interpolation of the boundary values     
 * ---------------------------------------------------------------------
 */
void initialize(){
	int i, j, k, m, ix, iy, iz;
	double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];
	/*
	 * ---------------------------------------------------------------------
	 * later (in compute_rhs) we compute 1/u for every element. a few of 
	 * the corner elements are not used, but it convenient (and faster) 
	 * to compute the whole thing with a simple loop. make sure those 
	 * values are nonzero by initializing the whole thing here. 
	 * ---------------------------------------------------------------------
	 */
	for(k=0; k<=grid_points[2]-1; k++){
		for(j=0; j<=grid_points[1]-1; j++){
			for(i=0; i<=grid_points[0]-1; i++){
				u[k][j][i][0]=1.0;
				u[k][j][i][1]=0.0;
				u[k][j][i][2]=0.0;
				u[k][j][i][3]=0.0;
				u[k][j][i][4]=1.0;
			}
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * first store the "interpolated" values everywhere on the grid    
	 * ---------------------------------------------------------------------
	 */
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(j=0; j<=grid_points[1]-1; j++){
			eta=(double)j*dnym1;
			for(i=0; i<=grid_points[0]-1; i++){
				xi=(double)i*dnxm1;
				for(ix=0; ix<2; ix++){
					Pxi=(double)ix;
					exact_solution(Pxi, eta, zeta, &Pface[ix][0][0]);
				}
				for(iy=0; iy<2; iy++){
					Peta=(double)iy;
					exact_solution(xi, Peta, zeta, &Pface[iy][1][0]);
				}
				for(iz=0; iz<2; iz++){
					Pzeta=(double)iz;
					exact_solution(xi, eta, Pzeta, &Pface[iz][2][0]);
				}
				for(m=0; m<5; m++){
					Pxi=xi*Pface[1][0][m]+(1.0-xi)*Pface[0][0][m];
					Peta=eta*Pface[1][1][m]+(1.0-eta)*Pface[0][1][m];
					Pzeta=zeta*Pface[1][2][m]+(1.0-zeta)*Pface[0][2][m];
					u[k][j][i][m]=Pxi+Peta+Pzeta- 
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
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(j=0; j<=grid_points[1]-1; j++){
			eta=(double)j*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
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
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(j=0; j<=grid_points[1]-1; j++){
			eta=(double)j*dnym1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
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
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(i=0; i<=grid_points[0]-1; i++){
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
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
	for(k=0; k<=grid_points[2]-1; k++){
		zeta=(double)k*dnzm1;
		for(i=0; i<=grid_points[0]-1; i++){
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
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
	for(j=0; j<=grid_points[1]-1; j++){
		eta=(double)j*dnym1;
		for(i=0; i<=grid_points[0]-1; i++){
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
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
	for(j=0; j<=grid_points[1]-1; j++){
		eta=(double)j*dnym1;
		for(i=0; i<=grid_points[0]-1; i++){
			xi=(double)i*dnxm1;
			exact_solution(xi, eta, zeta, temp);
			for(m=0; m<5; m++){
				u[k][j][i][m]=temp[m];
			}
		}
	}
}

void lhsinit(int ni, int nj){
	int j, m;
	/*
	 * ---------------------------------------------------------------------
	 * zap the whole left hand side for starters
	 * set all diagonal values to 1. This is overkill, but convenient
	 * ---------------------------------------------------------------------
	 */
	for(j=1; j<=nj; j++){
		for(m=0; m<5; m++){
			lhs[j][0][m]=0.0;
			lhsp[j][0][m]=0.0;
			lhsm[j][0][m]=0.0;
			lhs[j][ni][m]=0.0;
			lhsp[j][ni][m]=0.0;
			lhsm[j][ni][m]=0.0;
		}
		lhs[j][0][2]=1.0;
		lhsp[j][0][2]=1.0;
		lhsm[j][0][2]=1.0;
		lhs[j][ni][2]=1.0;
		lhsp[j][ni][2]=1.0;
		lhsm[j][ni][2]=1.0;
	}
}

void lhsinitj(int nj, int ni){
	int i, m;
	/*
	 * ---------------------------------------------------------------------
	 * zap the whole left hand side for starters
	 * set all diagonal values to 1. This is overkill, but convenient
	 * ---------------------------------------------------------------------
	 */
	for(i=1; i<=ni; i++){
		for(m=0; m<5; m++){
			lhs[0][i][m]=0.0;
			lhsp[0][i][m]=0.0;
			lhsm[0][i][m]=0.0;
			lhs[nj][i][m]=0.0;
			lhsp[nj][i][m]=0.0;
			lhsm[nj][i][m]=0.0;
		}
		lhs[0][i][2]=1.0;
		lhsp[0][i][2]=1.0;
		lhsm[0][i][2]=1.0;
		lhs[nj][i][2]=1.0;
		lhsp[nj][i][2]=1.0;
		lhsm[nj][i][2]=1.0;
	}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication              
 * ---------------------------------------------------------------------
 */
void ninvr(){
	int i, j, k;
	double r1, r2, r3, r4, r5, t1, t2;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_NINVR);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				r1=rhs[k][j][i][0];
				r2=rhs[k][j][i][1];
				r3=rhs[k][j][i][2];
				r4=rhs[k][j][i][3];
				r5=rhs[k][j][i][4];
				t1=bt*r3;
				t2=0.5*(r4+r5);
				rhs[k][j][i][0]=-r2;
				rhs[k][j][i][1]=r1;
				rhs[k][j][i][2]=bt*(r4-r5);
				rhs[k][j][i][3]=-t1+t2;
				rhs[k][j][i][4]=t1+t2;
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_NINVR);}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication                       
 * ---------------------------------------------------------------------
 */
void pinvr(){
	int i, j, k;
	double r1, r2, r3, r4, r5, t1, t2;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_PINVR);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				r1=rhs[k][j][i][0];
				r2=rhs[k][j][i][1];
				r3=rhs[k][j][i][2];
				r4=rhs[k][j][i][3];
				r5=rhs[k][j][i][4];
				t1=bt*r1;
				t2=0.5*(r4+r5);
				rhs[k][j][i][0]=bt*(r4-r5);
				rhs[k][j][i][1]=-r3;
				rhs[k][j][i][2]=r2;
				rhs[k][j][i][3]=-t1+t2;
				rhs[k][j][i][4]=t1+t2;
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_PINVR);}
}

void rhs_norm(double rms[]){
	int i, j, k, d, m;
	double add;
	for(m=0;m<5;m++){rms[m]=0.0;}
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				for(m=0; m<5; m++){
					add=rhs[k][j][i][m];
					rms[m]=rms[m]+add*add;
				} 
			} 
		} 
	} 
	for(m=0; m<5; m++){
		for(d=0; d<3; d++){
			rms[m]=rms[m]/(double)(grid_points[d]-2);
		}
		rms[m]=sqrt(rms[m]);
	}
}

void set_constants(){
	ce[0][0]=2.0;
	ce[1][0]=0.0;
	ce[2][0]=0.0;
	ce[3][0]=4.0;
	ce[4][0]=5.0;
	ce[5][0]=3.0;
	ce[6][0]=0.5;
	ce[7][0]=0.02;
	ce[8][0]=0.01;
	ce[9][0]=0.03;
	ce[10][0]=0.5;
	ce[11][0]=0.4;
	ce[12][0]=0.3;
	/* */
	ce[0][1]=1.0;
	ce[1][1]=0.0;
	ce[2][1]=0.0;
	ce[3][1]=0.0;
	ce[4][1]=1.0;
	ce[5][1]=2.0;
	ce[6][1]=3.0;
	ce[7][1]=0.01;
	ce[8][1]=0.03;
	ce[9][1]=0.02;
	ce[10][1]=0.4;
	ce[11][1]=0.3;
	ce[12][1]=0.5;
	/* */
	ce[0][2]=2.0;
	ce[1][2]=2.0;
	ce[2][2]=0.0;
	ce[3][2]=0.0;
	ce[4][2]=0.0;
	ce[5][2]=2.0;
	ce[6][2]=3.0;
	ce[7][2]=0.04;
	ce[8][2]=0.03;
	ce[9][2]=0.05;
	ce[10][2]=0.3;
	ce[11][2]=0.5;
	ce[12][2]=0.4;
	/* */
	ce[0][3]=2.0;
	ce[1][3]=2.0;
	ce[2][3]=0.0;
	ce[3][3]=0.0;
	ce[4][3]=0.0;
	ce[5][3]=2.0;
	ce[6][3]=3.0;
	ce[7][3]=0.03;
	ce[8][3]=0.05;
	ce[9][3]=0.04;
	ce[10][3]=0.2;
	ce[11][3]=0.1;
	ce[12][3]=0.3;
	/* */
	ce[0][4]=5.0;
	ce[1][4]=4.0;
	ce[2][4]=3.0;
	ce[3][4]=2.0;
	ce[4][4]=0.1;
	ce[5][4]=0.4;
	ce[6][4]=0.3;
	ce[7][4]=0.05;
	ce[8][4]=0.04;
	ce[9][4]=0.03;
	ce[10][4]=0.1;
	ce[11][4]=0.3;
	ce[12][4]=0.2;
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
	dxmax=max(dx3, dx4);
	dymax=max(dy2, dy4);
	dzmax=max(dz2, dz3);
	/* */
	dssp=0.25*max(dx1, max(dy1, dz1));
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
	int i, j, k;
	double t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, r4, r5, ac2inv;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_TXINVR);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				ru1=rho_i[k][j][i];
				uu=us[k][j][i];
				vv=vs[k][j][i];
				ww=ws[k][j][i];
				ac=speed[k][j][i];
				ac2inv=ac*ac;
				r1=rhs[k][j][i][0];
				r2=rhs[k][j][i][1];
				r3=rhs[k][j][i][2];
				r4=rhs[k][j][i][3];
				r5=rhs[k][j][i][4];
				t1=c2/ac2inv*(qs[k][j][i]*r1-uu*r2-vv*r3-ww*r4+r5);
				t2=bt*ru1*(uu*r1-r2);
				t3=(bt*ru1*ac)*t1;
				rhs[k][j][i][0]=r1-t1;
				rhs[k][j][i][1]=-ru1*(ww*r1-r4);
				rhs[k][j][i][2]=ru1*(vv*r1-r3);
				rhs[k][j][i][3]=-t2+t3;
				rhs[k][j][i][4]=t2+t3;
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_TXINVR);}
}

/*
 * ---------------------------------------------------------------------
 * block-diagonal matrix-vector multiplication                       
 * ---------------------------------------------------------------------
 */
void tzetar(){
	int i, j, k;
	double t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3, r4, r5, btuz, ac2u, uzik1;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_TZETAR);}
	#pragma omp for
	for(k=1; k<=nz2; k++){
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				xvel=us[k][j][i];
				yvel=vs[k][j][i];
				zvel=ws[k][j][i];
				ac=speed[k][j][i];
				ac2u=ac*ac;
				r1=rhs[k][j][i][0];
				r2=rhs[k][j][i][1];
				r3=rhs[k][j][i][2];
				r4=rhs[k][j][i][3];
				r5=rhs[k][j][i][4];
				uzik1=u[k][j][i][0];
				btuz=bt*uzik1;
				t1=btuz/ac*(r4+r5);
				t2=r3+t1;
				t3=btuz*(r4-r5);
				rhs[k][j][i][0]=t2;
				rhs[k][j][i][1]=-uzik1*r2+xvel*t2;
				rhs[k][j][i][2]=uzik1*r1+yvel*t2;
				rhs[k][j][i][3]=zvel*t2+t3;
				rhs[k][j][i][4]=uzik1*(-xvel*r2+yvel*r1) + 
					qs[k][j][i]*t2+c2iv*ac2u*t1+zvel*t3;
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_TZETAR);}
}

/*
 * ---------------------------------------------------------------------
 * verification routine                         
 * ---------------------------------------------------------------------
 */
void verify(int no_time_steps, char* class_npb, boolean* verified){
	double xcrref[5], xceref[5], xcrdif[5], xcedif[5], epsilon, xce[5], xcr[5], dtref;
	int m;
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
	for(m=0;m<5;m++){xcr[m]=xcr[m]/dt;}
	*class_npb='U';
	*verified=TRUE;
	for(m=0;m<5;m++){xcrref[m]=1.0;xceref[m]=1.0;}
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
	for(m=0; m<5; m++){
		xcrdif[m]=fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m]=fabs((xce[m]-xceref[m])/xceref[m]);
	}
	/*
	 * ---------------------------------------------------------------------
	 * output the comparison of computed results to known cases
	 * ---------------------------------------------------------------------
	 */
	if(*class_npb!='U'){
		printf(" Verification being performed for class %c\n",*class_npb);
		printf(" accuracy setting for epsilon = %20.13E\n",epsilon);
		*verified=(fabs(dt-dtref)<=epsilon);
		if(!(*verified)){  
			*class_npb='U';
			printf(" DT does not match the reference value of %15.8E\n",dtref);
		} 
	}else{
		printf(" Unknown class\n");
	}
	if(*class_npb!='U'){
		printf(" Comparison of RMS-norms of residual\n");
	}else{
		printf(" RMS-norms of residual\n");
	}
	for(m=0;m<5;m++){
		if(*class_npb=='U'){
			printf("          %2d%20.13E\n",m+1,xcr[m]);
		}else if(xcrdif[m]<=epsilon){
			printf("          %2d%20.13E%20.13E%20.13E\n",m+1,xcr[m],xcrref[m],xcrdif[m]);
		}else {
			*verified=FALSE;
			printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n",m+1,xcr[m],xcrref[m],xcrdif[m]);
		}
	}
	if(*class_npb!='U'){
		printf(" Comparison of RMS-norms of solution error\n");
	}else{
		printf(" RMS-norms of solution error\n");
	}
	for(m=0;m<5;m++){
		if(*class_npb=='U'){
			printf("          %2d%20.13E\n",m+1,xce[m]);
		}else if(xcedif[m]<=epsilon){
			printf("          %2d%20.13E%20.13E%20.13E\n",m+1,xce[m],xceref[m],xcedif[m]);
		}else{
			*verified = FALSE;
			printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n",m+1,xce[m],xceref[m],xcedif[m]);
		}
	}
	if(*class_npb=='U'){
		printf(" No reference values provided\n");
		printf(" No verification performed\n");
	}else if(*verified){
		printf(" Verification Successful\n");
	}else{
		printf(" Verification failed\n");
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
	int i, j, k, i1, i2, m;
	double ru1, fac1, fac2;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_XSOLVE);}

	#pragma omp for
	for(k=1; k<=nz2; k++){
		double cv[PROBLEM_SIZE], rhon[PROBLEM_SIZE];
		double lhs[IMAXP+1][IMAXP+1][5];
		double lhsp[IMAXP+1][IMAXP+1][5];
		double lhsm[IMAXP+1][IMAXP+1][5];

		for(j=1; j<=ny2; j++){
			for(m=0; m<5; m++){
				lhs[j][0][m]=0.0;
				lhsp[j][0][m]=0.0;
				lhsm[j][0][m]=0.0;
				lhs[j][nx2+1][m]=0.0;
				lhsp[j][nx2+1][m]=0.0;
				lhsm[j][nx2+1][m]=0.0;
			}
			lhs[j][0][2]=1.0;
			lhsp[j][0][2]=1.0;
			lhsm[j][0][2]=1.0;
			lhs[j][nx2+1][2]=1.0;
			lhsp[j][nx2+1][2]=1.0;
			lhsm[j][nx2+1][2]=1.0;
		}

		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three x-factors  
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue                   
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			for(i=0; i<=grid_points[0]-1; i++){
				ru1=c3c4*rho_i[k][j][i];
				cv[i]=us[k][j][i];
				rhon[i]=max(max(dx2+con43*ru1,dx5+c1c5*ru1), max(dxmax+ru1,dx1));
			}
			for(i=1; i<=nx2; i++){
				lhs[j][i][0]=0.0;
				lhs[j][i][1]=-dttx2*cv[i-1]-dttx1*rhon[i-1];
				lhs[j][i][2]=1.0+c2dttx1*rhon[i];
				lhs[j][i][3]=dttx2*cv[i+1]-dttx1*rhon[i+1];
				lhs[j][i][4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                             
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			i=1;
			lhs[j][i][2]=lhs[j][i][2]+comz5;
			lhs[j][i][3]=lhs[j][i][3]-comz4;
			lhs[j][i][4]=lhs[j][i][4]+comz1;
			lhs[j][i+1][1]=lhs[j][i+1][1]-comz4;
			lhs[j][i+1][2]=lhs[j][i+1][2]+comz6;
			lhs[j][i+1][3]=lhs[j][i+1][3]-comz4;
			lhs[j][i+1][4]=lhs[j][i+1][4]+comz1;
		}
		for(j=1; j<=ny2; j++){
			for(i=3; i<=grid_points[0]-4; i++){
				lhs[j][i][0]=lhs[j][i][0]+comz1;
				lhs[j][i][1]=lhs[j][i][1]-comz4;
				lhs[j][i][2]=lhs[j][i][2]+comz6;
				lhs[j][i][3]=lhs[j][i][3]-comz4;
				lhs[j][i][4]=lhs[j][i][4]+comz1;
			}
		}
		for(j=1; j<=ny2; j++){
			i=grid_points[0]-3;
			lhs[j][i][0]=lhs[j][i][0]+comz1;
			lhs[j][i][1]=lhs[j][i][1]-comz4;
			lhs[j][i][2]=lhs[j][i][2]+comz6;
			lhs[j][i][3]=lhs[j][i][3]-comz4;
			lhs[j][i+1][0]=lhs[j][i+1][0]+comz1;
			lhs[j][i+1][1]=lhs[j][i+1][1]-comz4;
			lhs[j][i+1][2]=lhs[j][i+1][2]+comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, fill the other factors (u+c), (u-c) by adding to 
		 * the first  
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			for(i=1; i<=nx2; i++){
				lhsp[j][i][0]=lhs[j][i][0];
				lhsp[j][i][1]=lhs[j][i][1]-dttx2*speed[k][j][i-1];
				lhsp[j][i][2]=lhs[j][i][2];
				lhsp[j][i][3]=lhs[j][i][3]+dttx2*speed[k][j][i+1];
				lhsp[j][i][4]=lhs[j][i][4];
				lhsm[j][i][0]=lhs[j][i][0];
				lhsm[j][i][1]=lhs[j][i][1]+dttx2*speed[k][j][i-1];
				lhsm[j][i][2]=lhs[j][i][2];
				lhsm[j][i][3]=lhs[j][i][3]-dttx2*speed[k][j][i+1];
				lhsm[j][i][4]=lhs[j][i][4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 * perform the thomas algorithm; first, FORWARD ELIMINATION     
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			for(i=0; i<=grid_points[0]-3; i++){
				i1=i+1;
				i2=i+2;
				fac1=1.0/lhs[j][i][2];
				lhs[j][i][3]=fac1*lhs[j][i][3];
				lhs[j][i][4]=fac1*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				}
				lhs[j][i1][2]=lhs[j][i1][2]-lhs[j][i1][1]*lhs[j][i][3];
				lhs[j][i1][3]=lhs[j][i1][3]-lhs[j][i1][1]*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhs[j][i1][1]*rhs[k][j][i][m];
				}
				lhs[j][i2][1]=lhs[j][i2][1]-lhs[j][i2][0]*lhs[j][i][3];
				lhs[j][i2][2]=lhs[j][i2][2]-lhs[j][i2][0]*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j][i2][m]=rhs[k][j][i2][m]-lhs[j][i2][0]*rhs[k][j][i][m];
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
		for(j=1; j<=ny2; j++){
			i=grid_points[0]-2;
			i1=grid_points[0]-1;
			fac1=1.0/lhs[j][i][2];
			lhs[j][i][3]=fac1*lhs[j][i][3];
			lhs[j][i][4]=fac1*lhs[j][i][4];
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			}
			lhs[j][i1][2]=lhs[j][i1][2]-lhs[j][i1][1]*lhs[j][i][3];
			lhs[j][i1][3]=lhs[j][i1][3]-lhs[j][i1][1]*lhs[j][i][4];
			for(m=0; m<3; m++){
				rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhs[j][i1][1]*rhs[k][j][i][m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			fac2 = 1.0/lhs[j][i1][2];
			for(m=0; m<3; m++){
				rhs[k][j][i1][m]=fac2*rhs[k][j][i1][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors                 
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			for(i=0; i<=grid_points[0]-3; i++){
				i1=i+1;
				i2=i+2;
				m=3;
				fac1=1.0/lhsp[j][i][2];
				lhsp[j][i][3]=fac1*lhsp[j][i][3];
				lhsp[j][i][4]=fac1*lhsp[j][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsp[j][i1][2]=lhsp[j][i1][2]-lhsp[j][i1][1]*lhsp[j][i][3];
				lhsp[j][i1][3]=lhsp[j][i1][3]-lhsp[j][i1][1]*lhsp[j][i][4];
				rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhsp[j][i1][1]*rhs[k][j][i][m];
				lhsp[j][i2][1]=lhsp[j][i2][1]-lhsp[j][i2][0]*lhsp[j][i][3];
				lhsp[j][i2][2]=lhsp[j][i2][2]-lhsp[j][i2][0]*lhsp[j][i][4];
				rhs[k][j][i2][m]=rhs[k][j][i2][m]-lhsp[j][i2][0]*rhs[k][j][i][m];
				m=4;
				fac1=1.0/lhsm[j][i][2];
				lhsm[j][i][3]=fac1*lhsm[j][i][3];
				lhsm[j][i][4]=fac1*lhsm[j][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsm[j][i1][2]=lhsm[j][i1][2]-lhsm[j][i1][1]*lhsm[j][i][3];
				lhsm[j][i1][3]=lhsm[j][i1][3]-lhsm[j][i1][1]*lhsm[j][i][4];
				rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhsm[j][i1][1]*rhs[k][j][i][m];
				lhsm[j][i2][1]=lhsm[j][i2][1]-lhsm[j][i2][0]*lhsm[j][i][3];
				lhsm[j][i2][2]=lhsm[j][i2][2]-lhsm[j][i2][0]*lhsm[j][i][4];
				rhs[k][j][i2][m]=rhs[k][j][i2][m]-lhsm[j][i2][0]*rhs[k][j][i][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			i=grid_points[0]-2;
			i1=grid_points[0]-1;
			m=3;
			fac1=1.0/lhsp[j][i][2];
			lhsp[j][i][3]=fac1*lhsp[j][i][3];
			lhsp[j][i][4]=fac1*lhsp[j][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsp[j][i1][2]=lhsp[j][i1][2]-lhsp[j][i1][1]*lhsp[j][i][3];
			lhsp[j][i1][3]=lhsp[j][i1][3]-lhsp[j][i1][1]*lhsp[j][i][4];
			rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhsp[j][i1][1]*rhs[k][j][i][m];
			m=4;
			fac1=1.0/lhsm[j][i][2];
			lhsm[j][i][3]=fac1*lhsm[j][i][3];
			lhsm[j][i][4]=fac1*lhsm[j][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsm[j][i1][2]=lhsm[j][i1][2]-lhsm[j][i1][1]*lhsm[j][i][3];
			lhsm[j][i1][3]=lhsm[j][i1][3]-lhsm[j][i1][1]*lhsm[j][i][4];
			rhs[k][j][i1][m]=rhs[k][j][i1][m]-lhsm[j][i1][1]*rhs[k][j][i][m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately
			 * ---------------------------------------------------------------------
			 */
			rhs[k][j][i1][3]=rhs[k][j][i1][3]/lhsp[j][i1][2];
			rhs[k][j][i1][4]=rhs[k][j][i1][4]/lhsm[j][i1][2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			i=grid_points[0]-2;
			i1=grid_points[0]-1;
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-lhs[j][i][3]*rhs[k][j][i1][m];
			}
			rhs[k][j][i][3]=rhs[k][j][i][3]-lhsp[j][i][3]*rhs[k][j][i1][3];
			rhs[k][j][i][4]=rhs[k][j][i][4]-lhsm[j][i][3]*rhs[k][j][i1][4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=ny2; j++){
			for(i=grid_points[0]-3; i>=0; i--){
				i1=i+1;
				i2=i+2;
				for(m=0; m<3; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]- 
						lhs[j][i][3]*rhs[k][j][i1][m]-
						lhs[j][i][4]*rhs[k][j][i2][m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[k][j][i][3]=rhs[k][j][i][3]- 
					lhsp[j][i][3]*rhs[k][j][i1][3] -
					lhsp[j][i][4]*rhs[k][j][i2][3];
				rhs[k][j][i][4]=rhs[k][j][i][4]- 
					lhsm[j][i][3]*rhs[k][j][i1][4]-
					lhsm[j][i][4]*rhs[k][j][i2][4];
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_XSOLVE);}
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
	int i, j, k, j1, j2, m;
	double ru1, fac1, fac2;
	int thread_id = omp_get_thread_num();

	if(timeron && thread_id==0){timer_start(T_YSOLVE);}
	#pragma omp for
	for(k=1; k<=grid_points[2]-2; k++){
		double cv[PROBLEM_SIZE], rhoq[PROBLEM_SIZE];
		double lhs[IMAXP+1][IMAXP+1][5];
		double lhsp[IMAXP+1][IMAXP+1][5];
		double lhsm[IMAXP+1][IMAXP+1][5];

		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				lhs[0][i][m]=0.0;
				lhsp[0][i][m]=0.0;
				lhsm[0][i][m]=0.0;
				lhs[ny2+1][i][m]=0.0;
				lhsp[ny2+1][i][m]=0.0;
				lhsm[ny2+1][i][m]=0.0;
			}
			lhs[0][i][2]=1.0;
			lhsp[0][i][2]=1.0;
			lhsm[0][i][2]=1.0;
			lhs[ny2+1][i][2]=1.0;
			lhsp[ny2+1][i][2]=1.0;
			lhsm[ny2+1][i][2]=1.0;
		}

		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three y-factors   
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue         
		 * ---------------------------------------------------------------------
		 */
		for(i=1; i<=grid_points[0]-2; i++){
			for(j=0; j<=grid_points[1]-1; j++){
				ru1=c3c4*rho_i[k][j][i];
				cv[j]=vs[k][j][i];
				rhoq[j]=max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
			}
			for(j=1; j<=grid_points[1]-2; j++){
				lhs[j][i][0]=0.0;
				lhs[j][i][1]=-dtty2*cv[j-1]-dtty1*rhoq[j-1];
				lhs[j][i][2]=1.0+c2dtty1*rhoq[j];
				lhs[j][i][3]=dtty2*cv[j+1]-dtty1*rhoq[j+1];
				lhs[j][i][4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                             
		 * ---------------------------------------------------------------------
		 */
		for(i=1; i<=grid_points[0]-2; i++){
			j=1;
			lhs[j][i][2]=lhs[j][i][2]+comz5;
			lhs[j][i][3]=lhs[j][i][3]-comz4;
			lhs[j][i][4]=lhs[j][i][4]+comz1;
			lhs[j+1][i][1]=lhs[j+1][i][1]-comz4;
			lhs[j+1][i][2]=lhs[j+1][i][2]+comz6;
			lhs[j+1][i][3]=lhs[j+1][i][3]-comz4;
			lhs[j+1][i][4]=lhs[j+1][i][4]+comz1;
		}
		for(j=3; j<=grid_points[1]-4; j++){
			for(i=1; i<=grid_points[0]-2; i++){
				lhs[j][i][0]=lhs[j][i][0]+comz1;
				lhs[j][i][1]=lhs[j][i][1]-comz4;
				lhs[j][i][2]=lhs[j][i][2]+comz6;
				lhs[j][i][3]=lhs[j][i][3]-comz4;
				lhs[j][i][4]=lhs[j][i][4]+comz1;
			}
		}
		for(i=1; i<=grid_points[0]-2; i++){
			j=grid_points[1]-3;
			lhs[j][i][0]=lhs[j][i][0]+comz1;
			lhs[j][i][1]=lhs[j][i][1]-comz4;
			lhs[j][i][2]=lhs[j][i][2]+comz6;
			lhs[j][i][3]=lhs[j][i][3]-comz4;
			lhs[j+1][i][0]=lhs[j+1][i][0]+comz1;
			lhs[j+1][i][1]=lhs[j+1][i][1]-comz4;
			lhs[j+1][i][2]=lhs[j+1][i][2]+comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, do the other two factors                    
		 * ---------------------------------------------------------------------
		 */
		for(j=1; j<=grid_points[1]-2; j++){
			for(i=1; i<=grid_points[0]-2; i++){
				lhsp[j][i][0]=lhs[j][i][0];
				lhsp[j][i][1]=lhs[j][i][1]-dtty2*speed[k][j-1][i];
				lhsp[j][i][2]=lhs[j][i][2];
				lhsp[j][i][3]=lhs[j][i][3]+dtty2*speed[k][j+1][i];
				lhsp[j][i][4]=lhs[j][i][4];
				lhsm[j][i][0]=lhs[j][i][0];
				lhsm[j][i][1]=lhs[j][i][1]+dtty2*speed[k][j-1][i];
				lhsm[j][i][2]=lhs[j][i][2];
				lhsm[j][i][3]=lhs[j][i][3]-dtty2*speed[k][j+1][i];
				lhsm[j][i][4]=lhs[j][i][4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 */
		for(j=0; j<=grid_points[1]-3; j++){
			j1=j+1;
			j2=j+2;
			for(i=1; i<=grid_points[0]-2; i++){
				fac1=1.0/lhs[j][i][2];
				lhs[j][i][3]=fac1*lhs[j][i][3];
				lhs[j][i][4]=fac1*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				}
				lhs[j1][i][2]=lhs[j1][i][2]-lhs[j1][i][1]*lhs[j][i][3];
				lhs[j1][i][3]=lhs[j1][i][3]-lhs[j1][i][1]*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhs[j1][i][1]*rhs[k][j][i][m];
				}
				lhs[j2][i][1]=lhs[j2][i][1]-lhs[j2][i][0]*lhs[j][i][3];
				lhs[j2][i][2]=lhs[j2][i][2]-lhs[j2][i][0]*lhs[j][i][4];
				for(m=0; m<3; m++){
					rhs[k][j2][i][m]=rhs[k][j2][i][m]-lhs[j2][i][0]*rhs[k][j][i][m];
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
		for(i=1; i<=grid_points[0]-2; i++){
			fac1=1.0/lhs[j][i][2];
			lhs[j][i][3]=fac1*lhs[j][i][3];
			lhs[j][i][4]=fac1*lhs[j][i][4];
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			}
			lhs[j1][i][2]=lhs[j1][i][2]-lhs[j1][i][1]*lhs[j][i][3];
			lhs[j1][i][3]=lhs[j1][i][3]-lhs[j1][i][1]*lhs[j][i][4];
			for(m=0; m<3; m++){
				rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhs[j1][i][1]*rhs[k][j][i][m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			fac2 = 1.0/lhs[j1][i][2];
			for(m=0; m<3; m++){
				rhs[k][j1][i][m]=fac2*rhs[k][j1][i][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors                 
		 * ---------------------------------------------------------------------
		 */
		for(j=0; j<=grid_points[1]-3; j++){
			j1=j+1;
			j2=j+2;
			for(i=1; i<=grid_points[0]-2; i++){
				m=3;
				fac1=1.0/lhsp[j][i][2];
				lhsp[j][i][3]=fac1*lhsp[j][i][3];
				lhsp[j][i][4]=fac1*lhsp[j][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsp[j1][i][2]=lhsp[j1][i][2]-lhsp[j1][i][1]*lhsp[j][i][3];
				lhsp[j1][i][3]=lhsp[j1][i][3]-lhsp[j1][i][1]*lhsp[j][i][4];
				rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhsp[j1][i][1]*rhs[k][j][i][m];
				lhsp[j2][i][1]=lhsp[j2][i][1]-lhsp[j2][i][0]*lhsp[j][i][3];
				lhsp[j2][i][2]=lhsp[j2][i][2]-lhsp[j2][i][0]*lhsp[j][i][4];
				rhs[k][j2][i][m]=rhs[k][j2][i][m]-lhsp[j2][i][0]*rhs[k][j][i][m];
				m=4;
				fac1=1.0/lhsm[j][i][2];
				lhsm[j][i][3]=fac1*lhsm[j][i][3];
				lhsm[j][i][4]=fac1*lhsm[j][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsm[j1][i][2]=lhsm[j1][i][2]-lhsm[j1][i][1]*lhsm[j][i][3];
				lhsm[j1][i][3]=lhsm[j1][i][3]-lhsm[j1][i][1]*lhsm[j][i][4];
				rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhsm[j1][i][1]*rhs[k][j][i][m];
				lhsm[j2][i][1]=lhsm[j2][i][1]-lhsm[j2][i][0]*lhsm[j][i][3];
				lhsm[j2][i][2]=lhsm[j2][i][2]-lhsm[j2][i][0]*lhsm[j][i][4];
				rhs[k][j2][i][m]=rhs[k][j2][i][m]-lhsm[j2][i][0]*rhs[k][j][i][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		j=grid_points[1]-2;
		j1=grid_points[1]-1;
		for(i=1; i<=grid_points[0]-2; i++){
			m=3;
			fac1=1.0/lhsp[j][i][2];
			lhsp[j][i][3]=fac1*lhsp[j][i][3];
			lhsp[j][i][4]=fac1*lhsp[j][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsp[j1][i][2]=lhsp[j1][i][2]-lhsp[j1][i][1]*lhsp[j][i][3];
			lhsp[j1][i][3]=lhsp[j1][i][3]-lhsp[j1][i][1]*lhsp[j][i][4];
			rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhsp[j1][i][1]*rhs[k][j][i][m];
			m=4;
			fac1=1.0/lhsm[j][i][2];
			lhsm[j][i][3]=fac1*lhsm[j][i][3];
			lhsm[j][i][4]=fac1*lhsm[j][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsm[j1][i][2]=lhsm[j1][i][2]-lhsm[j1][i][1]*lhsm[j][i][3];
			lhsm[j1][i][3]=lhsm[j1][i][3]-lhsm[j1][i][1]*lhsm[j][i][4];
			rhs[k][j1][i][m]=rhs[k][j1][i][m]-lhsm[j1][i][1]*rhs[k][j][i][m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately 
			 * ---------------------------------------------------------------------
			 */
			rhs[k][j1][i][3]=rhs[k][j1][i][3]/lhsp[j1][i][2];
			rhs[k][j1][i][4]=rhs[k][j1][i][4]/lhsm[j1][i][2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		j=grid_points[1]-2;
		j1=grid_points[1]-1;
		for(i=1; i<=grid_points[0]-2; i++){
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-lhs[j][i][3]*rhs[k][j1][i][m];
			}
			rhs[k][j][i][3]=rhs[k][j][i][3]-lhsp[j][i][3]*rhs[k][j1][i][3];
			rhs[k][j][i][4]=rhs[k][j][i][4]-lhsm[j][i][3]*rhs[k][j1][i][4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(j=grid_points[1]-3; j>=0; j--){
			j1=j+1;
			j2=j+2;
			for(i=1; i<=grid_points[0]-2; i++){
				for(m=0; m<3; m++){
					rhs[k][j][i][m]=rhs[k][j][i][m]- 
						lhs[j][i][3]*rhs[k][j1][i][m]-
						lhs[j][i][4]*rhs[k][j2][i][m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[k][j][i][3]=rhs[k][j][i][3]- 
					lhsp[j][i][3]*rhs[k][j1][i][3]-
					lhsp[j][i][4]*rhs[k][j2][i][3];
				rhs[k][j][i][4]=rhs[k][j][i][4]- 
					lhsm[j][i][3]*rhs[k][j1][i][4]-
					lhsm[j][i][4]*rhs[k][j2][i][4];
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_YSOLVE);}
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
	int i, j, k, k1, k2, m;
	double ru1, fac1, fac2;
	int thread_id = omp_get_thread_num();
	
	if(timeron && thread_id==0){timer_start(T_ZSOLVE);}
	#pragma omp for
	for(j=1; j<=ny2; j++){
		double cv[PROBLEM_SIZE], rhos[PROBLEM_SIZE];
		double lhs[IMAXP+1][IMAXP+1][5];
		double lhsp[IMAXP+1][IMAXP+1][5];
		double lhsm[IMAXP+1][IMAXP+1][5];

		for(i=1; i<=nx2; i++){
			for(m=0; m<5; m++){
				lhs[0][i][m]=0.0;
				lhsp[0][i][m]=0.0;
				lhsm[0][i][m]=0.0;
				lhs[nz2+1][i][m]=0.0;
				lhsp[nz2+1][i][m]=0.0;
				lhsm[nz2+1][i][m]=0.0;
			}
			lhs[0][i][2]=1.0;
			lhsp[0][i][2]=1.0;
			lhsm[0][i][2]=1.0;
			lhs[nz2+1][i][2]=1.0;
			lhsp[nz2+1][i][2]=1.0;
			lhsm[nz2+1][i][2]=1.0;
		}

		/*
		 * ---------------------------------------------------------------------
		 * computes the left hand side for the three z-factors   
		 * ---------------------------------------------------------------------
		 * first fill the lhs for the u-eigenvalue                          
		 * ---------------------------------------------------------------------
		 */
		for(i=1; i<=nx2; i++){
			for(k=0; k<=nz2+1; k++){
				ru1=c3c4*rho_i[k][j][i];
				cv[k]=ws[k][j][i];
				rhos[k]=max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
			}
			for(k=1; k<=nz2; k++){
				lhs[k][i][0]=0.0;
				lhs[k][i][1]=-dttz2*cv[k-1]-dttz1*rhos[k-1];
				lhs[k][i][2]=1.0+c2dttz1*rhos[k];
				lhs[k][i][3]=dttz2*cv[k+1]-dttz1*rhos[k+1];
				lhs[k][i][4]=0.0;
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * add fourth order dissipation                                  
		 * ---------------------------------------------------------------------
		 */
		for(i=1; i<=nx2; i++){
			k=1;
			lhs[k][i][2]=lhs[k][i][2]+comz5;
			lhs[k][i][3]=lhs[k][i][3]-comz4;
			lhs[k][i][4]=lhs[k][i][4]+comz1;
			k=2;
			lhs[k][i][1]=lhs[k][i][1]-comz4;
			lhs[k][i][2]=lhs[k][i][2]+comz6;
			lhs[k][i][3]=lhs[k][i][3]-comz4;
			lhs[k][i][4]=lhs[k][i][4]+comz1;
		}
		for(k=3; k<=nz2-2; k++){
			for(i=1; i<=nx2; i++){
				lhs[k][i][0]=lhs[k][i][0]+comz1;
				lhs[k][i][1]=lhs[k][i][1]-comz4;
				lhs[k][i][2]=lhs[k][i][2]+comz6;
				lhs[k][i][3]=lhs[k][i][3]-comz4;
				lhs[k][i][4]=lhs[k][i][4]+comz1;
			}
		}
		for(i=1; i<=nx2; i++){
			k=nz2-1;
			lhs[k][i][0]=lhs[k][i][0]+comz1;
			lhs[k][i][1]=lhs[k][i][1]-comz4;
			lhs[k][i][2]=lhs[k][i][2]+comz6;
			lhs[k][i][3]=lhs[k][i][3]-comz4;
			k=nz2;
			lhs[k][i][0]=lhs[k][i][0]+comz1;
			lhs[k][i][1]=lhs[k][i][1]-comz4;
			lhs[k][i][2]=lhs[k][i][2]+comz5;
		}
		/*
		 * ---------------------------------------------------------------------
		 * subsequently, fill the other factors (u+c), (u-c) 
		 * ---------------------------------------------------------------------
		 */
		for(k=1; k<=nz2; k++){
			for(i=1; i<=nx2; i++){
				lhsp[k][i][0]=lhs[k][i][0];
				lhsp[k][i][1]=lhs[k][i][1]-dttz2*speed[k-1][j][i];
				lhsp[k][i][2]=lhs[k][i][2];
				lhsp[k][i][3]=lhs[k][i][3]+dttz2*speed[k+1][j][i];
				lhsp[k][i][4]=lhs[k][i][4];
				lhsm[k][i][0]=lhs[k][i][0];
				lhsm[k][i][1]=lhs[k][i][1]+dttz2*speed[k-1][j][i];
				lhsm[k][i][2]=lhs[k][i][2];
				lhsm[k][i][3]=lhs[k][i][3]-dttz2*speed[k+1][j][i];
				lhsm[k][i][4]=lhs[k][i][4];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * FORWARD ELIMINATION  
		 * ---------------------------------------------------------------------
		 */
		for(k=0; k<=grid_points[2]-3; k++){
			k1=k+1;
			k2=k+2;
			for(i=1; i<=nx2; i++){
				fac1=1.0/lhs[k][i][2];
				lhs[k][i][3]=fac1*lhs[k][i][3];
				lhs[k][i][4]=fac1*lhs[k][i][4];
				for(m=0; m<3; m++){
					rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				}
				lhs[k1][i][2]=lhs[k1][i][2]-lhs[k1][i][1]*lhs[k][i][3];
				lhs[k1][i][3]=lhs[k1][i][3]-lhs[k1][i][1]*lhs[k][i][4];
				for(m=0; m<3; m++){
					rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhs[k1][i][1]*rhs[k][j][i][m];
				}
				lhs[k2][i][1]=lhs[k2][i][1]-lhs[k2][i][0]*lhs[k][i][3];
				lhs[k2][i][2]=lhs[k2][i][2]-lhs[k2][i][0]*lhs[k][i][4];
				for(m=0; m<3; m++){
					rhs[k2][j][i][m]=rhs[k2][j][i][m]-lhs[k2][i][0]*rhs[k][j][i][m];
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
		for(i=1; i<=nx2; i++){
			fac1=1.0/lhs[k][i][2];
			lhs[k][i][3]=fac1*lhs[k][i][3];
			lhs[k][i][4]=fac1*lhs[k][i][4];
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			}
			lhs[k1][i][2]=lhs[k1][i][2]-lhs[k1][i][1]*lhs[k][i][3];
			lhs[k1][i][3]=lhs[k1][i][3]-lhs[k1][i][1]*lhs[k][i][4];
			for(m=0; m<3; m++){
				rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhs[k1][i][1]*rhs[k][j][i][m];
			}
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately
			 * ---------------------------------------------------------------------
			 */
			fac2=1.0/lhs[k1][i][2];
			for(m=0; m<3; m++){
				rhs[k1][j][i][m]=fac2*rhs[k1][j][i][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * do the u+c and the u-c factors               
		 * ---------------------------------------------------------------------
		 */
		for(k=0; k<=grid_points[2]-3; k++){
			k1=k+1;
			k2=k+2;
			for(i=1; i<=nx2; i++){
				m=3;
				fac1=1.0/lhsp[k][i][2];
				lhsp[k][i][3]=fac1*lhsp[k][i][3];
				lhsp[k][i][4]=fac1*lhsp[k][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsp[k1][i][2]=lhsp[k1][i][2]-lhsp[k1][i][1]*lhsp[k][i][3];
				lhsp[k1][i][3]=lhsp[k1][i][3]-lhsp[k1][i][1]*lhsp[k][i][4];
				rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhsp[k1][i][1]*rhs[k][j][i][m];
				lhsp[k2][i][1]=lhsp[k2][i][1]-lhsp[k2][i][0]*lhsp[k][i][3];
				lhsp[k2][i][2]=lhsp[k2][i][2]-lhsp[k2][i][0]*lhsp[k][i][4];
				rhs[k2][j][i][m]=rhs[k2][j][i][m]-lhsp[k2][i][0]*rhs[k][j][i][m];
				m=4;
				fac1=1.0/lhsm[k][i][2];
				lhsm[k][i][3]=fac1*lhsm[k][i][3];
				lhsm[k][i][4]=fac1*lhsm[k][i][4];
				rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
				lhsm[k1][i][2]=lhsm[k1][i][2]-lhsm[k1][i][1]*lhsm[k][i][3];
				lhsm[k1][i][3]=lhsm[k1][i][3]-lhsm[k1][i][1]*lhsm[k][i][4];
				rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhsm[k1][i][1]*rhs[k][j][i][m];
				lhsm[k2][i][1]=lhsm[k2][i][1]-lhsm[k2][i][0]*lhsm[k][i][3];
				lhsm[k2][i][2]=lhsm[k2][i][2]-lhsm[k2][i][0]*lhsm[k][i][4];
				rhs[k2][j][i][m]=rhs[k2][j][i][m]-lhsm[k2][i][0]*rhs[k][j][i][m];
			}
		}
		/*
		 * ---------------------------------------------------------------------
		 * and again the last two rows separately
		 * ---------------------------------------------------------------------
		 */
		k=grid_points[2]-2;
		k1=grid_points[2]-1;
		for(i=1; i<=nx2; i++){
			m=3;
			fac1=1.0/lhsp[k][i][2];
			lhsp[k][i][3]=fac1*lhsp[k][i][3];
			lhsp[k][i][4]=fac1*lhsp[k][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsp[k1][i][2]=lhsp[k1][i][2]-lhsp[k1][i][1]*lhsp[k][i][3];
			lhsp[k1][i][3]=lhsp[k1][i][3]-lhsp[k1][i][1]*lhsp[k][i][4];
			rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhsp[k1][i][1]*rhs[k][j][i][m];
			m=4;
			fac1=1.0/lhsm[k][i][2];
			lhsm[k][i][3]=fac1*lhsm[k][i][3];
			lhsm[k][i][4]=fac1*lhsm[k][i][4];
			rhs[k][j][i][m]=fac1*rhs[k][j][i][m];
			lhsm[k1][i][2]=lhsm[k1][i][2]-lhsm[k1][i][1]*lhsm[k][i][3];
			lhsm[k1][i][3]=lhsm[k1][i][3]-lhsm[k1][i][1]*lhsm[k][i][4];
			rhs[k1][j][i][m]=rhs[k1][j][i][m]-lhsm[k1][i][1]*rhs[k][j][i][m];
			/*
			 * ---------------------------------------------------------------------
			 * scale the last row immediately (some of this is overkill
			 * if this is the last cell)
			 * ---------------------------------------------------------------------
			 */
			rhs[k1][j][i][3]=rhs[k1][j][i][3]/lhsp[k1][i][2];
			rhs[k1][j][i][4]=rhs[k1][j][i][4]/lhsm[k1][i][2];
		}
		/*
		 * ---------------------------------------------------------------------
		 * BACKSUBSTITUTION 
		 * ---------------------------------------------------------------------
		 */
		k=grid_points[2]-2;
		k1=grid_points[2]-1;
		for(i=1; i<=nx2; i++){
			for(m=0; m<3; m++){
				rhs[k][j][i][m]=rhs[k][j][i][m]-lhs[k][i][3]*rhs[k1][j][i][m];
			}
			rhs[k][j][i][3]=rhs[k][j][i][3]-lhsp[k][i][3]*rhs[k1][j][i][3];
			rhs[k][j][i][4]=rhs[k][j][i][4]-lhsm[k][i][3]*rhs[k1][j][i][4];
		}
		/*
		 * ---------------------------------------------------------------------
		 * whether or not this is the last processor, we always have
		 * to complete the back-substitution 
		 * ---------------------------------------------------------------------
		 * the first three factors
		 * ---------------------------------------------------------------------
		 */
		for(k=grid_points[2]-3; k>=0; k--){
			k1=k+1;
			k2=k+2;
			for(i=1; i<=nx2; i++){
				for (m = 0; m < 3; m++) {
					rhs[k][j][i][m]=rhs[k][j][i][m]- 
						lhs[k][i][3]*rhs[k1][j][i][m]-
						lhs[k][i][4]*rhs[k2][j][i][m];
				}
				/*
				 * ---------------------------------------------------------------------
				 * and the remaining two
				 * ---------------------------------------------------------------------
				 */
				rhs[k][j][i][3]=rhs[k][j][i][3]- 
					lhsp[k][i][3]*rhs[k1][j][i][3]-
					lhsp[k][i][4]*rhs[k2][j][i][3];
				rhs[k][j][i][4]=rhs[k][j][i][4]- 
					lhsm[k][i][3]*rhs[k1][j][i][4]-
					lhsm[k][i][4]*rhs[k2][j][i][4];
			}
		}
	}
	if(timeron && thread_id==0){timer_stop(T_ZSOLVE);}
	tzetar();
}
