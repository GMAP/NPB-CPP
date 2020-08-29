/**
 * NASA Advanced Supercomputing Parallel Benchmarks C++
 *
 * based on NPB 3.3.1
 *
 * original version and technical report:
 * http://www.nas.nasa.gov/Software/NPB/
 *
 * Authors:
 *     P. O. Frederickson
 *     D. H. Bailey
 *     A. C. Woo
 *
 * C++ version:
 *      Dalvan Griebler <dalvangriebler@gmail.com>
 *      Gabriell Alves de Araujo <hexenoften@gmail.com>
 *      Júnior Löff <loffjh@gmail.com>
 */

#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"

/*
 * --------------------------------------------------------------------
 * this is the serial version of the app benchmark 1,
 * the "embarassingly parallel" benchmark.
 * --------------------------------------------------------------------
 * M is the Log_2 of the number of complex pairs of uniform (0, 1) random
 * numbers. MK is the Log_2 of the size of each batch of uniform random
 * numbers.  MK can be set for convenience on a given system, since it does
 * not affect the results.
 * --------------------------------------------------------------------
 */
#define	MK 16
#define	MM (M - MK)
#define	NN (1 << MM)
#define	NK (1 << MK)
#define	NQ 10
#define EPSILON 1.0e-8
#define	A 1220703125.0
#define	S 271828183.0
#define NK_PLUS ((2*NK)+1)

/* global variables */
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
static double x[NK_PLUS];
static double q[NQ];
#else
static double (*x)=(double*)malloc(sizeof(double)*(NK_PLUS));
static double (*q)=(double*)malloc(sizeof(double)*(NQ));
#endif

/* ep */
int main(int argc, char **argv){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	double  Mops, t1, t2, t3, t4, x1, x2;
	double  sx, sy, tm, an, tt, gc;
	double  sx_verify_value, sy_verify_value, sx_err, sy_err;
	int     np;
	int     i, ik, kk, l, k, nit;
	int     k_offset, j;
	boolean verified, timers_enabled;
	double  dum[3] = {1.0, 1.0, 1.0};
	char    size[16];

	FILE* fp;
	if((fp = fopen("timer.flag", "r"))==NULL){
		timers_enabled = FALSE;
	}else{
		timers_enabled = TRUE;
		fclose(fp);
	}

	/*
	 * --------------------------------------------------------------------
	 * because the size of the problem is too large to store in a 32-bit
	 * integer for some classes, we put it into a string (for printing).
	 * have to strip off the decimal point put in there by the floating
	 * point print statement (internal file)
	 * --------------------------------------------------------------------
	 */    
	sprintf(size, "%15.0f", pow(2.0, M+1));
	j = 14;
	if(size[j]=='.'){j--;}
	size[j+1] = '\0';
	printf("\n\n NAS Parallel Benchmarks 4.1 Serial C++ version - EP Benchmark\n\n");
	printf(" Number of random numbers generated: %15s\n", size);

	verified = FALSE;

	/*
	 * --------------------------------------------------------------------
	 * compute the number of "batches" of random number pairs generated 
	 * per processor. Adjust if the number of processors does not evenly 
	 * divide the total number
	 * --------------------------------------------------------------------
	 */
	np = NN;

	/*
	 * call the random number generator functions and initialize
	 * the x-array to reduce the effects of paging on the timings.
	 * also, call all mathematical functions that are used. make
	 * sure these initializations cannot be eliminated as dead code.
	 */
	vranlc(0, &dum[0], dum[1], &dum[2]);
	dum[0] = randlc(&dum[1], dum[2]);
	for(i=0; i<NK_PLUS; i++){x[i] = -1.0e99;}
	Mops = log(sqrt(fabs(max(1.0, 1.0))));

	timer_clear(0);
	timer_clear(1);
	timer_clear(2);
	timer_start(0);

	t1 = A;
	vranlc(0, &t1, A, x);

	/* compute AN = A ^ (2 * NK) (mod 2^46) */

	t1 = A;

	for(i=0; i<MK+1; i++){
		t2 = randlc(&t1, t1);
	}

	an = t1;
	tt = S;
	gc = 0.0;
	sx = 0.0;
	sy = 0.0;

	for(i=0; i<=NQ-1; i++){
		q[i] = 0.0;
	}

	/*
	 * each instance of this loop may be performed independently. we compute
	 * the k offsets separately to take into account the fact that some nodes
	 * have more numbers to generate than others
	 */
	k_offset = -1;

	for(k=1; k<=np; k++){
		kk = k_offset + k;
		t1 = S;
		t2 = an;

		/* find starting seed t1 for this kk */
		for(i=1; i<=100; i++){
			ik = kk / 2;
			if((2*ik)!=kk){t3=randlc(&t1,t2);}
			if(ik==0){break;}
			t3=randlc(&t2,t2);
			kk=ik;
		}

		/* compute uniform pseudorandom numbers */
		if(timers_enabled){timer_start(2);}
		vranlc(2*NK, &t1, A, x);
		if(timers_enabled){timer_stop(2);}

		/*
		 * compute gaussian deviates by acceptance-rejection method and
		 * tally counts in concentric square annuli. this loop is not
		 * vectorizable.
		 */

		if(timers_enabled){timer_start(1);}

		for(i=0; i<NK; i++){
			x1 = 2.0 * x[2*i] - 1.0;
			x2 = 2.0 * x[2*i+1] - 1.0;
			t1 = pow2(x1) + pow2(x2);
			if(t1 <= 1.0){
				t2 = sqrt(-2.0 * log(t1) / t1);
				t3 = (x1 * t2);
				t4 = (x2 * t2);
				l = max(fabs(t3), fabs(t4));
				q[l] += 1.0;
				sx = sx + t3;
				sy = sy + t4;
			}
		}
		if(timers_enabled){timer_stop(1);}
	}

	for(i=0; i<=NQ-1; i++){
		gc = gc + q[i];
	}

	timer_stop(0);
	tm = timer_read(0);

	nit = 0;
	verified = TRUE;
	if(M == 24){
		sx_verify_value = -3.247834652034740e+3;
		sy_verify_value = -6.958407078382297e+3;
	}else if(M == 25){
		sx_verify_value = -2.863319731645753e+3;
		sy_verify_value = -6.320053679109499e+3;
	}else if(M == 28){
		sx_verify_value = -4.295875165629892e+3;
		sy_verify_value = -1.580732573678431e+4;
	}else if(M == 30){
		sx_verify_value =  4.033815542441498e+4;
		sy_verify_value = -2.660669192809235e+4;
	}else if(M == 32){
		sx_verify_value =  4.764367927995374e+4;
		sy_verify_value = -8.084072988043731e+4;
	}else if(M == 36){
		sx_verify_value =  1.982481200946593e+5;
		sy_verify_value = -1.020596636361769e+5;
	}else if (M == 40){
		sx_verify_value = -5.319717441530e+05;
		sy_verify_value = -3.688834557731e+05;
	}else{
		verified = FALSE;
	}
	if(verified){
		sx_err = fabs((sx - sx_verify_value) / sx_verify_value);
		sy_err = fabs((sy - sy_verify_value) / sy_verify_value);
		verified = ((sx_err <= EPSILON) && (sy_err <= EPSILON));
	}
	Mops = pow(2.0, M+1)/tm/1000000.0;

	printf("\n EP Benchmark Results:\n\n");
	printf(" CPU Time =%10.4f\n", tm);
	printf(" N = 2^%5d\n", M);
	printf(" No. Gaussian Pairs = %15.0f\n", gc);
	printf(" Sums = %25.15e %25.15e\n", sx, sy);
	printf(" Counts: \n");
	for(i=0; i<NQ-1; i++){
		printf("%3d%15.0f\n", i, q[i]);
	}

	c_print_results((char*)"EP",
			CLASS,
			M+1,
			0,
			0,
			nit,
			tm,
			Mops,
			(char*)"Random numbers generated",
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
			(char*)CS7);

	if(timers_enabled){
		if(tm <= 0.0){tm = 1.0;}
		tt = timer_read(0);
		printf("\nTotal time:     %9.3f (%6.2f)\n", tt, tt*100.0/tm);
		tt = timer_read(1);
		printf("Gaussian pairs: %9.3f (%6.2f)\n", tt, tt*100.0/tm);
		tt = timer_read(2);
		printf("Random numbers: %9.3f (%6.2f)\n", tt, tt*100.0/tm);
	}

	return 0;
}
