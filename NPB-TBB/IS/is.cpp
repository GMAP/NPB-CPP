/*
NASA Open Source Agreement (NOSA)

Copyright (c) 2019 NASA Advanced Supercomputing (NAS) Division 
	NASA website: http://www.nas.nasa.gov/Software/NPB/
	
	NAS Parallel Benchmarks Group
	NASA Ames Research Center
	Mail Stop: T27A-1
	Moffett Field, CA   94035-1000
	E-mail:  npb@nas.nasa.gov
	Fax:     (650) 604-3957

License agreements: https://opensource.org/licenses/nasa1.3.php

--------------------------------------------------------------------------

The original NPB 3.4.1 version belongs to: 
	http://www.nas.nasa.gov/Software/NPB/

Authors of the C code:
	M. Yarrow
	H. Jin

------------------------------------------------------------------------------

The Intel TBB version is a parallel implementation of the serial C version
Intel TBB version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-TBB

Authors of the Intel TBB code:
	Júnior Löff <loffjh@gmail.com>
	
*/

#include "thread"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/mutex.h"
#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"

#define T_BENCHMARKING (0)
#define T_INITIALIZATION (1)
#define T_SORTING (2)
#define T_TOTAL_EXECUTION (3)

/*****************************************************************/
/* For serial IS, buckets are not really req'd to solve NPB1 IS  */
/* spec, but their use on some machines improves performance, on */
/* other machines the use of buckets compromises performance,    */
/* probably because it is extra computation which is not req'd.  */
/* (Note: Mechanism not understood, probably cache related)      */
/* Example:  SP2-66MhzWN:  50% speedup with buckets              */
/* Example:  SGI Indy5000: 50% slowdown with buckets             */
/* Example:  SGI O2000:   400% slowdown with buckets (Wow!)      */
/*****************************************************************/
/* To disable the use of buckets, comment out the following line */
#define USE_BUCKETS

/******************/
/* default values */
/******************/
#ifndef CLASS
#define CLASS 'S'
#endif

/*************/
/*  CLASS S  */
/*************/
#if CLASS == 'S'
#define TOTAL_KEYS_LOG_2 16
#define MAX_KEY_LOG_2 11
#define NUM_BUCKETS_LOG_2 9
#endif

/*************/
/*  CLASS W  */
/*************/
#if CLASS == 'W'
#define TOTAL_KEYS_LOG_2 20
#define MAX_KEY_LOG_2 16
#define NUM_BUCKETS_LOG_2 10
#endif

/*************/
/*  CLASS A  */
/*************/
#if CLASS == 'A'
#define TOTAL_KEYS_LOG_2 23
#define MAX_KEY_LOG_2 19
#define NUM_BUCKETS_LOG_2 10
#endif

/*************/
/*  CLASS B  */
/*************/
#if CLASS == 'B'
#define TOTAL_KEYS_LOG_2 25
#define MAX_KEY_LOG_2 21
#define NUM_BUCKETS_LOG_2 10
#endif

/*************/
/*  CLASS C  */
/*************/
#if CLASS == 'C'
#define TOTAL_KEYS_LOG_2 27
#define MAX_KEY_LOG_2 23
#define NUM_BUCKETS_LOG_2 10
#endif

/*************/
/*  CLASS D  */
/*************/
#if CLASS == 'D'
#define TOTAL_KEYS_LOG_2 31
#define MAX_KEY_LOG_2 27
#define NUM_BUCKETS_LOG_2 10
#endif

#if CLASS == 'D'
#define TOTAL_KEYS (1L << TOTAL_KEYS_LOG_2)
#else
#define TOTAL_KEYS (1 << TOTAL_KEYS_LOG_2)
#endif
#define MAX_KEY (1 << MAX_KEY_LOG_2)
#define NUM_BUCKETS (1 << NUM_BUCKETS_LOG_2)
#define NUM_KEYS TOTAL_KEYS
#define SIZE_OF_BUFFERS NUM_KEYS                                           

#define MAX_ITERATIONS 10
#define TEST_ARRAY_SIZE 5

/*************************************/
/* Typedef: if necessary, change the */
/* size of int here by changing the  */
/* int type to, say, long            */
/*************************************/
#if CLASS == 'D'
typedef long INT_TYPE;
#else
typedef int INT_TYPE;
#endif

/********************/
/* Some global info */
/********************/
INT_TYPE* key_buff_ptr_global; /* used by full_verify to get copies of rank info */
int passed_verification;                                 

/************************************/
/* These are the three main arrays. */
/* See SIZE_OF_BUFFERS def above    */
/************************************/
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
INT_TYPE key_array[SIZE_OF_BUFFERS];
INT_TYPE key_buff1[MAX_KEY];
INT_TYPE key_buff2[SIZE_OF_BUFFERS];
INT_TYPE partial_verify_vals[TEST_ARRAY_SIZE];
INT_TYPE** key_buff1_aptr = NULL;
#else
INT_TYPE (*key_array)=(INT_TYPE*)malloc(sizeof(INT_TYPE)*(SIZE_OF_BUFFERS));    
INT_TYPE (*key_buff1)=(INT_TYPE*)malloc(sizeof(INT_TYPE)*(MAX_KEY));                
INT_TYPE (*key_buff2)=(INT_TYPE*)malloc(sizeof(INT_TYPE)*(SIZE_OF_BUFFERS));
INT_TYPE (*partial_verify_vals)=(INT_TYPE*)malloc(sizeof(INT_TYPE)*(TEST_ARRAY_SIZE));       
INT_TYPE** key_buff1_aptr = NULL;
#endif

#ifdef USE_BUCKETS
INT_TYPE** bucket_size; 
INT_TYPE bucket_ptrs[NUM_BUCKETS];
#endif

int num_workers;
tbb::mutex critical_section;

pthread_barrier_t worker_barrier; 

/**********************/
/* Partial verif info */
/**********************/
INT_TYPE test_index_array[TEST_ARRAY_SIZE],
	 test_rank_array[TEST_ARRAY_SIZE],

	 S_test_index_array[TEST_ARRAY_SIZE] = 
{48427,17148,23627,62548,4431},
	S_test_rank_array[TEST_ARRAY_SIZE] = 
{0,18,346,64917,65463},

	W_test_index_array[TEST_ARRAY_SIZE] = 
{357773,934767,875723,898999,404505},
	W_test_rank_array[TEST_ARRAY_SIZE] = 
{1249,11698,1039987,1043896,1048018},

	A_test_index_array[TEST_ARRAY_SIZE] = 
{2112377,662041,5336171,3642833,4250760},
	A_test_rank_array[TEST_ARRAY_SIZE] = 
{104,17523,123928,8288932,8388264},

	B_test_index_array[TEST_ARRAY_SIZE] = 
{41869,812306,5102857,18232239,26860214},
	B_test_rank_array[TEST_ARRAY_SIZE] = 
{33422937,10244,59149,33135281,99}, 

	C_test_index_array[TEST_ARRAY_SIZE] = 
{44172927,72999161,74326391,129606274,21736814},
	C_test_rank_array[TEST_ARRAY_SIZE] = 
{61147,882988,266290,133997595,133525895},

	D_test_index_array[TEST_ARRAY_SIZE] = 
{1317351170,995930646,1157283250,1503301535,1453734525},
	D_test_rank_array[TEST_ARRAY_SIZE] = 
{1,36538729,1978098519,2145192618,2147425337};

/***********************/
/* function prototypes */
/***********************/
void alloc_key_buff();
void* alloc_mem(size_t size);
void create_seq(double seed,
		double a);
double find_my_seed(int kn,
		int np,
		long nn,
		double s,
		double a );
void full_verify();    
void rank(int iteration);

/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/
int main(int argc, char** argv){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	int i, iteration, timer_on;
	double timecounter;
	FILE* fp;

	if(const char * nw = std::getenv("TBB_NUM_THREADS")) {
        num_workers = atoi(nw);
    } else {
        num_workers = 1;
    }
    
    tbb::task_scheduler_init init(num_workers);
	pthread_barrier_init(&worker_barrier,NULL,num_workers);
	
	/* Initialize timers */
	timer_on = 0;            
	if((fp = fopen("timer.flag", "r")) != NULL){
		fclose(fp);
		timer_on = 1;
	}
	timer_clear( T_BENCHMARKING );
	if(timer_on){
		timer_clear( T_INITIALIZATION );
		timer_clear( T_SORTING );
		timer_clear( T_TOTAL_EXECUTION );
	}

	if(timer_on)timer_start( T_TOTAL_EXECUTION );

	/* Initialize the verification arrays if a valid class */
	for( i=0; i<TEST_ARRAY_SIZE; i++ )
		switch( CLASS )
		{
			case 'S':
				test_index_array[i] = S_test_index_array[i];
				test_rank_array[i]  = S_test_rank_array[i];
				break;
			case 'A':
				test_index_array[i] = A_test_index_array[i];
				test_rank_array[i]  = A_test_rank_array[i];
				break;
			case 'W':
				test_index_array[i] = W_test_index_array[i];
				test_rank_array[i]  = W_test_rank_array[i];
				break;
			case 'B':
				test_index_array[i] = B_test_index_array[i];
				test_rank_array[i]  = B_test_rank_array[i];
				break;
			case 'C':
				test_index_array[i] = C_test_index_array[i];
				test_rank_array[i]  = C_test_rank_array[i];
				break;
			case 'D':
				test_index_array[i] = D_test_index_array[i];
				test_rank_array[i]  = D_test_rank_array[i];
				break;
		};        

	/* Printout initial NPB info */
	printf("\n\n NAS Parallel Benchmarks 4.1 Parallel C++ version with Intel TBB - IS Benchmark\n\n");
	printf(" Size:  %ld  (class %c)\n", (long)TOTAL_KEYS, CLASS);
	printf(" Iterations:   %d\n", MAX_ITERATIONS);
	printf( "\n" );

	if(timer_on)timer_start( T_INITIALIZATION );

	/* Generate random number sequence and subsequent keys on all procs */
	create_seq(314159265.00 /* Random number gen seed */, 
			1220703125.00 /* Random number gen mult */);                 

	alloc_key_buff();
	if(timer_on)timer_stop( T_INITIALIZATION );

	/* Do one interation for free (i.e., untimed) to guarantee initialization of */
	/* all data and code pages and respective tables */
	rank( 1 );  

	/* Start verification counter */
	passed_verification = 0;

	if( CLASS != 'S' ) printf( "\n   iteration\n" );

	/* Start timer */             
	timer_start( T_BENCHMARKING );

	/* This is the main iteration */
	for(iteration=1; iteration<=MAX_ITERATIONS; iteration++){
		if(CLASS != 'S')printf("        %d\n", iteration);
		rank( iteration );
	}

	/* End of timing, obtain maximum time of all processors */
	timer_stop( T_BENCHMARKING );
	timecounter = timer_read( T_BENCHMARKING );

	/* This tests that keys are in sequence: sorting of last ranked key seq */
	/* occurs here, but is an untimed operation */
	if(timer_on)timer_start( T_SORTING );
	full_verify();
	if(timer_on)timer_stop( T_SORTING );

	if(timer_on)timer_stop( T_TOTAL_EXECUTION );    

	/* The final printout */
	if(passed_verification != 5*MAX_ITERATIONS + 1){passed_verification = 0;}
	setenv("TBB_NUM_THREADS","1",0);
	c_print_results((char*)"IS",
			CLASS,
			(int)(TOTAL_KEYS/64),
			64,
			0,
			MAX_ITERATIONS,
			timecounter,
			((double)(MAX_ITERATIONS*TOTAL_KEYS))/timecounter/1000000.0,
			(char*)"keys ranked",
			passed_verification,
			(char*)NPBVERSION,
			(char*)COMPILETIME,
			(char*)COMPILERVERSION,
			(char*)LIBVERSION,
			std::getenv("TBB_NUM_THREADS"),
			(char*)CS1,
			(char*)CS2,
			(char*)CS3,
			(char*)CS4,
			(char*)CS5,
			(char*)CS6,
			(char*)CS7);

	/* Print additional timers */
	if(timer_on){
		double t_total, t_percent;
		t_total = timer_read( T_TOTAL_EXECUTION );
		printf("\nAdditional timers -\n");
		printf(" Total execution: %8.3f\n", t_total);
		if (t_total == 0.0) t_total = 1.0;
		timecounter = timer_read(T_INITIALIZATION);
		t_percent = timecounter/t_total * 100.;
		printf(" Initialization : %8.3f (%5.2f%%)\n", timecounter, t_percent);
		timecounter = timer_read(T_BENCHMARKING);
		t_percent = timecounter/t_total * 100.;
		printf(" Benchmarking   : %8.3f (%5.2f%%)\n", timecounter, t_percent);
		timecounter = timer_read(T_SORTING);
		t_percent = timecounter/t_total * 100.;
		printf(" Sorting        : %8.3f (%5.2f%%)\n", timecounter, t_percent);
	}

	return 0;
}

void alloc_key_buff(){
	INT_TYPE i;
	int num_procs;

	num_procs = num_workers;

#ifdef USE_BUCKETS
	bucket_size = (INT_TYPE**)alloc_mem(sizeof(INT_TYPE*)*num_procs);

	for(i = 0; i < num_procs; i++){
		bucket_size[i] = (INT_TYPE*)alloc_mem(sizeof(INT_TYPE)*NUM_BUCKETS);
	}

	for( i=0; i<NUM_KEYS; i++ )
		key_buff2[i] = 0;
#else /*USE_BUCKETS*/
	key_buff1_aptr = (INT_TYPE**)alloc_mem(sizeof(INT_TYPE*)*num_procs);

	key_buff1_aptr[0] = key_buff1;
	for(i = 1; i < num_procs; i++) {
		key_buff1_aptr[i] = (INT_TYPE *)alloc_mem(sizeof(INT_TYPE) * MAX_KEY);
	}
#endif /*USE_BUCKETS*/
}

/*****************************************************************/
/*****************    Allocate Working Buffer     ****************/
/*****************************************************************/
void* alloc_mem(size_t size){
	void* p;
	p = (void*)malloc(size);
	if(!p){
		perror("Memory allocation error");
		exit(1);
	}
	return p;
}

/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/
void create_seq(double seed,
		double a){

	tbb::parallel_for(tbb::blocked_range<size_t>(0,num_workers),[&](const tbb::blocked_range<size_t>& r){
		int worker_id;
		for(worker_id=r.begin(); worker_id != r.end(); worker_id++){
			double x, s;
			INT_TYPE i, k;

			INT_TYPE k1, k2;
			double an = a;
			int myid, num_procs;
			INT_TYPE mq;

			myid = worker_id;
			num_procs = num_workers;

			mq = (NUM_KEYS + num_procs - 1) / num_procs;
			k1 = mq * myid;
			k2 = k1 + mq;
			if ( k2 > NUM_KEYS ) k2 = NUM_KEYS;

			s = find_my_seed( myid, 
					num_procs,
					(long)4*NUM_KEYS,
					seed,
					an );

			k = MAX_KEY/4;

			for(i=k1; i<k2; i++){
				x = randlc(&s, an);
				x += randlc(&s, an);
				x += randlc(&s, an);
				x += randlc(&s, an);
				key_array[i] = k*x;
			}
		}
	});
}

/*****************************************************************/
/************   F  I  N  D  _  M  Y  _  S  E  E  D    ************/
/************                                         ************/
/************ returns parallel random number seq seed ************/
/*****************************************************************/
double find_my_seed(int kn, /* my processor rank, 0<=kn<=num procs */
		int np, /* np = num procs */
		long nn, /* total num of ran numbers, all procs */
		double s, /* Ran num seed, for ex.: 314159265.00 */
		double a){ /* Ran num gen mult, try 1220703125.00 */
	/*
	 * Create a random number sequence of total length nn residing
	 * on np number of processors.  Each processor will therefore have a
	 * subsequence of length nn/np.  This routine returns that random
	 * number which is the first random number for the subsequence belonging
	 * to processor rank kn, and which is used as seed for proc kn ran # gen.
	 */
	double t1,t2;
	long mq,nq,kk,ik;

	if ( kn == 0 ) return s;

	mq = (nn/4 + np - 1) / np;
	nq = mq * 4 * kn; /* number of rans to be skipped */

	t1 = s;
	t2 = a;
	kk = nq;
	while( kk > 1 ){
		ik = kk / 2;
		if(2 * ik ==  kk){
			(void)randlc( &t2, t2 );
			kk = ik;
		}
		else{
			(void)randlc( &t1, t2 );
			kk = kk - 1;
		}
	}
	(void)randlc( &t1, t2 );

	return( t1 );
}

/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/
void full_verify(){
	INT_TYPE i, j;
	INT_TYPE k, k1, k2;
	int myid, num_procs;

	myid = 0;
	num_procs = 1;

	/* Now, finally, sort the keys: */
	/* Copy keys into work array; keys in key_array will be reassigned. */

#ifdef USE_BUCKETS
	/* Buckets are already sorted. Sorting keys within each bucket */
	tbb::parallel_for(tbb::blocked_range<size_t>(0,NUM_BUCKETS),[&](const tbb::blocked_range<size_t>& r){
		INT_TYPE j,i,k,k1;
		for(j=r.begin(); j != r.end(); j++){
			k1 = (j > 0)? bucket_ptrs[j-1] : 0;
			for ( i = k1; i < bucket_ptrs[j]; i++ ) {
				k = --key_buff_ptr_global[key_buff2[i]];
				key_array[k] = key_buff2[i];
			}
		}
	});
#else    
	for( i=0; i<NUM_KEYS; i++ )
		key_buff2[i] = key_array[i];
	/* This is actual sorting. Each thread is responsible for a subset of key values */
	j = num_procs;
	j = (MAX_KEY + j - 1) / j;
	k1 = j * myid;
	k2 = k1 + j;
	if (k2 > MAX_KEY) k2 = MAX_KEY;
	for( i=0; i<NUM_KEYS; i++ ) {
		if (key_buff2[i] >= k1 && key_buff2[i] < k2) {
			k = --key_buff_ptr_global[key_buff2[i]];
			key_array[k] = key_buff2[i];
		}
	}
#endif

	/* Confirm keys correctly sorted: count incorrectly sorted keys, if any */
	j = 0;

	tbb::parallel_for(tbb::blocked_range<size_t>(1,NUM_KEYS),[&](const tbb::blocked_range<size_t>& r){
		INT_TYPE j_worker,i;
		j_worker = 0;

		for(i=r.begin(); i != r.end(); i++){
			if( key_array[i-1] > key_array[i] )
				j_worker++;
		}

		critical_section.lock();
			j += j_worker;
        critical_section.unlock();
	});

	if( j != 0 )
		printf( "Full_verify: number of keys out of sort: %ld\n", (long)j );
	else
		passed_verification++;
}

/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/
void rank(int iteration){

	INT_TYPE i, k;
	INT_TYPE *key_buff_ptr, *key_buff_ptr2;

#ifdef USE_BUCKETS
	int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
	INT_TYPE num_bucket_keys = (1L << shift);
#endif

	key_array[iteration] = iteration;
	key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;

	/* Determine where the partial verify test keys are, load into */
	/* top of array bucket_size */
	for( i=0; i<TEST_ARRAY_SIZE; i++ )
		partial_verify_vals[i] = key_array[test_index_array[i]];

	/* Setup pointers to key buffers */
#ifdef USE_BUCKETS
	key_buff_ptr2 = key_buff2;
#else
	key_buff_ptr2 = key_array;
#endif
	key_buff_ptr = key_buff1;

	INT_TYPE *work_buff, m, k1, k2;
	int myid = 0, num_procs = 1;

	/* Bucket sort is known to improve cache performance on some */
	/* cache based systems.  But the actual performance may depend */
	/* on cache size, problem size. */
#ifdef USE_BUCKETS
	tbb::parallel_for(tbb::blocked_range<size_t>(0,num_workers),[&](const tbb::blocked_range<size_t>& r){
		int worker_id;
		INT_TYPE i, k, j;
		INT_TYPE *work_buff, m, k1, k2;
		INT_TYPE bucket_ptrs_worker[NUM_BUCKETS];
	
		for(worker_id=r.begin(); worker_id != r.end(); worker_id++){
			work_buff = bucket_size[worker_id];

			/* Initialize */
			for( i=0; i<NUM_BUCKETS; i++ )  
				work_buff[i] = 0;

			int chunk = NUM_KEYS/num_workers;
			int residual = NUM_KEYS%num_workers;

			int ist_worker = worker_id*chunk + ((residual>1 && (worker_id<residual)) ? worker_id : residual);
			int iend_worker = worker_id*chunk + chunk + ((residual>1 && (worker_id<residual)) ? worker_id+1 : residual);


			if(worker_id==0) ist_worker = 0;
			if(worker_id==num_workers-1) iend_worker = NUM_KEYS;

			/* Determine the number of keys in each bucket */
			for(i=ist_worker; i<iend_worker; i++)
				work_buff[key_array[i] >> shift]++;

			pthread_barrier_wait(&worker_barrier);

			/* Accumulative bucket sizes are the bucket pointers. */
			/* These are global sizes accumulated upon to each bucket */
			bucket_ptrs_worker[0] = 0;
			for( k=0; k< worker_id; k++ )  
				bucket_ptrs_worker[0] += bucket_size[k][0];

			for( i=1; i< NUM_BUCKETS; i++ ) { 
				bucket_ptrs_worker[i] = bucket_ptrs_worker[i-1];
				for( k=0; k< worker_id; k++ )
					bucket_ptrs_worker[i] += bucket_size[k][i];
				for( k=worker_id; k< num_workers; k++ )
					bucket_ptrs_worker[i] += bucket_size[k][i-1];
			}

			/* Sort into appropriate bucket */
			for(i=ist_worker; i<iend_worker; i++){
				k = key_array[i];
				key_buff2[bucket_ptrs_worker[k >> shift]++] = k;
			}

			/* The bucket pointers now point to the final accumulated sizes */
			if (worker_id == 0) {
				for( i=0; i< NUM_BUCKETS; i++ ){
		        	bucket_ptrs[i] = bucket_ptrs_worker[i];
					for( k=worker_id+1; k< num_workers; k++ )
						bucket_ptrs[i] += bucket_size[k][i];
				}
			}
		}
	});

	/* Now, buckets are sorted.  We only need to sort keys inside */
	/* each bucket, which can be done in parallel.  Because the distribution */
	/* of the number of keys in the buckets is Gaussian, the use of */
	/* a dynamic schedule should improve load balance, thus, performance */
	tbb::parallel_for(tbb::blocked_range<size_t>(0,NUM_BUCKETS),[&](const tbb::blocked_range<size_t>& r){
		INT_TYPE i;
		INT_TYPE m, k1, k2, k;

		for(i=r.begin(); i != r.end(); i++){
			/* Clear the work array section associated with each bucket */
			k1 = i * num_bucket_keys;
			k2 = k1 + num_bucket_keys;
			for ( k = k1; k < k2; k++ )
				key_buff_ptr[k] = 0;
			/* Ranking of all keys occurs in this section: */
			/* In this section, the keys themselves are used as their */
			/* own indexes to determine how many of each there are: their */
			/* individual population */
			m = (i > 0)? bucket_ptrs[i-1] : 0;
			for ( k = m; k < bucket_ptrs[i]; k++ )
				key_buff_ptr[key_buff_ptr2[k]]++; /* Now they have individual key population */
			/* To obtain ranks of each key, successively add the individual key */
			/* population, not forgetting to add m, the total of lesser keys, */
			/* to the first key population */
			key_buff_ptr[k1] += m;
			for ( k = k1+1; k < k2; k++ )
				key_buff_ptr[k] += key_buff_ptr[k-1];
		}
	});
#else /*USE_BUCKETS*/
	work_buff = key_buff1_aptr[myid];
	/* Clear the work array */
	for( i=0; i<MAX_KEY; i++ )
		work_buff[i] = 0;
	/* Ranking of all keys occurs in this section: */
	/* In this section, the keys themselves are used as their */
	/* own indexes to determine how many of each there are: their */
	/* individual population */
	for( i=0; i<NUM_KEYS; i++ )
		work_buff[key_buff_ptr2[i]]++; /* Now they have individual key population */
	/* To obtain ranks of each key, successively add the individual key population */
	for( i=0; i<MAX_KEY-1; i++ )   
		work_buff[i+1] += work_buff[i];
	/* Accumulate the global key population */
	for( k=1; k<num_procs; k++ ){
		for( i=0; i<MAX_KEY; i++ )
			key_buff_ptr[i] += key_buff1_aptr[k][i];
	}
#endif /*USE_BUCKETS*/  

	/* This is the partial verify test section */
	/* Observe that test_rank_array vals are */
	/* shifted differently for different cases */
	for( i=0; i<TEST_ARRAY_SIZE; i++ ){                                             
		k = partial_verify_vals[i]; /* test vals were put here */
		if( 0 < k  &&  k <= NUM_KEYS-1 )
		{
			INT_TYPE key_rank = key_buff_ptr[k-1];
			int failed = 0;

			switch( CLASS )
			{
				case 'S':
					if( i <= 2 )
					{
						if( key_rank != test_rank_array[i]+iteration )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-iteration )
							failed = 1;
						else
							passed_verification++;
					}
					break;
				case 'W':
					if( i < 2 )
					{
						if( key_rank != test_rank_array[i]+(iteration-2) )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-iteration )
							failed = 1;
						else
							passed_verification++;
					}
					break;
				case 'A':
					if( i <= 2 )
					{
						if( key_rank != test_rank_array[i]+(iteration-1) )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-(iteration-1) )
							failed = 1;
						else
							passed_verification++;
					}
					break;
				case 'B':
					if( i == 1 || i == 2 || i == 4 )
					{
						if( key_rank != test_rank_array[i]+iteration )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-iteration )
							failed = 1;
						else
							passed_verification++;
					}
					break;
				case 'C':
					if( i <= 2 )
					{
						if( key_rank != test_rank_array[i]+iteration )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-iteration )
							failed = 1;
						else
							passed_verification++;
					}
					break;
				case 'D':
					if( i < 2 )
					{
						if( key_rank != test_rank_array[i]+iteration )
							failed = 1;
						else
							passed_verification++;
					}
					else
					{
						if( key_rank != test_rank_array[i]-iteration )
							failed = 1;
						else
							passed_verification++;
					}
					break;
			}
			if( failed == 1 )
				printf( "Failed partial verification: "
						"iteration %d, test key %d\n", 
						iteration, (int)i );
		}
	}

	/* Make copies of rank info for use by full_verify: these variables */
	/* in rank are local; making them global slows down the code, probably */
	/* since they cannot be made register by compiler */
	if( iteration == MAX_ITERATIONS ) 
		key_buff_ptr_global = key_buff_ptr;
}
