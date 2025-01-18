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
	
--------------------------------------------------------------------------

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>
	Arthur S. Bianchessi <arthur.bianchessi@edu.pucrs.br>
	Leonardo Mallmann <leonardo.mallmann@edu.pucrs.br>
*/


#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"
#include <barrier>

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
int kx2_const, kx2;
int num_procs = std::thread::hardware_concurrency();
std::barrier barrier(num_procs);
auto policy = std::execution::par;                           

/************************************/
/* These are the three main arrays. */
/* See SIZE_OF_BUFFERS def above    */
/************************************/
std::vector<INT_TYPE> key_array(SIZE_OF_BUFFERS);
std::vector<INT_TYPE> key_buff1(MAX_KEY);
std::vector<INT_TYPE> key_buff2(SIZE_OF_BUFFERS, 0);
std::vector<INT_TYPE> partial_verify_vals(TEST_ARRAY_SIZE);
std::vector<INT_TYPE> bucket_ptrs(NUM_BUCKETS);
std::vector<std::vector<INT_TYPE>> bucket_size(num_procs, std::vector<INT_TYPE>(NUM_BUCKETS));
std::vector<INT_TYPE> key_buff1_aptr(MAX_KEY*num_procs);

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
CountIterator iter(NUM_KEYS);

/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/
int main(int argc, char** argv){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	std::cout << " DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on" << std::endl;
#endif
	int timer_on;
	double timecounter;

	/* Initialize timers */
    timer_on = std::filesystem::exists("timer.flag");
	timer_clear( T_BENCHMARKING );
	if(timer_on){
		timer_clear( T_INITIALIZATION );
		timer_clear( T_SORTING );
		timer_clear( T_TOTAL_EXECUTION );
	}

	if(timer_on) timer_start( T_TOTAL_EXECUTION );

	/* Initialize the verification arrays if a valid class */
	for(int i=0; i<TEST_ARRAY_SIZE; i++ ) {
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
	}

	/* Printout initial NPB info */
	std::cout 
        << std::endl << std::endl 
        <<" NAS Parallel Benchmarks 4.1 Serial C++ version - IS Benchmark"
        << std::endl << std::endl;

	std::cout 
        << " Size:  " << (long)TOTAL_KEYS 
        << "  (class " <<  CLASS << ")" 
        << std::endl;

	std::cout << " Iterations:  " << MAX_ITERATIONS << std::endl << std::endl;

	if(timer_on)timer_start( T_INITIALIZATION );

	/* Generate random number sequence and subsequent keys on all procs */
	create_seq(314159265.00 /* Random number gen seed */,
			1220703125.00 /* Random number gen mult */);
	
	if(timer_on)timer_stop( T_INITIALIZATION );

	/* Do one interation for free (i.e., untimed) to guarantee initialization of */
	/* all data and code pages and respective tables */
	rank( 1 );

	/* Start verification counter */
	passed_verification = 0;

	if( CLASS != 'S' ) std::cout << std::endl << "iteration" << std::endl;

	/* Start timer */             
	timer_start( T_BENCHMARKING );

	/* This is the main iteration */
	for(int iteration=1; iteration<=MAX_ITERATIONS; iteration++){
		if(CLASS != 'S') std::cout << "\t\t" << iteration << std::endl;
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
	c_print_results(
        "IS",
		CLASS,
		(int)(TOTAL_KEYS/64),
		64,
		0,
		MAX_ITERATIONS,
		timecounter,
		((double)(MAX_ITERATIONS*TOTAL_KEYS))/timecounter/1000000.0,
		"keys ranked",
		passed_verification,
		NPBVERSION,
		COMPILETIME,
		COMPILERVERSION,
		CS1,
		CS2,
		CS3,
		CS4,
		CS5,
		CS6,
		CS7
    );

	/* Print additional timers */
	if(timer_on){
		double t_total, t_percent;

		t_total = timer_read( T_TOTAL_EXECUTION );
		std::cout << std::endl << "Additional timers -" << std::endl;
		std::cout << " Total execution:    " << std::setprecision(3) << t_total << std::endl;
		if (t_total == 0.0) t_total = 1.0;

		timecounter = timer_read(T_INITIALIZATION);
		t_percent = timecounter/t_total * 100.;
        std::cout 
            << " Initialization :    " 
            << std::setprecision(3) << timecounter << " (" 
            << std::setprecision(2) << t_percent << "%)"
            << std::endl;

		timecounter = timer_read(T_BENCHMARKING);
		t_percent = timecounter/t_total * 100.;
        std::cout 
            << " Benchmarking   :    " 
            << std::setprecision(3) << timecounter << " (" 
            << std::setprecision(2) << t_percent << "%)"
            << std::endl;

		timecounter = timer_read(T_SORTING);
		t_percent = timecounter/t_total * 100.;
        std::cout 
            << " Sorting        :    " 
            << std::setprecision(3) << timecounter << " (" 
            << std::setprecision(2) << t_percent << "%)"
            << std::endl;
	}

	return 0;
}

/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/
void create_seq(double seed,
		double a){
	std::for_each_n(policy, iter.front(), num_procs, [&](int myid){
		double x, s;
		INT_TYPE i, k;

		INT_TYPE k1, k2;
		double an = a;
		INT_TYPE mq;

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

		for(int i=k1; i<k2; i++){
			x = randlc(s, an);
			x += randlc(s, an);
			x += randlc(s, an);
			x += randlc(s, an);
			key_array[i] = k*x;
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
			(void)randlc( t2, t2 );
			kk = ik;
		}
		else{
			(void)randlc( t1, t2 );
			kk = kk - 1;
		}
	}
	(void)randlc( t1, t2 );

	return( t1 );
}

/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/
void full_verify(){
	INT_TYPE j;
	INT_TYPE k, k1, k2;

	/* Now, finally, sort the keys: */
	/* Copy keys into work array; keys in key_array will be reassigned. */

#ifdef USE_BUCKETS
	/* Buckets are already sorted. Sorting keys within each bucket */
	std::for_each_n(policy, iter.front(), NUM_BUCKETS, [&](int j){
		INT_TYPE k,k1;
		k1 = (j > 0) ? bucket_ptrs[j-1] : 0;
		for (int i = k1; i < bucket_ptrs[j]; i++ ) {
			k = --key_buff_ptr_global[key_buff2[i]];
			key_array[k] = key_buff2[i];
		}
	});
#else    
	int myid = 0, num_procs = 1;
	std::copy ( key_array.begin(), key_array.begin()+NUM_KEYS, key_buff2.begin() );

	/* This is actual sorting. Each thread is responsible for a subset of key values */
	j = num_procs;
	j = (MAX_KEY + j - 1) / j;
	k1 = j * myid;
	k2 = k1 + j;
	if (k2 > MAX_KEY) k2 = MAX_KEY;
	for(int i=0; i<NUM_KEYS; i++ ) {
		if (key_buff2[i] >= k1 && key_buff2[i] < k2) {
			k = --key_buff_ptr_global[key_buff2[i]];
			key_array[k] = key_buff2[i];
		}
	}
#endif

	/* Confirm keys correctly sorted: count incorrectly sorted keys, if any */
	j = std::count_if(policy, iter.front()+1, iter.front()+NUM_KEYS, [&](int i){return key_array[i-1] > key_array[i];});

	if( j != 0 )
		std::cout << "Full_verify: number of keys out of sort: " << (long)j << std::endl;
	else
		passed_verification++;
}


/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/
void rank(int iteration){

	INT_TYPE k;
	INT_TYPE *key_buff_ptr, *key_buff_ptr2;

#ifdef USE_BUCKETS
	int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
	INT_TYPE num_bucket_keys = (1L << shift);
#endif

	key_array[iteration] = iteration;
	key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;

	/* Determine where the partial verify test keys are, load into */
	/* top of array bucket_size */
	for(int i=0; i<TEST_ARRAY_SIZE; i++ ) {
		partial_verify_vals[i] = key_array[test_index_array[i]];
	}

	/* Setup pointers to key buffers */
#ifdef USE_BUCKETS
	key_buff_ptr2 = key_buff2.data();
#else
	key_buff_ptr2 = key_array.data();
#endif
	key_buff_ptr = key_buff1.data();

	INT_TYPE m, k1, k2;

	/* Bucket sort is known to improve cache performance on some */
	/* cache based systems.  But the actual performance may depend */
	/* on cache size, problem size. */
#ifdef USE_BUCKETS

	/* Determine the number of keys in each bucket */
	std::for_each_n(policy, iter.front(), num_procs, [&](int myid){

		std::fill(bucket_size[myid].begin(), bucket_size[myid].end(), 0);
		
		int chunk = NUM_KEYS/num_procs;
		int residual = NUM_KEYS%num_procs;

		int ist_worker = myid*chunk + ((residual>1 && (myid<residual)) ? myid : residual);
		int iend_worker = myid*chunk + chunk + ((residual>1 && (myid<residual)) ? myid+1 : residual);

		if(myid==0) ist_worker = 0;
		if(myid==num_procs-1) iend_worker = NUM_KEYS;

		for(int i=ist_worker; i<iend_worker; i++) {
			bucket_size[myid][key_array[i] >> shift]++;
		}

		barrier.arrive_and_wait();

		INT_TYPE m, k1, k2;
		int kx2;
		INT_TYPE i, k, j;
		std::vector<INT_TYPE> bucket_ptrs_worker(NUM_BUCKETS);

		/* Accumulative bucket sizes are the bucket pointers. */
		/* These are global sizes accumulated upon to each bucket */
		bucket_ptrs_worker[0] = 0;
		for(int k=0; k< myid; k++ )  {
			bucket_ptrs_worker[0] += bucket_size[k][0];
		}

		for(int i=1; i< NUM_BUCKETS; i++ ) { 
			bucket_ptrs_worker[i] = bucket_ptrs_worker[i-1];
			for(int k=0; k< myid; k++ ) {
				bucket_ptrs_worker[i] += bucket_size[k][i];
			}
			for(int k=myid; k< num_procs; k++ ) {
				bucket_ptrs_worker[i] += bucket_size[k][i-1];
			}
		}


		/* Sort into appropriate bucket */
		for(int i=ist_worker; i<iend_worker; i++){
			k = key_array[i];
			key_buff2[bucket_ptrs_worker[k >> shift]++] = k;
		}
		/* The bucket pointers now point to the final accumulated sizes */
		if (myid == 0) {
			for(int i=0; i< NUM_BUCKETS; i++ ){
		        bucket_ptrs[i] = bucket_ptrs_worker[i];
				for(int k=myid+1; k< num_procs; k++ ){
					bucket_ptrs[i] += bucket_size[k][i];
				}
			}
		}
	});
	/* Now, buckets are sorted.  We only need to sort keys inside */
	/* each bucket, which can be done in parallel.  Because the distribution */
	/* of the number of keys in the buckets is Gaussian, the use of */
	/* a dynamic schedule should improve load balance, thus, performance */
	std::for_each_n(policy, iter.front(), NUM_BUCKETS, [&](int i){
		INT_TYPE m, k1, k2, k;
		/* Clear the work array section associated with each bucket */
		k1 = i * num_bucket_keys;
		k2 = k1 + num_bucket_keys;
		std::fill(key_buff_ptr+k1, key_buff_ptr+k2, 0);
		
		/* Ranking of all keys occurs in this section: */
		/* In this section, the keys themselves are used as their */
		/* own indexes to determine how many of each there are: their */
		/* individual population */
		m = (i > 0)? bucket_ptrs[i-1] : 0;
		for (int k = m; k < bucket_ptrs[i]; k++ ) {
			key_buff_ptr[key_buff_ptr2[k]]++;
		}
		/* Now they have individual key population */
		/* To obtain ranks of each key, successively add the individual key */
		/* population, not forgetting to add m, the total of lesser keys, */
		/* to the first key population */
		key_buff_ptr[k1] += m;
		for (int k = k1+1; k < k2; k++ ) {
			key_buff_ptr[k] += key_buff_ptr[k-1];
		}
	});
#else /*USE_BUCKETS*/
	int myid = 0, num_procs = 1;
	work_buff = &key_buff1_aptr[myid];
	/* Clear the work array */
	std::fill(work_buff, work_buff+MAX_KEY, 0);
	
	/* Ranking of all keys occurs in this section: */
	/* In this section, the keys themselves are used as their */
	/* own indexes to determine how many of each there are: their */
	/* individual population */
	std::for_each(key_buff_ptr2, key_buff_ptr2+NUM_KEYS, [&](int k){work_buff[k]++;});
	/* Now they have individual key population */
	
	/* To obtain ranks of each key, successively add the individual key population */
	for(int i=0; i<MAX_KEY-1; i++ ) {
		work_buff[i+1] += work_buff[i];
	}

	/* Accumulate the global key population */
	for(int k=1; k<num_procs; k++ ){
		kx2 = k*MAX_KEY;
		for(int i=0; i<MAX_KEY; i++ ) {
			key_buff_ptr[i] += key_buff1_aptr[kx2+i];
		}
	}
#endif /*USE_BUCKETS*/  

	/* This is the partial verify test section */
	/* Observe that test_rank_array vals are */
	/* shifted differently for different cases */
	for(int i=0; i<TEST_ARRAY_SIZE; i++ ) {                                           
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
				std::cout
                    << "Failed partial verification: "
                    << "iteration" << iteration 
                    << ", test key " << (int)i 
                    << std::endl;
		}
	}

	/* Make copies of rank info for use by full_verify: these variables */
	/* in rank are local; making them global slows down the code, probably */
	/* since they cannot be made register by compiler */
	if( iteration == MAX_ITERATIONS ) 
		key_buff_ptr_global = key_buff_ptr;
}