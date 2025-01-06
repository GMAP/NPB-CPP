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
	P. O. Frederickson
	D. H. Bailey
	A. C. Woo

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
#include <execution>

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

struct Returnable{
	double sx;
	double sy;
	std::vector<double> q;
	
	Returnable() : sx(0.0), sy(0.0), q(std::vector<double>(NQ)) {}
	Returnable(const double &sx, const double &sy, const std::vector<double> &q) : sx(sx), sy(sy), q(q) {}

	Returnable operator+(Returnable other) const {
		std::transform(q.begin(), q.end(), 	other.q.begin(), other.q.begin(), std::plus<double>());
		other.sx += sx;
		other.sy += sy;
		return other;
	}
};

auto parallel_policy = std::execution::par;

static std::vector<double> x(NK_PLUS);
static std::vector<double> q(NQ);

/* ep */
int main(int argc, char **argv){
	double  Mops, t3, t4, x1, x2;
	double t1, t2;
	double  sx, sy, tm, an, tt, gc;
	double  sx_verify_value, sy_verify_value, sx_err, sy_err;
	int     i, ik, kk, l, k;
	boolean verified, timers_enabled;

	timers_enabled = std::filesystem::exists("timer.flag");

	std::cout << std::endl << std::endl;
	std::cout << " NAS Parallel Benchmarks 4.1 Serial C++ version - EP Benchmark";
	std::cout << std::endl << std::endl;
	std::cout 
		<< " Number of random numbers generated: "
		<< std::setprecision(0) << std::fixed
		<< pow(2.0, M+1)
		<< std::endl;

	timer_clear(0);
	timer_clear(1);
	timer_clear(2);
	timer_start(0);

	t1 = A;
	vranlc(0, t1, A, x.data());

	/* compute AN = A ^ (2 * NK) (mod 2^46) */
	t1 = A;

	for(i=0; i<=MK; i++) {
		t2 = randlc(t1, t1);
	}

	an = t1;

	/*
	 * each instance of this loop may be performed independently. we compute
	 * the k offsets separately to take into account the fact that some nodes
	 * have more numbers to generate than others
	 */

	CountIterator iter(NN);
	Returnable result = std::transform_reduce(parallel_policy, iter.front(), iter.front()+NN, Returnable(), std::plus<Returnable>(), [&] (int index) {
		double t1, t2;
    	std::vector<std::pair<double, double>> pairs_x(NK_PLUS/2);
		std::vector<double> q(NQ);
		double sx=0.0, sy=0.0;
		double t3, t4, x1, x2;
		int kk, ik, l;

		kk = index;
		t1 = S;
		t2 = an;

		auto breakPointEl = std::find_if(iter.front()+1, iter.front()+100, [&kk] (int el) {kk/=2;return kk==0;});
		kk=index;
		
		std::for_each_n(iter.front()+1, *breakPointEl+1, [&] (int i) {
			ik = kk / 2;
			if((2 * ik) != kk)
				t3=randlc(t1,t2);

			t3 = randlc(t2,t2);
			kk = ik;
		});

		/* compute uniform pseudorandom numbers */
		if(timers_enabled) timer_start(2);
		vranlc(2*NK, t1, A, (double*)pairs_x.data());
    	if(timers_enabled) timer_stop(2);

		/*
		 * compute gaussian deviates by acceptance-rejection method and
		 * tally counts in concentric square annuli. this loop is not
		 * vectorizable.
		 */

		if(timers_enabled) timer_start(1);

		/* 
		 * Iterate through the pairs of x vector
		*/
		std::for_each(pairs_x.begin(), pairs_x.end(), [&] (std::pair<double,double>& p) {
			x1 = 2.0 * p.first - 1.0;
			x2 = 2.0 * p.second - 1.0;
			t1 = pow2(x1) + pow2(x2);
			if (t1 <= 1) {
				t2 = sqrt(-2.0 * log(t1) / t1);
				t3 = (x1 * (t2));
				t4 = (x2 * (t2));
				l = std::max(fabs(t3), fabs(t4));
				q.at(l) += 1.0;
				sx += t3;
				sy += t4;
			}
		});

		if(timers_enabled) timer_stop(1);
		return Returnable(sx, sy, q);
	});
	q  = result.q;
	sx = result.sx;
	sy = result.sy;

	gc = std::accumulate(q.begin(), q.end(), 0);

	timer_stop(0);
	tm = timer_read(0);

	verified = TRUE;
	switch (M){
		case 24:
			sx_verify_value = -3.247834652034740e+3;
			sy_verify_value = -6.958407078382297e+3;
			break;
		case 25:
			sx_verify_value = -2.863319731645753e+3;
			sy_verify_value = -6.320053679109499e+3;
			break;
		case 28:
			sx_verify_value = -4.295875165629892e+3;
			sy_verify_value = -1.580732573678431e+4;
			break;
		case 30:
			sx_verify_value =  4.033815542441498e+4;
			sy_verify_value = -2.660669192809235e+4;
			break;
		case 32:
			sx_verify_value =  4.764367927995374e+4;
			sy_verify_value = -8.084072988043731e+4;
			break;
		case 36:
			sx_verify_value =  1.982481200946593e+5;
			sy_verify_value = -1.020596636361769e+5;
			break;
		case 40:
			sx_verify_value = -5.319717441530e+05;
			sy_verify_value = -3.688834557731e+05;
			break;
		default:
			verified = FALSE;
	}

	if(verified){
		sx_err = fabs((sx - sx_verify_value) / sx_verify_value);
		sy_err = fabs((sy - sy_verify_value) / sy_verify_value);
		verified = ((sx_err <= EPSILON) && (sy_err <= EPSILON));
	}
	
	Mops = pow(2.0, M+1)/tm/1000000.0;

	std::cout << std::endl << " EP Benchmark Results:";
	std::cout << std::endl << std::endl;

	std::cout 
		<< " CPU Time ="
		<< std::setw(10) << std::setprecision(4)
		<< std::fixed << tm << std::endl;

	std::cout << " N = 2^" << std::setw(5) << M << std::endl;

	std::cout 
		<< " No. Gaussian Pairs ="
		<< std::setw(10) << std::setprecision(0)
		<< std::fixed << gc << std::endl;

	std::cout 
		<< " Sums = "
		<< std::setw(25) << std::setprecision(15) << std::scientific
		<< sx << std::setw(25) << sy << std::endl;

	std::cout << " Counts:" << std::endl;
	for(int i=0; i<NQ-1; i++){
		std::cout
			<< std::setw(3) << (int) i
			<< std::setw(15) << std::setprecision(0)
			<< std::fixed << q.at(i) << std::endl;
	}

	c_print_results(
		"EP",
		CLASS,
		M+1,
		0,
		0,
		0,
		tm,
		Mops,
		"Random numbers generated",
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
		CS7
	);

	if(timers_enabled){
		if(tm <= 0.0)
			tm = 1.0;

		tt = timer_read(0);
		std::cout 
			<< "Total time:" << std::setw(14)
			<< std::setprecision(3) << std::fixed << tt
			<< "(" << std::setw(6)
			<< std::setprecision(2) << std::fixed << (tt*100.0/tm)
			<< ")" << std::endl;

		tt = timer_read(1);
		std::cout 
			<< "Gaussian pairs:" << std::setw(10)
			<< std::setprecision(3) << std::fixed << tt
			<< "(" << std::setw(6)
			<< std::setprecision(2) << std::fixed << (tt*100.0/tm)
			<< ")" << std::endl;

		tt = timer_read(2);
		std::cout 
			<< "Random numbers:" << std::setw(10)
			<< std::setprecision(3) << std::fixed << tt
			<< "(" << std::setw(6)
			<< std::setprecision(2) << std::fixed << (tt*100.0/tm)
			<< ")" << std::endl;
	}

	return 0;
}
