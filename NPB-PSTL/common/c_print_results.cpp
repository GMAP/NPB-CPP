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

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>
*/ 

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <iomanip>
#include <iostream>

/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/
void c_print_results(std::string name,
		char class_npb,
		int n1, 
		int n2,
		int n3,
		int niter,
		double t,
		double mops,
		std::string optype,
		int passed_verification,
		std::string npbversion,
		std::string compiletime,
		std::string compilerversion,
		std::string cc,
		std::string clink,
		std::string c_lib,
		std::string c_inc,
		std::string cflags,
		std::string clinkflags,
		std::string rand){
	std::cout << "\n\n" << name << " Benchmark Completed" << std::endl;
	std::cout << " class_npb       =                        " << class_npb << std::endl;
	if((name[0]=='I')&&(name[1]=='S')){
		if(n3==0){
			long nn = n1;
			if(n2!=0){nn*=n2;}
			std::cout << " Size            =             "
				<< std::setw(12) << nn << std::endl;
		}else{
			std::cout << " Size            =           "
				<< std::setw(4) << n1 << "x"
				<< std::setw(4) << n2 << "x"
				<< std::setw(4) << n3
				<< std::endl;
		}
	}else{
		if((n2==0) && (n3==0)){
			if((name[0]=='E')&&(name[1]=='P')){
				std::cout << " Size            =          "
					<< std::setw(15) << pow(2.0, n1)
					<< std::endl;
			}else{
				std::cout << " Size            =             "
				<< std::setw(12) << n1
				<< std::endl;
			}
		}else{
			std::cout << " Size            =           "
				<< std::setw(4) << n1 << "x"
				<< std::setw(4) << n2 << "x"
				<< std::setw(4) << n3
				<< std::endl;
		}
	}	
	std::cout << " Iterations      =             "
		<< std::setw(12) << niter
		<< std::endl;
	std::cout << " Time in seconds =             "
		<< std::setw(12)
		<< std::setprecision(2) << std::fixed << t
		<< std::endl;
	std::cout << " Mop/s total     =             "
		<< std::setw(12)
		<< std::setprecision(2) << std::fixed << mops
		<< std::endl;
	std::cout << " Operation type  = "
		<< std::setw(24) << optype << std::endl;
	if(passed_verification < 0){
		std::cout << " Verification    =            NOT PERFORMED\n";
	}else if(passed_verification){
		std::cout << " Verification    =               SUCCESSFUL\n";
	}else{
		std::cout << " Verification    =             UNSUCCESSFUL\n";
	}
	std::cout << " Version         =             "
		<< std::setw(12) << npbversion
		<< std::endl;
	std::cout << " Compiler ver    =             "
		<< std::setw(12) << compilerversion
		<< std::endl;
	std::cout << " Compile date    =             "
		<< std::setw(12) << compiletime
		<< std::endl;
	std::cout << "\n Compile options:\n";
	std::cout << "    CC           = " << cc << std::endl;
	std::cout << "    CLINK        = " << clink << std::endl;
	std::cout << "    C_LIB        = " << c_lib << std::endl;
	std::cout << "    C_INC        = " << c_inc << std::endl;
	std::cout << "    CFLAGS       = " << cflags << std::endl;
	std::cout << "    CLINKFLAGS   = " << clinkflags << std::endl;
	std::cout << "    RAND         = " << rand << std::endl;
#ifdef SMP
	evalue = getenv("MP_SET_NUMTHREADS");
	printf("   MULTICPUS = %s\n", evalue);
#endif    
	std::cout << "\n\n";
	
	std::cout << "----------------------------------------------------------------------\n";
	std::cout << "    NPB-CPP is developed by: \n";
	std::cout << "        Dalvan Griebler\n";
	std::cout << "        Gabriell Araujo (Sequential Porting)\n";
	std::cout << "        Júnior Löff (Parallel Implementation)\n";
	std::cout << "\n";
	std::cout << "    In case of questions or problems, please send an e-mail to us:\n";	
	std::cout << "        dalvan.griebler; gabriell.araujo; junior.loff@edu.pucrs.br\n";
	std::cout << "----------------------------------------------------------------------\n";
	std::cout << "\n";
}
