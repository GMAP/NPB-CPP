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

#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/syssgi.h>
#include <sys/immu.h>
#include <cerrno>
#include <cstdio>

/* the following works on SGI Power Challenge systems */
typedef unsigned long iotimer_t;
unsigned int cycleval;
volatile iotimer_t *iotimer_addr, base_counter;
double resolution;

/* address_t is an integer type big enough to hold an address */
typedef unsigned long address_t;

void timer_init(){  
	int fd;
	char *virt_addr;
	address_t phys_addr, page_offset, pagemask, pagebase_addr;

	pagemask = getpagesize() - 1;
	errno = 0;
	phys_addr = syssgi(SGI_QUERY_CYCLECNTR, &cycleval);
	if (errno != 0) {
		perror("SGI_QUERY_CYCLECNTR");
		exit(1);
	}
	/* rel_addr = page offset of physical address */
	page_offset = phys_addr & pagemask;
	pagebase_addr = phys_addr - page_offset;
	fd = open("/dev/mmem", O_RDONLY);

	virt_addr = mmap(0, pagemask, PROT_READ, MAP_PRIVATE, fd, pagebase_addr);
	virt_addr = virt_addr + page_offset;
	iotimer_addr = (iotimer_t *)virt_addr;
	/* cycleval in picoseconds to this gives resolution in seconds */
	resolution = 1.0e-12*cycleval; 
	base_counter = *iotimer_addr;
}

void wtime_(double *time){
	static int initialized = 0;
	volatile iotimer_t counter_value;
	if(!initialized){ 
		timer_init();
		initialized = 1;
	}
	counter_value = *iotimer_addr - base_counter;
	*time = (double)counter_value * resolution;
}

void wtime(double *time){
	static int initialized = 0;
	volatile iotimer_t counter_value;
	if (!initialized) { 
		timer_init();
		initialized = 1;
	}
	counter_value = *iotimer_addr - base_counter;
	*time = (double)counter_value * resolution;
}
