#ifndef TICK_H
#define TICK_H

#include <sys/time.h>

long getTickCount();
long tic();
long toc(bool chatty = false);
void anti_tic();
void anti_toc();

/*
static inline void clflush(volatile void *__p)
{
	asm volatile("clflush %0" : "+m" (*(char __force *)__p));
}
*/

#endif
		
