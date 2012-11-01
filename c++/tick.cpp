#include "tick.h"

#include <iostream>
using namespace std;

static long tic_time = 0;
static long anti_tic_time = 0;

long tic()
{
	tic_time = getTickCount();
	return tic_time;
}

long toc(bool chatty)
{
	long dt = getTickCount() - tic_time;
	if (chatty) 
		cout << "Elapsed time: " << dt << "ms" << endl;
	return dt;
}

long getTickCount()
{
	timeval ts;
	gettimeofday(&ts, 0);
	return (long)(ts.tv_sec * 1000 + (ts.tv_usec / 1000));
}

void anti_tic()
{
	anti_tic_time = getTickCount();
}

void anti_toc()
{
	long anti_dt = getTickCount() - anti_tic_time;
	tic_time += anti_dt;
}

