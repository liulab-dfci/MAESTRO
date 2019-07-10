/*****************************************************************************
timer.c
(c) 2012 - Ryan M. Layer
Hall Laboratory
Quinlan Laboratory
Department of Computer Science
Department of Biochemistry and Molecular Genetics
Department of Public Health Sciences and Center for Public Health Genomics,
University of Virginia
rl6sf@virginia.edu

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <sys/time.h>
#include "timer.h"

static struct timeval _start, _stop;

void start()
{
	gettimeofday(&_start,0);
}

void stop()
{
	gettimeofday(&_stop,0);
}

unsigned long report()
{
	return (_stop.tv_sec - _start.tv_sec) * 1000000 +  //seconds to microseconds
		_stop.tv_usec - _start.tv_usec;
}

struct timeval in()
{
	struct timeval i;
	gettimeofday(&i,0);
	return i;
}

unsigned long out(struct timeval i)
{
	struct timeval o;
	gettimeofday(&o,0);

	return (o.tv_sec - i.tv_sec) * 1000000 +  //seconds to microseconds
		o.tv_usec - i.tv_usec;
}
