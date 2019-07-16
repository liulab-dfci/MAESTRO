/*****************************************************************************
timer.h
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
#ifndef __TIMER_H__
#define __TIMER_H__

void start();

void stop();

unsigned long report();

struct timeval in();

unsigned long out(struct timeval i);

#endif
