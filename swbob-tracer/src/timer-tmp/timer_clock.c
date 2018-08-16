/**********************************************************************
 *
 *                           Copyright (C) 1995   
 *
 *             University Corporation for Atmospheric Research
 *
 *                          All rights reserved.
 *
 *     No part of this  code may  be reproduced,  stored in  a retrieval
 *     system,  translated,  transcribed,  transmitted or distributed in
 *     any form or by any means, manual, electric, electronic, chemical,
 *     electo-magnetic, mechanical, optical, photocopying, recording, or
 *     otherwise, without the prior  explicit written permission  of the
 *     author(s)  or their  designated proxies.  In no event  shall this
 *     copyright notice be removed or altered in any way.
 *
 *     This code is provided "as is",  without any warranty of any kind,
 *     either  expressed or implied,  including but not limited to,  any
 *     implied warranty of  merchantibility or fitness  for any purpose.
 *     In no event will any party who distributed the code be liable for
 *     damages or for any claim(s) by any other party, including but not
 *     limited to,  any lost  profits,  lost monies,  lost data  or data
 *     rendered inaccurate,  losses sustained  by third parties,  or any
 *     other special, incidental or consequential damages arising out of
 *     the use or inability to use the program,  even if the possibility
 *     of such damages has  been advised against.  The entire risk as to
 *     the quality,  the performance, and the fitness of the program for
 *     any particular purpose lies with the party using the code.
 *
 *     *****************************************************************
 *     ANY USE OF  THIS CODE CONSTITUES  ACCEPTANCE OF  THE TERMS OF THE
 *                              ABOVE STATEMENTS
 *     *****************************************************************
 *
***********************************************************************/

/*
 *	timer_clock.c
 *
 *	this file contains vendor defined timer routines that returns
 *      a 64 bit floating point clock time.
 *
 *
 */

#if defined(HP)

#include <time.h>
#include <stdio.h>
#include <math.h>
#if defined(MPI)
#include "mpi.h"
#endif

/*
 *   Hewlett-Packard timer routine
 */
double timer_clock()
{
  double t;
#if defined(MPI)
  t = MPI_Wtime();
  printf(" timer_clock (debug): t = %f\n",t);
#else
  struct timeval buffer;
  struct timezone dummy;

  gettimeofday (&buffer, &dummy);
  t = (double)buffer.tv_sec + ((double)buffer.tv_usec*1.0e-6);
#endif

  return (t);

}

#endif

#if defined (SunOS)

/*
 *   SUN nanosecond timer routine
 */
#include <sys/time.h>

static int first_call = 1;

double timer_clock_(){
   hrtime_t nsec;

/*
   if (first_call) {
      first_call = 0;
      init_ecache_();
   }
*/

   nsec = gethrtime();
   return(((double)nsec)*1.0e-09);
}

#endif

#if defined(Linux)

/*
 *   Linux timer routine
 */

#include <linux/time.h>

double timer_clock_()
{
  double t;

  struct timeval buffer;
  struct timezone dummy;

  gettimeofday (&buffer, &dummy);
  t = (double)buffer.tv_sec + ((double)buffer.tv_usec*1.0e-6);

  return (t);
}
#endif

#if defined(AIX)
extern double rtc();
static double start = 0.0;
double timer_clock()
{
        extern double start;
        if ( start != 0.0 )
                return (rtc() - start);
        else {
                start = rtc();
                return 0.0;
        }
}

#endif

#if defined(IRIX64)

#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

struct rusage
   rusage1, rusage2;

double timer_clock_()
{

/*
 *	get some the user and system time statistics and
 *       sum them
 */

   getrusage(RUSAGE_SELF, &rusage2);

   return( (double) rusage2.ru_utime.tv_sec +
          ((double) rusage2.ru_utime.tv_usec) * 1.e-6 +
           (double) rusage2.ru_stime.tv_sec +
          ((double) rusage2.ru_stime.tv_usec) * 1.e-6);
   }

#endif

#if defined(CRAY)

#if defined(T3D)
#define CLK_TCK 151499000
#endif

/*
 *  Cray timer routine
 */

#include <time.h>
double TIMER_CLOCK()
{
	return((double)rtclock()/CLK_TCK);
}

#endif

#if defined(HP_EXEMPLAR)
/*
 * A lightweight and thread-safe procedure to return the
 * wall-clock time on the HP/Convex Exemplar with
 * microsecond resolution.
 *
 * timer_clock usage is:
 *
 *   double begin_time , end_time , elapsed_time;
 *   double timer_clock();
 *
 *   begin_time = timer_clock();
 *      < code to be timed >
 *   end_time = timer_clock();
 *   elapsed_time = end_time - begin_time;
 *
 */

#include <spp_prog_model.h>
#include <sys/types.h>
#include <sys/cnx_ail.h>

#define SECS_PER_MICROSEC 0.000001
#define TWO_TO_32 4294967296.0

double timer_clock()
{
   cnx_toc_t temp_toc;
   unsigned int *ptr;
   double temp_secs;

   ptr = (unsigned int *) &temp_toc;

/*
 * Read the time-of-century clock to get a 64-bit integer.
 */

   temp_toc = toc_read();

/*
 * Convert the 64-bit integer to a double that
 * represents the absolute time since the epoch.
 */

   temp_secs  = ( (double) ptr[0] ) * TWO_TO_32;
   temp_secs += ( (double) ptr[1] );
   temp_secs *= SECS_PER_MICROSEC;

   return( temp_secs );
}
#endif


