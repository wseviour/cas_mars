#include <params.h>

/* error_dup.c 
 * redirect standard error from process into a file 
 * also redirect standard output to a file 
 */

#include <stdio.h>
#include <fcntl.h>

#define PP_ONE_FOR_ALL_ERROR
#define STANDARD_ERROR 2

#define STANDARD_OUTPUT 1

void error_dup(int *me)
{
    int newfd ;
    char filename[256] ;

/* redirect standard out*/
    if ( *me == 0 )
    {
        sprintf(filename,"mpi.out.%03d",*me) ;
    }
    else
    {
        sprintf(filename,"/dev/null") ; 
    }
    if ((newfd = open( filename, O_CREAT | O_WRONLY, 0666 )) < 0 )
    {
	perror("error_dup: cannot open mpi.out.xxx; ") ;
	fprintf(stderr,"sending error to standard error and continuing.\n") ;
	return ;
    }
    if(*me!=0)
    if( dup2( newfd, STANDARD_OUTPUT ) < 0 )
    {
	perror("error_dup: dup2 fails to change output descriptor; ") ;
	fprintf(stderr,"sending error to standard error and continuing.\n") ;
	close(newfd) ;
	return ;
    }

/* redirect standard error */
#ifdef PP_ONE_FOR_ALL_ERROR
    if ( *me == 0 )
    {
        sprintf(filename,"mpi.err.%03d",*me) ;
    }
    else
    {
        sprintf(filename,"/dev/null") ;
    }
#else
        sprintf(filename,"mpi.err.%03d",*me) ;
#endif
    if ((newfd = open( filename, O_CREAT | O_WRONLY, 0666 )) < 0 )
    {
	perror("error_dup: cannot open mpi.err; ") ;
	fprintf(stderr,"sending error to standard error and continuing.\n") ;
	return ;
    }
    if( dup2( newfd, STANDARD_ERROR ) < 0 )
    {
	perror("error_dup: dup2 fails to change error descriptor; ") ;
	fprintf(stderr,"sending error to standard error and continuing.\n") ;
	close(newfd) ;
	return ;
    }

}

void error_dup_(int *me){
    error_dup(me);
}
