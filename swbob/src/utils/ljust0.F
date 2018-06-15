      subroutine ljust0(chi     ,cho     )
      implicit none
C-----------------------------------------------------------------------
C
C ljust left justifies a character string (removes leading blanks)
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C
C------------------------------Arguments--------------------------------
C
C Input
C
      character*(*) chi       ! Input character string
C
C Output
C
      character*(*) cho       ! Returned left-justified string
C
C---------------------------Local variables-----------------------------
C
      integer i               ! Loop index
C
C-----------------------------------------------------------------------
C
C Find first nonblank character
C
      do i = 1,len( chi )
        if (chi(i:i).eq.' ') then
           cho(i:i)='0'
        else
           cho(i:i)=chi(i:i)
        endif
      end do
C
      return
      end
