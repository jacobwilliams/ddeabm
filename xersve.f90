!DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,&
         ICOUNT)
!***BEGIN PROLOGUE  XERSVE
!***SUBSIDIARY
!***PURPOSE  Record that an error has occurred.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
      IF (KFLAG<=0) THEN
!
!        Dump the table.
!
         IF (NMSG==0) RETURN
!
!        Print to each unit.
!
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT==0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
            WRITE (IUNIT,9000)
!
!           Print body of table.
!
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),&
                  NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
            IF (KOUNTX/=0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
         IF (KFLAG==0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB==LIBTAB(I) .AND. SUB==SUBTAB(I) .AND.&
               MES==MESTAB(I) .AND. NERR==NERTAB(I) .AND.&
               LEVEL==LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
!
         IF (NMSG<LENTAB) THEN
!
!           Empty slot found for new message.
!
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
!
!           Table is full.
!
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /&
         ' LIBRARY    SUBROUTINE MESSAGE START             NERR',&
         '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
