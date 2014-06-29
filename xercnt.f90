!DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!***BEGIN PROLOGUE  XERCNT
!***SUBSIDIARY
!***PURPOSE  Allow user control over handling of errors.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCNT-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
!***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
