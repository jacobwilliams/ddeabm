!DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!***BEGIN PROLOGUE  XERMSG
!***PURPOSE  Process error messages for SLATEC and other libraries.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR /= 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
      IF (NERR<-9999999 .OR. NERR>99999999 .OR. NERR==0 .OR.&
         LEVEL<-1 .OR. LEVEL>2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //&
            'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//&
            'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
!
!       RECORD THE MESSAGE.
!
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
      IF (LEVEL==-1 .AND. KOUNT>1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
      IF (LEVEL<2 .AND. LKNTRL==0) GO TO 30
      IF (LEVEL==0 .AND. KOUNT>MAXMES) GO TO 30
      IF (LEVEL==1 .AND. KOUNT>MAXMES .AND. MKNTRL==1) GO TO 30
      IF (LEVEL==2 .AND. KOUNT>MAX(1,MAXMES)) GO TO 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!       IS NOT ZERO.
!
      IF (LKNTRL /= 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
      IF (LKNTRL > 0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
         IF (LEVEL <= 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL == 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
         IF ((MKNTRL==2 .AND. LEVEL>=1) .OR.&
             (MKNTRL==1 .AND. LEVEL==2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
         IF (LKNTRL > 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       NOW SEND OUT THE MESSAGE.
!
      CALL XERPRN (' *  ', -1, MESSG, 72)
!
!       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
      IF (LKNTRL > 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) /= ' ') GO TO 20
   10    CONTINUE
!
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
!
!       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
      IF (LKNTRL /= 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
!
!       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 IF (LEVEL<=0 .OR. (LEVEL==1 .AND. MKNTRL<=1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
      IF (LKNTRL>0 .AND. KOUNT<MAX(1,MAXMES)) THEN
         IF (LEVEL == 1) THEN
            CALL XERPRN&
               (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
