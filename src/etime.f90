! --------------------------------------------------------------------------
! dummy function ETIME to be used with compilers like ICF that don't have
! a useful ETIME definition
! --------------------------------------------------------------------------

REAL(4) FUNCTION ETIME(TIME)
  REAL(4) TIME, MyTIME
  SAVE MyTIME
  MyTIME= MyTIME+1.0
  ETIME= MyTIME
  RETURN
END
