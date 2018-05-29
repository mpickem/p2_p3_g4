program bandspin
  implicit none

  integer            :: nbands = 5
  integer, parameter :: spins = 2

  integer            :: i1,i2,i3,i4,j1,j2,j3,j4
  integer            :: ind
  integer            :: cnt

  logical            :: ch_pp = .true.
  logical            :: int_d = .true.
  character(len=150) :: tmp_string


  if (iargc() .ne. 1) then
    stop 'program has to be executed with 1 arguments -- bands'
  endif

  call getarg(1,tmp_string)
  read(tmp_string,'(I1)') nbands

  open(unit=123, file='channels.dat')
  open(unit=124, file='pythonlist.dat')
  write(124,'("mylist=[")', advance='no')

  cnt = 0

  do i1=1,nbands
  do i2=1,nbands
  do i3=1,nbands
  do i4=1,nbands
    do j1=1,spins
    do j2=1,spins
    do j3=1,spins
    do j4=1,spins


! SU2
      if ((j1 .eq. j2) .and. (j3 .eq. j4)) goto 100 ! sigma sigma'
      if ((j1 .eq. j4) .and. (j2 .eq. j3)) goto 100 ! sigma sigma' quer
      cycle
100 continue

! Kanamori
      if (i1.eq. i2 .and. i3 .eq. i4)  then
          goto 101
      endif
      if (i1.eq. i3 .and. i2 .eq. i4) then
        goto 101
      endif
      if (i1.eq. i4 .and. i2 .eq. i3) then
        goto 101
      endif
      cycle
101 continue

      call component2index(nbands, ind, i1, j1, i2, j2, i3, j3, i4, j4)

      if (ind .ne. 0) then
        write(123,'(I5.5,5X,A4,4I4,5X,A4,4I4)') ind, '----', i1, i2, i3, i4, '----', j1, j2, j3, j4
        ! if (mod(cnt,10)==9) then
        !   write(124,'(I5.5,", ")', advance='no') ind
        if (mod(cnt,10)==0 .and. cnt>0) then
          write(124,'(" \")')
          write(124,'(8X,I5.5,", ")', advance='no') ind
        else if(ind < 2**4*nbands**4) then
          write(124,'(I5.5,", ")', advance='no') ind
        else
          write(124,'(I5.5,"]")', advance='no') ind
        endif

        cnt = cnt+1
      endif

    enddo
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  enddo

  close(123)

  contains

  subroutine component2index(Nbands, ind, b1, s1, b2, s2, b3, s3, b4, s4)
    implicit none

    integer,intent(in) :: Nbands
    integer,intent(in) :: b1, s1, b2, s2, b3, s3, b4, s4
    integer,intent(inout) :: ind
    integer :: g1, g2, g3, g4

    g1=2*(b1-1) + s1
    g2=2*(b2-1) + s2
    g3=2*(b3-1) + s3
    g4=2*(b4-1) + s4

    ind =  8*Nbands**3*(g1-1) + 4*Nbands**2*(g2-1) + 2*Nbands*(g3-1) + g4

  end subroutine component2index

end program
