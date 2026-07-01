! ===========================================================================
! msd-com.f90 — Mean Squared Displacement resolved by z-bin
! ===========================================================================
!
! PURPOSE:
!   Computes the lateral (xy-plane) MSD as a function of time, binned by
!   each atom's z-position. Uses non-overlapping sliding windows so that
!   the result can be averaged over independent time origins.
!
! expects unifies dump format
!   dump  d1 all custom 1000 traj.lammpstrj id type mol xu yu zu
!   dump_modify d1 sort id
!   Columns: id  type  mol  xu  yu  zu   (unwrapped coords — required for MSD)
!
! INPUT FILE:
!   path to LAMMPS dump file
!   path for output file
!   nbins     : number of z-bins
!   nsteps    : total number of frames in the dump file
!   windowsize: frames per analysis window
!   dt_output : time between dumped frames in fs (dump_every * timestep)
!
! OUTPUT:
!   ASCII file: time(fs)  MSD_bin1  MSD_bin2 ... MSD_binN   [units: fs, Å²]

program msdcom
IMPLICIT NONE

double precision, allocatable :: coord(:,:)   !current positions
double precision, allocatable :: coord0(:,:)  !reference positions for current window
double precision, allocatable :: msd(:,:)     !accumulated MSD
integer,          allocatable :: counter(:,:) !atom count per bin per lag

double precision :: zlo, zhi      ! z-extent of the box
double precision :: binsize       ! z width of each bin
double precision :: dist(3)       ! displacement vector for one atom
double precision :: lx, ly        ! box lengths in x and y (read each frame for PBC)
double precision :: dt_output     ! time between output frames in fs

integer :: i, j, k, l, step_count
integer :: nsteps, natoms, nbins, windowsize, nwindows

character(LEN=50)  :: inputfile, inputtraj, outputfile
character(LEN=200) :: dummy
integer            :: dummyint
integer            :: itype, imol

! reads input file
CALL getarg(1, inputfile)
OPEN(9, file=inputfile)
READ(9,*) inputtraj    ! path to LAMMPS dump file
READ(9,*) outputfile   ! path for MSD output
READ(9,*) nbins        ! number of z-bins
READ(9,*) nsteps       ! total frames in trajectory
READ(9,*) windowsize   ! frames per sliding window
READ(9,*) dt_output    ! fs between frames (dump_freq * timestep_fs)
CLOSE(9)

OPEN(10, file=inputtraj)
zhi = -1.0d10
zlo =  1.0d10

! Read first frame to get natoms and z-range
READ(10,*) dummy          ! ITEM: TIMESTEP
READ(10,*) dummy          ! timestep value
READ(10,*) dummy          ! ITEM: NUMBER OF ATOMS
READ(10,*) natoms         ! atom count

!allocate all arrays now that natoms is known
allocate(coord(3,natoms), coord0(3,natoms), &
         msd(nbins,windowsize), counter(nbins,windowsize))

READ(10,*) dummy          ! ITEM: BOX BOUNDS
READ(10,*) dummy          ! xlo xhi  
READ(10,*) dummy          ! ylo yhi
READ(10,*) dummy          ! zlo zhi  
READ(10,*) dummy          ! ITEM: ATOMS

do i = 1, natoms
  !id  type  mol  xu  yu  zu
  READ(10,*) dummyint, itype, imol, coord0(1,i), coord0(2,i), coord0(3,i)
  if (coord0(3,i) > zhi) zhi = coord0(3,i)
  if (coord0(3,i) < zlo) zlo = coord0(3,i)
enddo

! Rewind the file to start from the beginning
REWIND(10)

!processing
nwindows = nsteps / windowsize
binsize  = (zhi - zlo) / nbins

write(*,'(A,F10.4,A,F10.4,A,I6)') &
  ' z-range: ', (zhi-zlo), ' Å  |  bin size: ', binsize, ' Å  |  windows: ', nwindows

msd     = 0.0d0
counter = 0

! Main processing loop
step_count = 0

do l = 1, nwindows          ! loop over time windows
  
  do j = 1, windowsize      ! loop over lag times within the window
    step_count = step_count + 1
    
    ! --- Read frame header ---
    READ(10,*) dummy          ! ITEM: TIMESTEP
    READ(10,*) dummy          ! timestep value
    READ(10,*) dummy          ! ITEM: NUMBER OF ATOMS
    READ(10,*) dummy          ! atom count (should match natoms)
    READ(10,*) dummy          ! ITEM: BOX BOUNDS
    READ(10,*) dummy          ! xlo xhi  
    READ(10,*) dummy          ! ylo yhi
    READ(10,*) dummy          ! zlo zhi  
    READ(10,*) dummy          ! ITEM: ATOMS

    write(*,'(A,I5,A,I5,A,I6)') ' Window ', l, '  Step ', j, '  Frame ', step_count

    if (j == 1) then
      !first frame of window: store as reference positions
      do i = 1, natoms
        READ(10,*) dummyint, itype, imol, coord0(1,i), coord0(2,i), coord0(3,i)
      enddo
    else
      ! subsequent frames: read current positions and compute MSD
      do i = 1, natoms
        READ(10,*) dummyint, itype, imol, coord(1,i), coord(2,i), coord(3,i)
      enddo
      
      !Accumulate MSD for each atom
      do i = 1, natoms
        ! Find the bin this atom belongs to
        do k = 1, nbins
          if (coord0(3,i) >= (zlo + binsize*(k-1)) .and. &
              coord0(3,i) <  (zlo + binsize*k)) then

            dist(1) = coord(1,i) - coord0(1,i)
            dist(2) = coord(2,i) - coord0(2,i)

            msd(k,j)     = msd(k,j) + dist(1)**2 + dist(2)**2
            counter(k,j) = counter(k,j) + 1
          endif
        enddo
      enddo
    endif

  enddo

enddo 
CLOSE(10)

! Normalize and write output
OPEN(20, file=outputfile)

do j = 1, windowsize
  do k = 1, nbins
    ! if a bin had fewer than 10 atom-frames at any lag,
    ! zero it out for all lags (avoids noisy sparse-bin data)
    if (counter(k,j) < 10) then
      msd(k,:) = 0.0d0
    endif
    ! Normalize by number of contributing atoms
    if (counter(k,j) /= 0) msd(k,j) = msd(k,j) / counter(k,j)
  enddo
  ! Write: time in fs (lag index * time between frames), then MSD per bin
  write(20,'(101F16.8)') (j * dt_output), msd(:,j)
enddo

CLOSE(20)

deallocate(coord, coord0, msd, counter)

end program msdcom
