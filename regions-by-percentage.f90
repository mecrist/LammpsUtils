! ===========================================================================
! regions-by-percentage.f90 — Split trajectory by molecular residence time
! ===========================================================================
!
! PURPOSE:
!   Divides the simulation box along z into nbins equal slabs and tracks
!   each molecule's centre-of-mass (COM) z-position across all frames.
!   A molecule is assigned to a slab if its COM spends more than
!   'timeinregionthr' fraction of the trajectory there.
!   The trajectory is then re-read and atoms belonging to each slab are
!   written to separate LAMMPS dump files.
!
!   useful for separating p.e. bulk water / interfacial water / adsorbed
!   molecules, or analysing populations of ions at different depths.
!
! DUMP FORMAT EXPECTED: unifies dump
!   dump  d1 all custom 1000 traj.lammpstrj id type mol xu yu zu
!   dump_modify d1 sort id
!   Columns: id  type  mol  xu  yu  zu
!
! *** THINGS YOU MUST EDIT (marked USER EDIT) ***
!   nstep, nbins, timeinregionthr

program regions_by_percentage

IMPLICIT NONE

double precision, allocatable :: coord(:,:)       ! atom coords (natm, 3)
double precision, allocatable :: molcm(:,:)       ! molecule COM (maxmol, 3)
double precision, allocatable :: timeinregion(:,:)! fraction of time mol spends in each bin (maxmol, nbins)
integer,          allocatable :: idx(:), mol(:), itype(:)
integer,          allocatable :: molatmnum(:)     ! number of atoms per molecule
integer,          allocatable :: mark(:,:)        ! which bin mol is in at each step (maxmol, nstep)
integer,          allocatable :: numatmbin(:)     ! total atoms in each bin (for output header)
integer,          allocatable :: markporenumber(:)

double precision, dimension(3,2) :: box   ! box(dim, lo/hi)
double precision :: binsize               ! z width per bin

character(LEN=30)  :: dummy1, dm(30), inputname
character(LEN=50),  allocatable :: filename(:)
integer :: filenum

integer :: i, j, k, l, m, natm, maxmol

! ===========================================================================
! USER EDIT — Analysis parameters===========================================
! ===========================================================================
integer,          parameter :: nstep           = 4000  ! total frames in trajectory
integer,          parameter :: nbins           = 80    ! number of z-slabs
double precision, parameter :: timeinregionthr = 0.5d0 ! residence fraction threshold (0–1)
! ===========================================================================

!dead dump file path from command line
CALL getarg(1, inputname)

open(20, file=inputname)

read(20,*) dummy1   ! ITEM: TIMESTEP
read(20,*) dummy1   ! timestep value
read(20,*) dummy1   ! ITEM: NUMBER OF ATOMS
read(20,*) natm
read(20,*) dummy1   ! ITEM: BOX BOUNDS
read(20,*) box(1,1), box(1,2)
read(20,*) box(2,1), box(2,2)
read(20,*) box(3,1), box(3,2)
read(20,*) dummy1   ! ITEM: ATOMS

allocate(coord(natm,3), itype(natm), mol(natm))
allocate(numatmbin(nbins), filename(nbins+20))
allocate(markporenumber(nbins))

!read atom data of frame 1 to find maxmol (highest molecule id)
maxmol = 0
do i = 1, natm
  ! Unified columns: id  type  mol  xu  yu  zu
  ! we only need type, mol here; coords are dummy-read
  read(20,*) dm(1), itype(i), mol(i), dm(2), dm(3), dm(4)
  if (mol(i) > maxmol) maxmol = mol(i)
enddo

write(*,'(I8,A)') natm,   ' atoms'
write(*,'(I8,A)') maxmol, ' molecules'

allocate(molatmnum(maxmol), molcm(maxmol,3))
allocate(timeinregion(maxmol,nbins))
allocate(mark(maxmol,nstep))

rewind(20)

binsize = (box(3,2) - box(3,1)) / nbins

! read all frames, track which z-bin each molecule's COM is in
mark = 1   ! default: bin 1

do j = 1, nstep

  molatmnum = 0
  molcm     = 0.0d0

  !frame header
  read(20,*) dummy1
  read(20,*) dummy1
  read(20,*) dm(1)
  read(20,*) natm
  read(20,*) dummy1
  read(20,*) box(1,1), box(1,2)
  read(20,*) box(2,1), box(2,2)
  read(20,*) box(3,1), box(3,2)
  read(20,*) dm(1)

  !atom positions; accumulate into molecule COM
  do i = 1, natm
    ! Unified columns: id  type  mol  xu  yu  zu
    read(20,*) dm(1), dm(2), mol(i), coord(i,1), coord(i,2), coord(i,3)
  enddo

  ! Compute COM for each molecule (OpenMP parallelised)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,i,k) &
  !$OMP& REDUCTION(+:molcm, molatmnum)
  do l = 1, maxmol
    do i = 1, natm
      if (mol(i) == l) then
        molcm(l,1) = molcm(l,1) + coord(i,1)
        molcm(l,2) = molcm(l,2) + coord(i,2)
        molcm(l,3) = molcm(l,3) + coord(i,3)
        molatmnum(l) = molatmnum(l) + 1
      endif
    enddo
    !finalise COM (safe outside REDUCTION since each l is independent)
    if (molatmnum(l) > 0) then
      molcm(l,1) = molcm(l,1) / molatmnum(l)
      molcm(l,2) = molcm(l,2) / molatmnum(l)
      molcm(l,3) = molcm(l,3) / molatmnum(l)
    endif
    ! Assign molecule to its z-bin based on COM
    do k = 1, nbins
      if (molcm(l,3) >= (k-1)*binsize .and. molcm(l,3) < k*binsize) then
        mark(l,j) = k
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO

  write(*,'(A,I6,A,I6)') 'Pass 1 frame ', j, ' / ', nstep

enddo

! Computes fraction of time each molecule spends in each bin
timeinregion = 0.0d0

do j = 1, nstep
  do l = 1, maxmol
    do k = 1, nbins
      if (mark(l,j) == k) timeinregion(l,k) = timeinregion(l,k) + 1.0d0
    enddo
  enddo
enddo

!normalising to fraction (0–1)
timeinregion = timeinregion / nstep

! atoms per bin (for the output dump)
numatmbin = 0
do l = 1, maxmol
  do k = 1, nbins
    if (timeinregion(l,k) > timeinregionthr) then
      numatmbin(k) = numatmbin(k) + molatmnum(l)
    endif
  enddo
enddo

! open one output file per populated bin
do k = 1, nbins
  filenum = 20 + k
  write(filename(filenum),'(A28,I4.4,A10)') &
    './regions-by-percentage/bin-', k, '.lammpstrj'
  if (numatmbin(k) > 0) open(filenum, file=filename(filenum))
enddo

! Pass 2: re-read trajectory, write selected atoms to per-bin dump files
rewind(20)

do j = 1, nstep

  write(*,'(A,I6,A,I6)') 'Writing frame ', j, ' / ', nstep

  ! Read the frame header as raw tokens (to re-write them)
  read(20,*) dm(1), dm(2)                              ! ITEM: TIMESTEP
  read(20,*) dm(3)                                      ! timestep value
  read(20,*) dm(4), dm(5), dm(6), dm(7)                ! ITEM: NUMBER OF ATOMS
  read(20,*) natm
  read(20,*) dm(9), dm(10), dm(11), dm(12), dm(13), dm(14)  ! ITEM: BOX BOUNDS
  read(20,*) box(1,1), box(1,2)
  read(20,*) box(2,1), box(2,2)
  read(20,*) box(3,1), box(3,2)
  read(20,*) dm(15), dm(16), dm(17), dm(18), dm(19), dm(20), dm(21), dm(22)  ! ITEM: ATOMS ...

  ! Write header to each active bin file
  do k = 1, nbins
    if (numatmbin(k) > 0) then
      filenum = 20 + k
      write(filenum,'(A,1X,A)')   dm(1), dm(2)
      write(filenum,'(A)')        dm(3)
      write(filenum,'(4A)')       dm(4), ' ', dm(5), ' '
      write(filenum,'(I6)')       numatmbin(k)      ! atom count for this bin
      write(filenum,'(6A)')       dm(9), ' ', dm(10), ' ', dm(11), ' '
      write(filenum,'(2F18.8)')   box(1,1), box(1,2)
      write(filenum,'(2F18.8)')   box(2,1), box(2,2)
      write(filenum,'(2F18.8)')   box(3,1), box(3,2)
      write(filenum,'(8A)')       dm(15),' ',dm(16),' ',dm(17),' ',dm(18),' '
    endif
  enddo

  ! read and distribute atoms
  do i = 1, natm
    ! Unified columns: id  type  mol  xu  yu  zu
    read(20,*) dm(1), dm(2), mol(i), coord(i,1), coord(i,2), coord(i,3)
    do k = 1, nbins
      if (numatmbin(k) > 0) then
        if (timeinregion(mol(i),k) > timeinregionthr) then
          filenum = 20 + k
          write(filenum,'(A,1X,A,I8,3F18.8)') &
            trim(dm(1)), trim(dm(2)), mol(i), coord(i,1), coord(i,2), coord(i,3)
        endif
      endif
    enddo
  enddo

enddo

close(20)
do k = 1, nbins
  if (numatmbin(k) > 0) close(20+k)
enddo

deallocate(coord, molcm, timeinregion, idx, mol, itype, molatmnum, mark, numatmbin)

end program regions_by_percentage
