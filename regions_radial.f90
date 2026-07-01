! ===========================================================================
! regions-by-radius.f90 — Split trajectory into inner/outer radial regions
! ===========================================================================
!
! PURPOSE:
!   Divides the simulation box into two regions based on radial distance
!   from the pore axis (assumed to be along z, centred in XY).
!   Each molecule's centre of mass (COM) is tracked across all frames.
!   Molecules are assigned to the inner region if their COM is within
!   R_inner of the axis, otherwise to the outer region.
!   Only molecules that spend more than 'timeinregionthr' fraction of the
!   trajectory in a given region are written to that region's dump file.
!
! DUMP FORMAT EXPECTED: unified dump
!   dump d1 all custom 1000 traj.lammpstrj id type mol xu yu zu
!   dump_modify d1 sort id
!   Columns: id  type  mol  xu  yu  zu
!
! *** THINGS YOU MUST EDIT (marked USER EDIT) ***
!   nstep, R_inner, timeinregionthr

program regions_by_radius

IMPLICIT NONE

double precision, allocatable :: coord(:,:)       ! atom coords (natm, 3)
double precision, allocatable :: molcm(:,:)       ! molecule COM (maxmol, 3)
double precision, allocatable :: timeinregion(:,:)! fraction of time mol spends in each bin (maxmol, 2)
integer,          allocatable :: idx(:), mol(:), itype(:)
integer,          allocatable :: molatmnum(:)     ! number of atoms per molecule
integer,          allocatable :: mark(:,:)        ! which region mol is in at each step (maxmol, nstep)
integer,          allocatable :: numatmbin(:)     ! total atoms in each region (for output header)

double precision, dimension(3,2) :: box   ! box(dim, lo/hi)
double precision :: cx, cy                   ! centre of the pore in XY (box centre)
double precision :: r2, R_inner2            ! squared radii for comparison

character(LEN=30)  :: dummy1, dm(30), inputname
character(LEN=50),  allocatable :: filename(:)
integer :: filenum

integer :: i, j, k, l, m, natm, maxmol

! ===========================================================================
! USER EDIT — Analysis parameters
! ===========================================================================
integer,          parameter :: nstep           = 2000    ! total frames in trajectory
double precision, parameter :: R_inner         = 7.0d0  ! inner radius (Å)
double precision, parameter :: timeinregionthr = 0.9d0   ! residence fraction threshold (0–1)
! ===========================================================================

! dump file path from command line
CALL getarg(1, inputname)

open(20, file="traj.lammpstrj")

read(20,*) dummy1   ! ITEM: TIMESTEP
read(20,*) dummy1   ! timestep value
read(20,*) dummy1   ! ITEM: NUMBER OF ATOMS
read(20,*) natm
read(20,*) dummy1   ! ITEM: BOX BOUNDS
read(20,*) box(1,1), box(1,2)
read(20,*) box(2,1), box(2,2)
read(20,*) box(3,1), box(3,2)
read(20,*) dummy1   ! ITEM: ATOMS

! compute box centre in XY (assumes pore axis at centre)
cx = 0.5d0 * (box(1,1) + box(1,2))
cy = 0.5d0 * (box(2,1) + box(2,2))
R_inner2 = R_inner * R_inner

allocate(coord(natm,3), itype(natm), mol(natm))
allocate(numatmbin(2), filename(2+20))

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
allocate(timeinregion(maxmol,2))
allocate(mark(maxmol,nstep))

rewind(20)

! read all frames, track which region each molecule's COM is in
mark = 1   ! default: inner

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
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,i) &
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
    !finalise COM
    if (molatmnum(l) > 0) then
      molcm(l,1) = molcm(l,1) / molatmnum(l)
      molcm(l,2) = molcm(l,2) / molatmnum(l)
      molcm(l,3) = molcm(l,3) / molatmnum(l)
    endif
    ! Assign molecule to radial region based on COM
    r2 = (molcm(l,1) - cx)**2 + (molcm(l,2) - cy)**2
    if (r2 <= R_inner2) then
      mark(l,j) = 1
    else
      mark(l,j) = 2
    endif
  enddo
  !$OMP END PARALLEL DO

  write(*,'(A,I6,A,I6)') 'Pass 1 frame ', j, ' / ', nstep

enddo

! Compute fraction of time each molecule spends in each region
timeinregion = 0.0d0

do j = 1, nstep
  do l = 1, maxmol
    do k = 1, 2
      if (mark(l,j) == k) timeinregion(l,k) = timeinregion(l,k) + 1.0d0
    enddo
  enddo
enddo

!normalising to fraction (0–1)
timeinregion = timeinregion / nstep

! atoms per region (for output header)
numatmbin = 0
do l = 1, maxmol
  do k = 1, 2
    if (timeinregion(l,k) > timeinregionthr) then
      numatmbin(k) = numatmbin(k) + molatmnum(l)
    endif
  enddo
enddo

! open one output file per region
filenum = 21
open(filenum, file='inner.lammpstrj')
filenum = 22
open(filenum, file='outer.lammpstrj')

! Pass 2: re-read trajectory, write selected atoms to region dump files
rewind(20)

do j = 1, nstep

  write(*,'(A,I6,A,I6)') 'Writing frame ', j, ' / ', nstep

  ! Read the frame header as raw tokens
  read(20,*) dm(1), dm(2)                              ! ITEM: TIMESTEP
  read(20,*) dm(3)                                      ! timestep value
  read(20,*) dm(4), dm(5), dm(6), dm(7)                ! ITEM: NUMBER OF ATOMS
  read(20,*) natm
  read(20,*) dm(9), dm(10), dm(11), dm(12), dm(13), dm(14)  ! ITEM: BOX BOUNDS
  read(20,*) box(1,1), box(1,2)
  read(20,*) box(2,1), box(2,2)
  read(20,*) box(3,1), box(3,2)
  read(20,*) dm(15), dm(16), dm(17), dm(18), dm(19), dm(20), dm(21), dm(22)  ! ITEM: ATOMS ...

  ! Write header to each region file
  do k = 1, 2
    filenum = 20 + k
    if (numatmbin(k) > 0) then
      write(filenum,'(A,1X,A)')   dm(1), dm(2)
      write(filenum,'(A)')        dm(3)
      write(filenum,'(4A)')       dm(4), ' ', dm(5), ' '
      write(filenum,'(I6)')       numatmbin(k)      ! atom count for this region
      write(filenum,'(6A)')       dm(9), ' ', dm(10), ' ', dm(11), ' '
      write(filenum,'(2F18.8)')   box(1,1), box(1,2)
      write(filenum,'(2F18.8)')   box(2,1), box(2,2)
      write(filenum,'(2F18.8)')   box(3,1), box(3,2)
      write(filenum,'(8A)')       dm(15),' ',dm(16),' ',dm(17),' ',dm(18),' '
    endif
  enddo

  ! read and distribute atoms
  do i = 1, natm
    read(20,*) dm(1), dm(2), mol(i), coord(i,1), coord(i,2), coord(i,3)
    do k = 1, 2
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
do k = 1, 2
  if (numatmbin(k) > 0) close(20+k)
enddo

deallocate(coord, molcm, timeinregion, idx, mol, itype, molatmnum, mark, numatmbin)

end program regions_by_radius
