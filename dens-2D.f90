! ===========================================================================
! dens-2D.f90 — mass density profile along z, resolved by molecule type
! ===========================================================================
!
!   reads a LAMMPS trajectory and computes the time-averaged mass density
!   profile ρ(z) [g/cm³] for each molecule type separately.
!   useful for visualising how different species distribute at interfaces.
!
! DUMP FORMAT EXPECTED (unified dump, shared with all other tools, except S_tensions):
!   dump  d1 all custom 1000 traj.lammpstrj id type mol xu yu zu
!   dump_modify d1 sort id
!   Columns: id  type  mol  xu  yu  zu
!   Note: xu/yu/zu (unwrapped) are used here only for z-binning;
!         wrapped z would give identical results for NVT (box fixed).
!
! *** INPUTS (marked with USER EDIT) ***
!   1. nconf, nbins         — frame count and z-resolution
!   2. nspecies, nmoltype   — number of LAMMPS atom types and molecule classes
!   3. atmmass(1..nspecies) — atomic mass for each LAMMPS type [g/mol]
!   4. moltype assignments  — which atom type identifies each molecule class
!   5. inputtraj, outputfile — file paths
! ===========================================================================

program dens2d
implicit none

integer :: ii, jj, kk, ll

integer :: nnumber, ntype, nmol   ! atom id, LAMMPS type, molecule id
integer :: natoms, dummyint

double precision :: lox, hix, loy, hiy, lo, hi, lo1, hi1
double precision :: dr            ! bin width in z (Å)
double precision :: low, high     ! bin lower/upper bounds
double precision :: d, r          ! z-coordinate of atom; bin centre
double precision :: denominador   ! bin volume for density normalisation

double precision, dimension(:), allocatable :: xx, yy, zz
integer,          dimension(:), allocatable :: nntype, nnnumber, nnmol, moltype

! --- density arrays ---
!   num(moltype, species)        — atom count per type per species in current bin
!   nr(moltype, species, bin)    — same, saved per bin
!   massnr(moltype, bin)         — total mass per moltype per bin
!   dens(frame, moltype, bin)    — density per frame
!   densm(moltype)               — time-averaged density (final output)
integer,          allocatable :: num(:,:), nr(:,:,:)
double precision, allocatable :: massnr(:,:), dens(:,:,:), densm(:)
double precision, allocatable :: atmmass(:)

! --- other variables
integer          :: nmol_dummy
integer          :: atomsinside
character(50)    :: dummychar

! ===========================================================================
! USER EDIT — System parameters
! ===========================================================================
integer, parameter :: nbins    = 218  ! number of z-slabs; adjust to box_z / desired_resolution
integer, parameter :: nconf    = 800  ! total frames in trajectory
integer, parameter :: nspecies = 38   ! number of distinct LAMMPS atom types in your system
integer, parameter :: nmoltype = 11   ! number of molecule classes you want to track

character(len=100), parameter :: inputtraj  = 'traj.lammpstrj'  ! unified dump file
character(len=100), parameter :: outputfile = 'dens-2D.dat'
! ===========================================================================

allocate(atmmass(nspecies))
allocate(nntype(1), nnnumber(1), nnmol(1))   ! temporary; reallocated after natoms known

!Read first frame header to get natoms, then allocate position arrays
open(10, file=inputtraj, status='old')

read(10,*) dummychar          ! ITEM: TIMESTEP
read(10,*) dummyint           ! timestep value
read(10,*) dummychar          ! ITEM: NUMBER OF ATOMS
read(10,*) natoms             ! atom count
read(10,*) dummychar          ! ITEM: BOX BOUNDS
read(10,*) lo, hi             ! xlo xhi
read(10,*) lo, hi             ! ylo yhi
read(10,*) lo1, hi1           ! zlo zhi  — used to set bin range
read(10,*) dummychar          ! ITEM: ATOMS

!skip the rest of frame 1
do ii = 1, natoms
  ! id  type  mol  xu  yu  zu
  read(10,*) nnumber, ntype, nmol_dummy, xx(1), xx(1), xx(1)
enddo

deallocate(nntype, nnnumber, nnmol)
allocate(xx(natoms), yy(natoms), zz(natoms))
allocate(nntype(natoms), nnnumber(natoms), nnmol(natoms), moltype(natoms))
allocate(num(nmoltype, nspecies), nr(nmoltype, nspecies, nbins))
allocate(massnr(nmoltype, nbins), dens(nconf, nmoltype, nbins), densm(nmoltype))

dr = (hi1 - lo1) / nbins   ! z-bin width in Å

! ===========================================================================
! USER EDIT — Atomic masses
! one entry per LAMMPS atom type (in order 1..nspecies).
! formula: (1E24/6.02E23) * atomic_mass_g_per_mol
! converts to g (the density normalisation uses Å³, giving g/cm³).
! ===========================================================================
atmmass(1)  = (1E24/6.02E23) * 12.011   
atmmass(2)  = (1E24/6.02E23) * 15.9994  
atmmass(3)  = (1E24/6.02E23) * 14.007   
atmmass(4)  = (1E24/6.02E23) * 12.011
atmmass(5)  = (1E24/6.02E23) * 1.0079
atmmass(6)  = (1E24/6.02E23) * 12.011
atmmass(7)  = (1E24/6.02E23) * 1.0079
atmmass(8)  = (1E24/6.02E23) * 12.011
atmmass(9)  = (1E24/6.02E23) * 1.0079
atmmass(10) = (1E24/6.02E23) * 12.011
atmmass(11) = (1E24/6.02E23) * 1.0079
atmmass(12) = (1E24/6.02E23) * 12.011
atmmass(13) = (1E24/6.02E23) * 1.0079
atmmass(14) = (1E24/6.02E23) * 15.9994  
atmmass(15) = (1E24/6.02E23) * 1.0079   
atmmass(16) = (1E24/6.02E23) * 35.453   
atmmass(17) = (1E24/6.02E23) * 22.990   
atmmass(18) = (1E24/6.02E23) * 32.065   
atmmass(19) = (1E24/6.02E23) * 15.9994
atmmass(20) = (1E24/6.02E23) * 24.305   
atmmass(21) = (1E24/6.02E23) * 40.078   
atmmass(22) = (1E24/6.02E23) * 39.098   
atmmass(23) = (1E24/6.02E23) * 12.011
atmmass(24) = (1E24/6.02E23) * 32.065
atmmass(25) = (1E24/6.02E23) * 15.9994
atmmass(26) = (1E24/6.02E23) * 1.0079
atmmass(27) = (1E24/6.02E23) * 15.9994
atmmass(28) = (1E24/6.02E23) * 1.0079
atmmass(29) = (1E24/6.02E23) * 15.9994
atmmass(30) = (1E24/6.02E23) * 15.9994
atmmass(31) = (1E24/6.02E23) * 12.011
atmmass(32) = (1E24/6.02E23) * 12.011
atmmass(33) = (1E24/6.02E23) * 15.9994
atmmass(34) = (1E24/6.02E23) * 15.9994
atmmass(35) = (1E24/6.02E23) * 1.0079
atmmass(36) = (1E24/6.02E23) * 32.065
atmmass(37) = (1E24/6.02E23) * 12.011
atmmass(38) = (1E24/6.02E23) * 1.0079
! ===========================================================================

rewind(10)

! read all frames and accumulate density
do ll = 1, nconf

  read(10,*) dummychar           ! ITEM: TIMESTEP
  read(10,*) dummyint            ! timestep value
  read(10,*) dummychar           ! ITEM: NUMBER OF ATOMS
  read(10,*) natoms
  read(10,*) dummychar           ! ITEM: BOX BOUNDS
  read(10,*) lox, hix            ! x bounds (used for bin volume)
  read(10,*) loy, hiy            ! y bounds
  read(10,*) lo, hi              ! z bounds
  read(10,*) dummychar           ! ITEM: ATOMS

  do ii = 1, natoms
    read(10,*) nnnumber(ii), nntype(ii), nnmol(ii), xx(ii), yy(ii), zz(ii)

    ! ===========================================================================
    ! USER EDIT — Molecule type assignment
    ! map the LAMMPS atom type to a molecule class (1..nmoltype).
    ! Use the atom type that uniquely identifies each molecule (e.g. the oxygen
    ! of water, the central ion, a unique heavy atom).

    if (nntype(ii) == 14) moltype(nnmol(ii)) = 1   
    if (nntype(ii) == 34) moltype(nnmol(ii)) = 2   
    if (nntype(ii) == 31) moltype(nnmol(ii)) = 3   
    if (nntype(ii) == 32) moltype(nnmol(ii)) = 4   
    if (nntype(ii) == 16) moltype(nnmol(ii)) = 5   
    if (nntype(ii) == 17) moltype(nnmol(ii)) = 6   
    if (nntype(ii) == 18) moltype(nnmol(ii)) = 7   
    if (nntype(ii) == 19) moltype(nnmol(ii)) = 8   
    if (nntype(ii) == 20) moltype(nnmol(ii)) = 9   
    if (nntype(ii) ==  8) moltype(nnmol(ii)) = 10  
    if (nntype(ii) == 22) moltype(nnmol(ii)) = 11  
  enddo

  massnr = 0.0d0

  do jj = 1, nbins
    num  = 0
    low  = lo + (jj-1) * dr    ! bin lower bound in z
    high = lo +  jj    * dr    ! bin upper bound in z

    do ii = 1, natoms
      d = zz(ii)
      if (d >= low .and. d <= high) then
        ! Count atom under molecule type and species
        do kk = 1, nmoltype
          if (moltype(nnmol(ii)) == kk) then
            if (nntype(ii) >= 1 .and. nntype(ii) <= nspecies) then
              num(kk, nntype(ii)) = num(kk, nntype(ii)) + 1
            endif
          endif
        enddo
      endif
    enddo

    !save counts
    do kk = 1, nmoltype
      nr(kk, :, jj) = num(kk, :)
    enddo
  enddo

  ! conerting atom counts to mass density ---
  do jj = 1, nbins
    denominador = dr * (hix - lox) * (hiy - loy)   ! bin volume in Å³
    do kk = 1, nmoltype
      massnr(kk, jj) = 0.0d0
      do ii = 1, nspecies
        massnr(kk, jj) = massnr(kk, jj) + nr(kk, ii, jj) * atmmass(ii)
      enddo
      dens(ll, kk, jj) = massnr(kk, jj) / denominador   ! g/cm³
    enddo
  enddo

  write(*,*) 'Frame ', ll, ' / ', nconf

enddo

close(10)

! Time-average and write output
open(13, file=outputfile)

do jj = 1, nbins
  densm = 0.0d0
  r = lo1 + (jj - 1) * dr 

  do ll = 1, nconf
    do kk = 1, nmoltype
      densm(kk) = densm(kk) + dens(ll, kk, jj)
    enddo
  enddo

  densm = densm / nconf

  write(13, "(20F12.3)") r, densm(:)
enddo

close(13)

deallocate(xx, yy, zz, nntype, nnnumber, nnmol, moltype)
deallocate(num, nr, massnr, dens, densm, atmmass)

end program dens2d
