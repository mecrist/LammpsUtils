! ===========================================================================
! dens-radial.f90 — radial mass density profile ρ(r) [g/cm³] by molecule type
!   Bins atoms by distance r = sqrt(x²+y²) from the pore axis (z-axis),
!   restricted to a z-range that excludes entrance effects.
! ===========================================================================

program dens_radial
  implicit none

  integer :: ii, jj, kk, ll

  integer :: natoms, dummyint
  integer :: ntype_i, nmol_i, nnum_i   ! per-atom scalars for trajectory read

  double precision :: lox, hix, loy, hiy, lo, hi, lo1, hi1
  double precision :: dr                   ! bin width in r (Å)
  double precision :: low, high            ! bin lower/upper r bounds
  double precision :: d, r_center          ! radial coordinate; bin centre
  double precision :: rx, ry, rr           ! x, y, radial distance
  double precision :: denominador          ! bin volume for density normalisation
  double precision :: z_lo_cut, z_hi_cut   ! z-range after trimming entrance effects
  double precision :: L_z                  ! effective z length used in normalisation

  double precision, dimension(:), allocatable :: xx, yy, zz
  integer,          dimension(:), allocatable :: nntype, nnnumber, nnmol, moltype

  ! massnr(kk, jj)  — mass per molecule class per bin, current frame
  ! dens(ll, kk, jj) — density per frame per class per bin (for time-average)
  ! densm(kk)        — time-averaged density per class for output
  double precision, allocatable :: massnr(:,:), dens(:,:,:), densm(:)
  double precision, allocatable :: atmmass(:)

  integer          :: bin_idx              ! direct bin index for an atom
  character(50)    :: dummychar

  ! ===========================================================================
  ! USER EDIT — System parameters
  ! ===========================================================================
  integer, parameter :: nbins    = 42
  integer, parameter :: nconf    = 2000
  integer, parameter :: nspecies = 8      ! number of distinct LAMMPS atom types
  integer, parameter :: nmoltype = 2      ! number of molecule classes

  ! Pore axis is z; pore spans z = -46 to +46 Å
  ! Trim 15 Å (1.5 nm) from each end to remove entrance effects
  double precision, parameter :: z_pore_lo =  -46.0d0
  double precision, parameter :: z_pore_hi =   46.0d0
  double precision, parameter :: z_trim    =   15.0d0   ! 1.5 nm in Å

  double precision, parameter :: r_max     =   42.0d0   ! max radial distance (Å)

  character(len=100), parameter :: inputtraj  = 'traj.lammpstrj'
  character(len=100), parameter :: outputfile = 'dens-radial.dat'
  ! ===========================================================================

  allocate(atmmass(nspecies))

  ! -------------------------------------------------------------------
  ! First pass: read one header to get natoms, then allocate arrays
  ! -------------------------------------------------------------------
  open(10, file=inputtraj, status='old')

  read(10,*) dummychar
  read(10,*) dummyint
  read(10,*) dummychar
  read(10,*) natoms
  read(10,*) dummychar
  read(10,*) lo, hi
  read(10,*) lo, hi
  read(10,*) lo1, hi1
  read(10,*) dummychar

  do ii = 1, natoms
    read(10,*) nnum_i, ntype_i, nmol_i, rx, ry, d
  enddo

  allocate(xx(natoms), yy(natoms), zz(natoms))
  allocate(nntype(natoms), nnnumber(natoms), nnmol(natoms), moltype(natoms))
  allocate(massnr(nmoltype, nbins))
  allocate(dens(nconf, nmoltype, nbins))
  allocate(densm(nmoltype))

  ! Radial bin width
  dr = r_max / dble(nbins)

  ! Effective z-length after trimming (used in cylindrical shell volume)
  z_lo_cut = z_pore_lo + z_trim    ! = -31.0 Å
  z_hi_cut = z_pore_hi - z_trim    ! = +31.0 Å
  L_z      = z_hi_cut - z_lo_cut   ! =  62.0 Å

  ! ===========================================================================
  ! USER EDIT — Atomic masses [g] (= amu * 1e24 / N_A)
  ! ===========================================================================
  atmmass(1)  = (1.0d24/6.02d23) * 1.0079d0    ! H  (water)
  atmmass(2)  = (1.0d24/6.02d23) * 15.9994d0   ! O  (water)
  atmmass(3)  = (1.0d24/6.02d23) * 15.9994d0   ! O  (silica)
  atmmass(4)  = (1.0d24/6.02d23) * 28.0855d0   ! Si
  atmmass(5)  = (1.0d24/6.02d23) * 28.0855d0   ! Si
  atmmass(6)  = (1.0d24/6.02d23) * 15.9994d0   ! O
  atmmass(7)  = (1.0d24/6.02d23) * 15.9994d0   ! O  (water, type 7)
  atmmass(8)  = (1.0d24/6.02d23) * 1.0079d0    ! H  (water, type 8)
  ! ===========================================================================

  rewind(10)

  ! ===================================================================
  ! Main loop over frames
  ! ===================================================================
  do ll = 1, nconf

    read(10,*) dummychar
    read(10,*) dummyint
    read(10,*) dummychar
    read(10,*) natoms
    read(10,*) dummychar
    read(10,*) lox, hix
    read(10,*) loy, hiy
    read(10,*) lo, hi
    read(10,*) dummychar

    ! Read atom data (id type mol xu yu zu)
    do ii = 1, natoms
      read(10,*) nnnumber(ii), nntype(ii), nnmol(ii), xx(ii), yy(ii), zz(ii)
      ! Assign molecule class by atom type
      if (nntype(ii) == 7 .or. nntype(ii) == 8) then
        moltype(ii) = 1   ! Water
      else
        moltype(ii) = 2   ! Silica
      endif
    enddo

    massnr = 0.0d0

    ! -------------------------------------------------------------------
    ! O(natoms) binning: compute bin index directly, no loop over bins
    ! -------------------------------------------------------------------
    do ii = 1, natoms

      ! z filter: only atoms in the trimmed pore region
      if (zz(ii) < z_lo_cut .or. zz(ii) > z_hi_cut) cycle

      ! Radial distance from pore axis
      rr = sqrt(xx(ii)**2 + yy(ii)**2)
      if (rr >= r_max) cycle

      ! Direct bin index (1-based)
      bin_idx = int(rr / dr) + 1
      if (bin_idx > nbins) cycle   ! guard against rr exactly == r_max

      kk = moltype(ii)   ! FIX: was moltype(nnmol(ii)) — nnmol was unpopulated

      if (nntype(ii) >= 1 .and. nntype(ii) <= nspecies) then
        massnr(kk, bin_idx) = massnr(kk, bin_idx) + atmmass(nntype(ii))
      endif

    enddo

    ! Convert mass to density: divide by cylindrical shell volume
    do jj = 1, nbins
      low  = (jj-1) * dr
      high =  jj    * dr
      denominador = 3.141592653589793d0 * (high**2 - low**2) * L_z

      do kk = 1, nmoltype
        dens(ll, kk, jj) = massnr(kk, jj) / denominador
      enddo
    enddo

    write(*,*) 'Frame ', ll, ' / ', nconf

  enddo

  close(10)

  ! ===================================================================
  ! Time-average and write output
  ! ===================================================================
  open(13, file=outputfile)
  write(13, "(A)") "# r(Angstrom)   dens_water   dens_silica   [g/cm^3]"
  write(13, "(A)") "# z trimmed to: z_lo_cut=-31 Ang, z_hi_cut=+31 Ang"

  do jj = 1, nbins
    densm    = 0.0d0
    r_center = (jj - 0.5d0) * dr

    do ll = 1, nconf
      do kk = 1, nmoltype
        densm(kk) = densm(kk) + dens(ll, kk, jj)
      enddo
    enddo

    densm = densm / dble(nconf)

    write(13, "(F10.4, 20F12.5)") r_center, densm(:)
  enddo

  close(13)

  deallocate(xx, yy, zz, nntype, nnnumber, nnmol, moltype)
  deallocate(massnr, dens, densm, atmmass)

end program dens_radial
