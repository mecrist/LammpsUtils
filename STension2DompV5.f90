! ===========================================================================
! STension2DompV5.f90 — interfacial/surface tension from per-atom stress
! ===========================================================================
!
! computes the interfacial (surface) tension γ from the per-atom pressure
! tensor using the mechanical (Irving-Kirkwood) definition:
!
!     γ = (1/2) × ∫ [Pzz(z) - (Pxx(z) + Pyy(z))/2] dz
!
! the factor 1/2 accounts for two interfaces (top and bottom).
! also outputs the time-averaged local pressure anisotropy profile P(z).
!
! units:
!   lammps stress/atom outputs in pressure × volume [bar·å³] (real/metal units).
!   dividing by bin volume gives [bar].
!   integration over z gives [bar·å].
!   conversion: bar·å × 1e5 (pa/bar) × 1e-10 (m/å) × 1e3 (mn/m per n/m) = × 1e-2
!   final unit: mn/m  (= dyn/cm in cgs, the common experimental unit for γ).
!
! sign convention:
!   lammps stress/atom stores the negative of the pressure tensor (compressive
!   positive). we apply a global sign flip when binning so that the result
!   represents the standard mechanical pressure (tensile positive).
!
! dump format expected:
!   compute   stress all stress/atom NULL virial
!   dump      d2 all custom 5000 traj_stress.lammpstrj &
!             id type mol xu yu zu c_stress[1] c_stress[2] c_stress[3]
!   dump_modify d2 sort id
!   columns: id  type  mol  xu  yu  zu  pxx  pyy  pzz
!
!
! THINGS YOU MUST EDIT! (marked user edit):
!   nconf, nbins, jump, inputtraj, outputfile
!
! output:
!   screen (each frame): frame  γ_this_frame(mn/m)  γ_running_avg(mn/m)
!   screen (final):      γ_avg  std_dev  std_err   (all in mn/m)
!   file pressures.dat:  z(å)  p(z)(bar)           (time-averaged profile)
! ===========================================================================
program stension_final
  implicit none

  ! ---------------------------------------------------------------------------
  ! user edit — analysis parameters
  ! ---------------------------------------------------------------------------
  integer,            parameter :: nconf      = 800    ! production frames to process
  integer,            parameter :: nbins      = 210    ! z-slabs along interface normal
  integer,            parameter :: jump       = 0      ! equilibration frames to skip; 0 if not present
  character(len=100), parameter :: inputtraj  = 'traj_stress.lammpstrj'
  character(len=100), parameter :: outputfile = 'pressures.dat'
  ! ---------------------------------------------------------------------------

  integer :: i, j, b, l, k

  ! per-atom arrays
  integer,          allocatable :: nnumber(:), nntype(:)
  real*8,           allocatable :: zz(:)
  real*8,           allocatable :: pxx_a(:), pyy_a(:), pzz_a(:)
  real*8,           allocatable :: ptot(:)       ! local pressure anisotropy per atom

  ! box bounds (updated each frame)
  real*8 :: lox, hix, loy, hiy, loz, hiz
  real*8 :: boxx, boxy, boxz, dr

  ! binned pressure and surface tension
  real*8, allocatable :: soma(:,:)    ! pressure anisotropy per bin, per frame (nbins, nconf)
  real*8, allocatable :: temp1(:)     ! time-averaged pressure per bin
  real*8, allocatable :: sigma(:)     ! γ per frame (mn/m)
  real*8, allocatable :: dpdr(:)      ! integral of p(z) per frame (bar·å)

  real*8 :: sigmabig                  ! running sum of γ
  real*8 :: gamma_mean, std, ste      ! final statistics

  integer :: dummyint, natoms
  real*8  :: xu_dummy, yu_dummy, zu_dummy
  integer :: ntype_dummy, nmol_dummy
  integer          :: t0_count, t1_count, count_rate_int, count_max_int
  real*8           :: wall_start, wall_end

  ! ---------------------------------------------------------------------------
  ! timing start
  ! ---------------------------------------------------------------------------
  call system_clock(t0_count, count_rate_int, count_max_int)
  wall_start = dble(t0_count) / dble(count_rate_int)

  ! ---------------------------------------------------------------------------
  ! open dump file and skip 'jump' equilibration frames
  ! IF PRODUCTION ONLY, set jump = 0
  ! ---------------------------------------------------------------------------
  open(9, file=inputtraj, status='old')

  do i = 1, jump
    write(*,'(a,i6)') 'skipping frame ', i

    read(9,*) dummyint   ! item: timestep (read header as int to be safe)
    read(9,*) dummyint   ! timestep value
    read(9,*) dummyint   ! item: number of atoms
    read(9,*) natoms
    read(9,*) dummyint   ! item: box bounds
    read(9,*) lox, hix
    read(9,*) loy, hiy
    read(9,*) loz, hiz
    read(9,*) dummyint   ! item: atoms

    ! allocate temporary arrays on first skip frame, reuse for subsequent ones
    if (i == 1) then
      allocate(pxx_a(natoms), pyy_a(natoms), pzz_a(natoms))
    endif

    do j = 1, natoms
      read(9,*) dummyint, ntype_dummy, nmol_dummy, &
                xu_dummy, yu_dummy, zu_dummy, &
                pxx_a(j), pyy_a(j), pzz_a(j)
    enddo
  enddo

  ! free skip arrays before main allocation (if jump > 0)
  if (jump > 0) deallocate(pxx_a, pyy_a, pzz_a)

  ! ---------------------------------------------------------------------------
  ! read first production frame: get natoms and box for bin setup
  ! ---------------------------------------------------------------------------
  read(9,*) dummyint   ! item: timestep
  read(9,*) dummyint   ! timestep value
  read(9,*) dummyint   ! item: number of atoms
  read(9,*) natoms
  read(9,*) dummyint   ! item: box bounds
  read(9,*) lox, hix
  read(9,*) loy, hiy
  read(9,*) loz, hiz
  read(9,*) dummyint   ! item: atoms

  boxx = hix - lox
  boxy = hiy - loy
  boxz = hiz - loz
  dr   = boxz / nbins

  write(*,'(a,f10.4,a)') ' z-bin width: ', dr, ' å'

  ! ---------------------------------------------------------------------------
  ! allocate all working arrays
  ! ---------------------------------------------------------------------------
  allocate(nnumber(natoms), nntype(natoms))
  allocate(zz(natoms))
  allocate(pxx_a(natoms), pyy_a(natoms), pzz_a(natoms))
  allocate(ptot(natoms))
  allocate(soma(nbins, nconf), temp1(nbins))
  allocate(sigma(nconf), dpdr(nconf))

  ! ---------------------------------------------------------------------------
  ! read first production frame atom data
  ! ---------------------------------------------------------------------------
  do j = 1, natoms
    read(9,*) nnumber(j), nntype(j), nmol_dummy, &
              xu_dummy, yu_dummy, zz(j), &
              pxx_a(j), pyy_a(j), pzz_a(j)
    ! shift z so that loz maps to 0 for bin indexing
    zz(j)   = zz(j) - loz
    ! local pressure anisotropy integrand: pzz - (pxx + pyy)/2
    ! units: bar·å³
    ptot(j) = pzz_a(j) - 0.5d0 * (pxx_a(j) + pyy_a(j))
  enddo

  ! bin first frame
  do l = 1, nbins
    soma(l,1) = 0.0d0
    do k = 1, natoms
      ! lower bound inclusive on bin 1; upper bound inclusive on last bin
      if (zz(k) >= (l-1)*dr .and. zz(k) < l*dr) then
        soma(l,1) = soma(l,1) + ptot(k)
      endif
    enddo
    ! divide by bin volume to get bar; flip sign (lammps compressive convention)
    soma(l,1) = -soma(l,1) / (boxx * boxy * dr)
  enddo

  ! integrate with trapezoidal rule
  dpdr(1) = 0.0d0
  do i = 2, nbins
    dpdr(1) = dpdr(1) + dr * (soma(i,1) + soma(i-1,1)) * 0.5d0
  enddo

  ! convert to mn/m:
  !   bar·å × 1e5 (pa/bar) × 1e-10 (m/å) × 1e3 (mn/m) × 0.5 (two interfaces)
  !   = bar·å × 1e-2
  sigma(1)  = dpdr(1) * 0.5d0 * 1.0d-2
  sigmabig  = sigma(1)

  write(*,'(a6,2a16)') 'frame', 'gamma(mn/m)', 'gamma_avg(mn/m)'
  write(*,'(i6,2f16.5)') 1, sigma(1), sigmabig

  ! ---------------------------------------------------------------------------
  ! main loop: remaining production frames
  ! ---------------------------------------------------------------------------
  do b = 2, nconf

    read(9,*) dummyint   ! item: timestep
    read(9,*) dummyint   ! timestep value
    read(9,*) dummyint   ! item: number of atoms
    read(9,*) dummyint   ! natoms (same system)
    read(9,*) dummyint   ! item: box bounds
    read(9,*) lox, hix
    read(9,*) loy, hiy
    read(9,*) loz, hiz
    read(9,*) dummyint   ! item: atoms

    boxx = hix - lox
    boxy = hiy - loy
    boxz = hiz - loz

    do j = 1, natoms
      read(9,*) nnumber(j), nntype(j), nmol_dummy, &
                xu_dummy, yu_dummy, zz(j), &
                pxx_a(j), pyy_a(j), pzz_a(j)
      zz(j)   = zz(j) - loz
      ptot(j) = pzz_a(j) - 0.5d0 * (pxx_a(j) + pyy_a(j))
    enddo

    ! bin pressure anisotropy
    do l = 1, nbins
      soma(l,b) = 0.0d0
      do k = 1, natoms
        if (zz(k) >= (l-1)*dr .and. zz(k) < l*dr) then
          soma(l,b) = soma(l,b) + ptot(k)
        endif
      enddo
      soma(l,b) = -soma(l,b) / (boxx * boxy * dr)
    enddo

    ! trapezoidal integration
    dpdr(b) = 0.0d0
    do i = 2, nbins
      dpdr(b) = dpdr(b) + dr * (soma(i,b) + soma(i-1,b)) * 0.5d0
    enddo

    sigma(b)  = dpdr(b) * 0.5d0 * 1.0d-2
    sigmabig  = sigmabig + sigma(b)

    write(*,'(i6,2f16.5)') b, sigma(b), sigmabig / dble(b)

  enddo   ! b (frame loop)

  close(9)

  ! ---------------------------------------------------------------------------
  ! final statistics
  ! ---------------------------------------------------------------------------
  gamma_mean = sigmabig / dble(nconf)

  std = 0.0d0
  do i = 1, nconf
    std = std + (sigma(i) - gamma_mean)**2
  enddo
  std = sqrt(std / dble(nconf - 1))   ! sample standard deviation
  ste = std / sqrt(dble(nconf))

  call system_clock(t1_count)
  wall_end = dble(t1_count) / dble(count_rate_int)

  write(*,'(a,f10.3,a)') ' wall time: ', wall_end - wall_start, ' s'
  write(*,*)
  write(*,'(3a16)') 'gamma(mn/m)', 'std_dev', 'std_err'
  write(*,'(3f16.5)') gamma_mean, std, ste

  ! ---------------------------------------------------------------------------
  ! write time-averaged pressure profile p(z) to file
  ! ---------------------------------------------------------------------------
  open(47, file=outputfile, status='replace')
  write(47,'(a)') '# z(a)    p(z)(bar)'

  temp1 = 0.0d0
  do l = 1, nbins
    do b = 1, nconf
      temp1(l) = temp1(l) + soma(l,b)
    enddo
    temp1(l) = temp1(l) / dble(nconf)
    write(47,'(2f14.6)') (l-1)*dr + 0.5d0*dr, temp1(l)   ! bin centre
  enddo

  close(47)

  ! ---------------------------------------------------------------------------
  ! cleanup
  ! ---------------------------------------------------------------------------
  deallocate(nnumber, nntype, zz, pxx_a, pyy_a, pzz_a, ptot)
  deallocate(soma, temp1, sigma, dpdr)

end program stension_final
