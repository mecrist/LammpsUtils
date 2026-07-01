! ===========================================================================
! hbonds-openmp.f90 — hydrogen bond count resolved by z-bin
! ===========================================================================
!
! counts the average number of hydrogen bonds per water molecule as a
! function of z-position, accumulated over all frames.
!
! geometric h-bond criterion (donor-h···acceptor):
!   O-O distance  <= oodistcut   (default 4.0 å)
!   O-H-O angle   >= thetacut   (default 120°)
!
! the H atom that is closer to the acceptor o is selected for the angle.
!
! the outer loop over donor oxygens is parallelised with openmp.
!
! DUMP format expected (unified dump, columns: id type mol x y z):
!
! THINGS YOU MIUST EDIT(marked user edit):
!   nstep, wateroxygentype, waterhydrogentype, oodistcut, thetacut_deg, binsize
!
program hbonds
  use omp_lib
  implicit none

  ! ---------------------------------------------------------------------------
  ! user edit — analysis parameters
  ! ---------------------------------------------------------------------------
  integer,            parameter :: nstep            = 700
  integer,            parameter :: wateroxygentype  = 14
  integer,            parameter :: waterhydrogentype = 15
  double precision,   parameter :: oodistcut        = 4.0d0    ! A
  double precision,   parameter :: thetacut_deg     = 120.0d0  ! degrees
  double precision,   parameter :: binsize          = 1.0d0    ! A
  character(len=100), parameter :: inputtraj        = 'traj.lammpstrj'
  ! ---------------------------------------------------------------------------

  ! atomic data
  integer,          allocatable :: id(:), itype(:), mol(:)
  double precision, allocatable :: pos(:,:)   ! pos(3,natoms)

  ! box bounds (updated every frame)
  double precision :: lox, loy, loz, hix, hiy, hiz
  double precision :: boxx, boxy, boxz

  ! distances and angle
  double precision :: dist
  double precision :: distx, disty, distz
  double precision :: d1x, d1y, d1z, d1
  double precision :: d2x, d2y, d2z, d2
  double precision :: v1(3), v2(3), v1m, v2m, v1v2, costheta, theta

  ! H-bond accumulation (grand total across all frames)
  integer, allocatable :: h2obinnum(:)   ! total water-O count per bin
  integer, allocatable :: hbbinnum(:)    ! total HB count per bin

  ! per-frame accumulators (reduced by openmp)
  integer, allocatable :: h2obin_frame(:)
  integer, allocatable :: hbbin_frame(:)

  integer :: nbins
  integer :: i, j, k, l, hmol, h1, h2
  integer :: natoms, h2onum, hbnum_frame
  double precision :: thetacut   ! angle cutoff in radians

  ! timing
  integer          :: t0_count, t1_count, count_rate_int, count_max_int
  double precision :: wall_start, wall_end

  ! ---------------------------------------------------------------------------
  ! derived parameters
  ! ---------------------------------------------------------------------------
  thetacut = thetacut_deg * 0.017453292519943d0   ! degrees to radians

  ! ---------------------------------------------------------------------------
  ! read first frame: get natoms, box, count water oxygens
  ! ---------------------------------------------------------------------------
  open(9, file=inputtraj, status='old')

  read(9,*)   ! item: timestep
  read(9,*)   ! timestep value
  read(9,*)   ! item: number of atoms
  read(9,*) natoms
  read(9,*)   ! item: box bounds
  read(9,*) lox, hix
  read(9,*) loy, hiy
  read(9,*) loz, hiz
  read(9,*)   ! item: atoms

  allocate(pos(3,natoms), id(natoms), itype(natoms), mol(natoms))

  h2onum = 0
  do i = 1, natoms
    read(9,*) id(i), itype(i), mol(i), pos(1,i), pos(2,i), pos(3,i)
    if (itype(i) == wateroxygentype) h2onum = h2onum + 1
  enddo

  rewind(9)

  ! ---------------------------------------------------------------------------
  ! set up z-bins
  ! ---------------------------------------------------------------------------
  nbins = 1 + int((hiz - loz) / binsize)

  allocate(h2obinnum(nbins),    hbbinnum(nbins))
  allocate(h2obin_frame(nbins), hbbin_frame(nbins))

  h2obinnum = 0
  hbbinnum  = 0

  write(*,'(I8,A)') h2onum, ' water molecules found in first frame'
  write(*,'(I8,A)') nbins,  ' z-bins'

  ! ----------------------------------------------------------------------------
  ! main loop over frames
  ! ---------------------------------------------------------------------------
  do i = 1, nstep

    call system_clock(t0_count, count_rate_int, count_max_int)
    wall_start = dble(t0_count) / dble(count_rate_int)

    write(*,'(A,I6)') 'frame ', i

    ! read frame header
    read(9,*)   ! item: TIMESTEP
    read(9,*)   ! timestep value
    read(9,*)   ! item: number of atoms
    read(9,*) natoms
    read(9,*)   ! item: box bounds
    read(9,*) lox, hix
    read(9,*) loy, hiy
    read(9,*) loz, hiz
    read(9,*)   ! item: atoms

    boxx = hix - lox
    boxy = hiy - loy
    boxz = hiz - loz

    ! read atomic positions; shift z so bin index is measured from loz
    do j = 1, natoms
      read(9,*) id(j), itype(j), mol(j), pos(1,j), pos(2,j), pos(3,j)
      pos(3,j) = pos(3,j) - loz
    enddo

    ! reset per-frame accumulators
    h2obin_frame = 0
    hbbin_frame  = 0
    hbnum_frame  = 0

    ! -------------------------------------------------------------------------
    ! OPENMP loop over all donor oxygens (j)
    ! -------------------------------------------------------------------------
    !$OMP PARALLEL DO                                            &
    !$OMP& DEFAULT(NONE)                                         &
    !$OMP& REDUCTION(+:hbbin_frame, h2obin_frame, hbnum_frame)  &
    !$OMP& PRIVATE(j,k,l,h1,h2,hmol)                           &
    !$OMP& PRIVATE(dist,distx,disty,distz)                      &
    !$OMP& PRIVATE(d1x,d1y,d1z,d1,d2x,d2y,d2z,d2)             &
    !$OMP& PRIVATE(v1,v2,v1m,v2m,v1v2,costheta,theta)          &
    !$OMP& SHARED(natoms,itype,mol,pos,thetacut,oodistcut)      &
    !$OMP& SHARED(wateroxygentype,waterhydrogentype)            &
    !$OMP& SHARED(nbins,binsize,boxx,boxy,boxz)
    do j = 1, natoms

      if (itype(j) /= wateroxygentype) cycle   ! j must be a donor oxygen

      ! bin this donor o for the denominator
      block
        integer :: m
        do m = 1, nbins
          if (pos(3,j) >= (m-1)*binsize .and. pos(3,j) < m*binsize) then
            h2obin_frame(m) = h2obin_frame(m) + 1
            exit
          endif
        enddo
      end block

      ! loop over potential acceptor oxygens (k)
      do k = 1, natoms

        if (itype(k) /= wateroxygentype) cycle
        if (k == j) cycle

        ! O-O distance with minimum image convention (mic)
        distx = pos(1,j) - pos(1,k)
        distx = distx - boxx * anint(distx / boxx)
        disty = pos(2,j) - pos(2,k)
        disty = disty - boxy * anint(disty / boxy)
        distz = pos(3,j) - pos(3,k)
        distz = distz - boxz * anint(distz / boxz)
        dist  = sqrt(distx**2 + disty**2 + distz**2)

        if (dist > oodistcut) cycle   ! O-O distance criterion not met

        ! find the two H atoms that belong to the donor (mol j)
        ! and measure their distance to the acceptor O (k)
        hmol = 0
        h1   = 0
        h2   = 0
        d1   = 1.0d10
        d2   = 1.0d10

        do l = 1, natoms
          if (mol(l) /= mol(j)) cycle
          if (itype(l) /= waterhydrogentype) cycle

          hmol = hmol + 1

          if (hmol == 1) then
            h1  = l
            d1x = pos(1,k) - pos(1,l)
            d1x = d1x - boxx * anint(d1x / boxx)
            d1y = pos(2,k) - pos(2,l)
            d1y = d1y - boxy * anint(d1y / boxy)
            d1z = pos(3,k) - pos(3,l)
            d1z = d1z - boxz * anint(d1z / boxz)
            d1  = sqrt(d1x**2 + d1y**2 + d1z**2)

          else if (hmol == 2) then
            h2  = l
            d2x = pos(1,k) - pos(1,l)
            d2x = d2x - boxx * anint(d2x / boxx)
            d2y = pos(2,k) - pos(2,l)
            d2y = d2y - boxy * anint(d2y / boxy)
            d2z = pos(3,k) - pos(3,l)
            d2z = d2z - boxz * anint(d2z / boxz)
            d2  = sqrt(d2x**2 + d2y**2 + d2z**2)

          endif
        enddo   ! l (hydrogen search)

        ! need exactly 2 H atoms on the donor; skip otherwise
        if (hmol /= 2) cycle

        ! select the h closer to the acceptor O (k) for the angle calculation.
        ! vectors are: v1 = h -> donor-o(j),  v2 = h -> acceptor-o(k)
        ! this gives the o-h···o angle at h.
        if (d1 <= d2) then
          ! h1 is closer to acceptor
          v1(1) = pos(1,j) - pos(1,h1)
          v1(1) = v1(1) - boxx * anint(v1(1) / boxx)
          v1(2) = pos(2,j) - pos(2,h1)
          v1(2) = v1(2) - boxy * anint(v1(2) / boxy)
          v1(3) = pos(3,j) - pos(3,h1)
          v1(3) = v1(3) - boxz * anint(v1(3) / boxz)

          v2(1) = pos(1,k) - pos(1,h1)
          v2(1) = v2(1) - boxx * anint(v2(1) / boxx)
          v2(2) = pos(2,k) - pos(2,h1)
          v2(2) = v2(2) - boxy * anint(v2(2) / boxy)
          v2(3) = pos(3,k) - pos(3,h1)
          v2(3) = v2(3) - boxz * anint(v2(3) / boxz)
        else
          ! h2 is closer to acceptor
          v1(1) = pos(1,j) - pos(1,h2)
          v1(1) = v1(1) - boxx * anint(v1(1) / boxx)
          v1(2) = pos(2,j) - pos(2,h2)
          v1(2) = v1(2) - boxy * anint(v1(2) / boxy)
          v1(3) = pos(3,j) - pos(3,h2)
          v1(3) = v1(3) - boxz * anint(v1(3) / boxz)

          v2(1) = pos(1,k) - pos(1,h2)
          v2(1) = v2(1) - boxx * anint(v2(1) / boxx)
          v2(2) = pos(2,k) - pos(2,h2)
          v2(2) = v2(2) - boxy * anint(v2(2) / boxy)
          v2(3) = pos(3,k) - pos(3,h2)
          v2(3) = v2(3) - boxz * anint(v2(3) / boxz)
        endif

        ! compute angle at h
        v1m     = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
        v2m     = sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2)
        v1v2    = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        costheta = v1v2 / (v1m * v2m)
        ! clamp to [-1,1] to guard against floating-point rounding
        costheta = max(-1.0d0, min(1.0d0, costheta))
        theta    = acos(costheta)

        if (theta < thetacut) cycle   ! angle criterion not met

        ! count this h-bond in the donor's bin
        block
          integer :: m
          do m = 1, nbins
            if (pos(3,j) >= (m-1)*binsize .and. pos(3,j) < m*binsize) then
              hbbin_frame(m) = hbbin_frame(m) + 1
              exit
            endif
          enddo
        end block

        hbnum_frame = hbnum_frame + 1

      enddo   ! k (acceptor search)
    enddo     ! j (donor search)
    !$OMP END PARALLEL DO

    ! accumulate into grand totals
    h2obinnum = h2obinnum + h2obin_frame
    hbbinnum  = hbbinnum  + hbbin_frame

    ! per-frame summary
    write(*,'(A,I8,A,I8,A,F10.6,A)') &
      '  hbonds= ', hbnum_frame, '  water-o= ', h2onum, &
      '  hb/mol= ', dble(hbnum_frame)/dble(h2onum), ' (this frame)'

    call system_clock(t1_count)
    wall_end = dble(t1_count) / dble(count_rate_int)
    write(*,'(A,I6,A,F10.3,A)') '  frame ', i, '  wall: ', wall_end - wall_start, ' s'

  enddo   ! i (frame loop)

  ! ---------------------------------------------------------------------------
  ! print grand-average profile
  ! ---------------------------------------------------------------------------
  write(*,'(/,A)') '--- grand-average h-bond profile ---'
  write(*,'(A5,A14,A10,A10)') 'bin', 'hb/mol', 'n_water', 'n_hb'
  do k = 1, nbins
    if (h2obinnum(k) > 0) then
      write(*,'(I5,F14.6,2I10)') k, dble(hbbinnum(k))/dble(h2obinnum(k)), &
                                     h2obinnum(k), hbbinnum(k)
    else
      write(*,'(I5,A,2I10)') k, '      0.000000', 0, 0
    endif
  enddo

  deallocate(pos, id, itype, mol)
  deallocate(h2obinnum, hbbinnum, h2obin_frame, hbbin_frame)

end program hbonds
