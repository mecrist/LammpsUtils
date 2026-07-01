! ===========================================================================
! dump2xyz.f90 — convert LAMMPS dump to multi-frame XYZ format
! ===========================================================================
!
! DUMP FORMAT EXPECTED (unified dump, shared with all other tools, except S_tensions):
!   dump  d1 all custom 1000 traj.lammpstrj id type mol xu yu zu
!   dump_modify d1 sort id
!   Columns: id  type  mol  xu  yu  zu
!
! *** INPUT DATA (marked USER EDIT) ***
!   1. nconfs             — number of frames to convert
!   2. SELECT CASE block  — your LAMMPS type numbers mapping to element symbols

program dump2xyz

IMPLICIT NONE

CHARACTER(255) :: inputfile, outputfile

CHARACTER(3) :: atomtype
INTEGER :: jj, kk
INTEGER :: ntype, nnumber, mol   ! LAMMPS type, atom id, molecule id
REAL*8  :: x, y, z              ! coordinates

INTEGER :: natoms, nconfs, dummyint
CHARACTER(255) :: dummychar

REAL*8 :: xlo, xhi, ylo, yhi, zlo, zhi

! ===========================================================================
! USER EDIT — Number of frames to convert ===================================
nconfs = 100

CALL getarg(1, inputfile)
CALL getarg(2, outputfile)

!reads firsst frame
OPEN(9, FILE=trim(inputfile), STATUS='old')
READ(9,*) dummychar          ! ITEM: TIMESTEP
READ(9,*) dummyint           ! timestep value
READ(9,*) dummychar          ! ITEM: NUMBER OF ATOMS
READ(9,*) natoms
READ(9,*) dummychar          ! ITEM: BOX BOUNDS
READ(9,*) xlo, xhi
READ(9,*) ylo, yhi
READ(9,*) zlo, zhi
READ(9,*) dummychar          ! ITEM: ATOMS

! skiping atom data of first frame
do jj = 1, natoms
  READ(9,*) nnumber, ntype, mol, x, y, z
enddo
CLOSE(9)

! read and convert each frame
OPEN(9,  FILE=trim(inputfile),  STATUS='old')
OPEN(10, FILE=trim(outputfile))

DO kk = 1, nconfs

  write(*,'(A,I6,A,I6)') 'Converting frame ', kk, ' / ', nconfs

  READ(9,*) dummychar          ! ITEM: TIMESTEP
  READ(9,*) dummyint           ! timestep value
  READ(9,*) dummychar          ! ITEM: NUMBER OF ATOMS
  READ(9,*) natoms
  READ(9,*) dummychar          ! ITEM: BOX BOUNDS
  READ(9,*) xlo, xhi
  READ(9,*) ylo, yhi
  READ(9,*) zlo, zhi
  READ(9,*) dummychar          ! ITEM: ATOMS

  !write XYZ frame header
  WRITE(10,*) natoms
  WRITE(10,*) (xhi-xlo), (yhi-ylo), (zhi-zlo)   ! box lengths for PBC

  !Read and write each atom
  DO jj = 1, natoms

    READ(9,*) nnumber, ntype, mol, x, y, z

    ! =========================================================================
    ! USER EDIT — Atom type to element symbol mapping=========================
    ! replace the CASE values with your own LAMMPS type numbers===============
    ! =========================================================================
    SELECT CASE (ntype)
      CASE(1,4,6,8,10,12,23,31,32,37)
        atomtype = 'C  '
      CASE(5,7,9,11,13,15,26,28,35,38)
        atomtype = 'H  '
      CASE(2,14,19,25,27,29,30,33,34)
        atomtype = 'O  '
      CASE(3)
        atomtype = 'N  '
      CASE(16)
        atomtype = 'Cl '
      CASE(17)
        atomtype = 'Na '
      CASE(18,24,36)
        atomtype = 'S  '
      CASE(20)
        atomtype = 'Mg '
      CASE(21)
        atomtype = 'Ca '
      CASE(22)
        atomtype = 'K  '
      CASE DEFAULT
        !unknown type — write 'X' so the file is still valid XYZ
        atomtype = 'X  '
        write(*,'(A,I4,A)') 'WARNING: unknown atom type ', ntype, ' written as X'
    END SELECT
    ! =========================================================================

    WRITE(10,'(A3,3F12.4)') atomtype, x, y, z

  ENDDO

ENDDO  

CLOSE(9)
CLOSE(10)

write(*,*) 'Done. Output written to ', trim(outputfile)

STOP
END program dump2xyz
