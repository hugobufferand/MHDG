MODULE domain_decomposition_module
  USE MPI_OMP
  USE globals
  USE adaptivity_common_module, ONLY: unique_1D, unique_stable, quicksort_int
  USE preprocess, ONLY: GetFaces_mod
  IMPLICIT NONE

  ! Partitioning methods (divType argument of split_mesh)
  !   1, 2, 3 : geometric strips (radial / X / Y), rotated 5 deg until valid
  !   4       : graph partitioning of the element connectivity with SCOTCH
  !             (falls back to Y strips if SCOTCH is unavailable or fails)
  INTEGER, PARAMETER :: DD_STRIPS_Y = 3
  INTEGER, PARAMETER :: DD_SCOTCH   = 4

CONTAINS

  !> Partitioning method to use: SCOTCH when compiled with -DWITH_SCOTCH,
  !> geometric Y strips otherwise. The environment variable
  !> MHDG_PARTITIONER=strips|scotch overrides the default at run time.
  INTEGER FUNCTION default_partition_type() RESULT(divtype)
    CHARACTER(LEN=32) :: env
    INTEGER           :: env_len, env_status

#ifdef WITH_SCOTCH
    divtype = DD_SCOTCH
#else
    divtype = DD_STRIPS_Y
#endif
    CALL GET_ENVIRONMENT_VARIABLE('MHDG_PARTITIONER', env, env_len, env_status)
    IF (env_status .EQ. 0 .AND. env_len .GT. 0) THEN
       SELECT CASE (TRIM(ADJUSTL(env)))
       CASE ('strips', 'STRIPS')
          divtype = DD_STRIPS_Y
       CASE ('scotch', 'SCOTCH')
          divtype = DD_SCOTCH
       END SELECT
    ENDIF
  END FUNCTION default_partition_type


  SUBROUTINE split_mesh(ndiv, divtype, verbose)

    INTEGER, INTENT(IN)                     :: nDiv, divType
    LOGICAL, INTENT(IN)                     :: verbose

    TYPE DD_mesh
       INTEGER, ALLOCATABLE                 :: TT(:,:), intfaces_loc(:,:), extfaces_loc(:,:), TTb_CON(:,:), Tb(:,:)
       REAL*8, ALLOCATABLE                  :: XX(:,:)
       INTEGER, ALLOCATABLE                 :: loc2glob_fa(:), loc2glob_el(:), reptFaces(:), ghostFaces(:), ghostElems(:), loc2glob_no(:), ghostFlp(:),&
            & mapextfa_loc(:), Tb_sizes(:), GhelsPro(:), GhelsLoc(:), GhostPro(:), GhostLoc(:), boundaryFlag(:)
       INTEGER                              :: Nfaces, Nnodes, Nelems, Ndim, Nnodesperelem, Nnodesperface, Nintfaces, NextFaces, elemType
    ENDTYPE DD_mesh

    TYPE(DD_mesh)                           :: M(nDiv)

    INTEGER, POINTER, DIMENSION(:,:)        :: intfaces_glob => NULL(), extfaces_glob => NULL(), T => NULL(), TTb => NULL()
    INTEGER, POINTER, DIMENSION(:)          :: boundaryFlag => NULL()
    REAL*8, POINTER, DIMENSION(:,:)         :: X => NULL()

    INTEGER                                 :: elemType, Nvert, n_ext_faces, n_int_faces, Nelems, Nfaces, Nnodes, Nnodes_mean, Nelems_mean, Nel_tot
    INTEGER                                 :: Nel_loc, n_boundaries
    INTEGER                                 :: iel, it, i , j, ii, jj, counter, clg, offset, ifa, iface, ib, n, ierr, b1, b2, n1, n2, s
    INTEGER, PARAMETER                      :: imax = 1000
    INTEGER                                 :: Efaces(refElPol%Nfaces, refElPol%Nfacenodeslin)
    INTEGER, ALLOCATABLE, DIMENSION(:)      :: indgh, innext, ghel, ghel_temp, invmap, loc2glob_fa, loc2glob_el, reptFaces, ghostFaces, faceglo, nodes, unique_boundaryFlag, Tb_prov_unique, mapextf_unique, checkFaces
    INTEGER, ALLOCATABLE, DIMENSION(:)      :: nodesfix, mapTfix, Tb_sizes, mapextf, flipFace, indices, indices2, auxext, auxflp, boundaryFlag_loc
    INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: auxgh, Tprov, Tfix, Tp, Tbprov_temp, intfaces, extfaces, Tb_bound, Tb_prov, T_CON, Tb_loc, loc2glob_faces
    REAL*8, ALLOCATABLE, DIMENSION(:,:)     :: Xprov
    LOGICAL, ALLOCATABLE, DIMENSION(:)      :: ind
    REAL*8                                  :: xcenter, ycenter, step, theta0, theta1, angle, tol = 1e-6, pi_loc
    REAL*8, ALLOCATABLE                     :: xb(:), yb(:), thetab(:), theta(:)
    LOGICAL                                 :: redo, check
    REAL*8                                  :: RotMat(2,2)
    ! partition of the global elements (1..nDiv) and face ownership
    INTEGER, ALLOCATABLE                    :: part(:), tbidx(:)
    INTEGER                                 :: dtype, owner, gfa, nbal
    LOGICAL                                 :: part_ok

    NULLIFY(intfaces_glob, extfaces_glob, T, TTb, boundaryFlag, X)
    dtype = divType

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "*************************************************"
       WRITE(*,*) "            DOMAIN DECOMPOSITION                 "
       WRITE(*,*) "*************************************************"
    ENDIF

    pi_loc = 4.D0*DATAN(1.D0)

    elemType     = refElPol%elemType
    Nvert        = refElPol%Nvertices
    n_ext_faces  = Mesh%Nextfaces
    n_int_faces  = Mesh%Nintfaces

    intfaces_glob => Mesh%intfaces
    extfaces_glob => Mesh%extfaces

    Nelems = Mesh%Nelems
    Nfaces = Mesh%Nfaces
    Nnodes = Mesh%Nnodes

    Nnodes_mean = FLOOR(REAL(Nnodes/ndiv))
    Nelems_mean = FLOOR(REAL(Nelems/ndiv))

    ! rescale back to physical values
    Mesh%X = Mesh%X*phys%lscale

    X => Mesh%X
    T => Mesh%T

    ! Apply shift back if axisymmetric case
    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. MINVAL(Mesh%X(:,1)) < tol)) THEN
       X(:, 1) = X(:, 1) - geom%R0
    END IF

    boundaryFlag => Mesh%boundaryFlag
    TTb => Mesh%Tb

    CALL UNIQUE_STABLE(boundaryFlag,unique_boundaryFlag)
    n_boundaries = SIZE(unique_boundaryFlag)

    ALLOCATE(Tb_sizes(n_boundaries))
    Tb_sizes = 0
    counter = 1

    DO i = 1, n_boundaries
       Tb_sizes(i) = COUNT(boundaryFlag .EQ. unique_boundaryFlag(i))
    ENDDO

    angle = 0
    redo = .TRUE.

    xcenter = 0.5*(MAXVAL(X(:,1))+MINVAL(X(:,1)))
    ycenter = 0.5*(MAXVAL(X(:,2))+MINVAL(X(:,2)))

    ALLOCATE(xb(Nelems))
    ALLOCATE(yb(Nelems))
    ALLOCATE(thetab(Nelems))
    ALLOCATE(theta(nDiv+1))
    ALLOCATE(auxgh(Nfaces,2))
    ALLOCATE(indgh(Nfaces))
    ALLOCATE(ind(Nelems))

    xb = 0.
    yb = 0.
    thetab = 0.
    theta = 0.
    auxgh = 0
    indgh = 0
    ind = .FALSE.

    ALLOCATE(part(Nelems))
    part = 0
    IF (dtype .EQ. DD_SCOTCH) THEN
       CALL graph_partition(nDiv, Nelems, intfaces_glob, part, part_ok)
       IF (.NOT. part_ok) THEN
          IF (MPIvar%glob_id .EQ. 0) THEN
             WRITE(*,*) "Graph partitioning unavailable or failed: using Y strips."
          ENDIF
          dtype = DD_STRIPS_Y
       ENDIF
    ENDIF

    ! keep on rotating the mesh if no good partition is found
    DO WHILE(redo)

       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

       IF (dtype .NE. DD_SCOTCH) THEN
          IF (angle .GE. 360.) THEN
             IF(MPIvar%glob_id .EQ. 0) THEN
                WRITE(*,*) "Domain decomposition: no valid strip partition found for any rotation."
                WRITE(*,*) "Use fewer MPI processes or the SCOTCH graph partitioner. STOP."
             ENDIF
             CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          ENDIF
          IF(MPIvar%glob_id .EQ. 0) THEN
             WRITE(*,*) "Angle of rotation: ", angle
          ENDIF
       ENDIF


       !*********************** Mesh's info ************************!


       xb = 0.
       yb = 0.
       thetab = 0.
       theta = 0.
       auxgh = 0
       indgh = 0
       ind = .FALSE.

       DO iel = 1, Nelems
          xb(iel) = SUM(X(T(iel,1:Nvert),1))/REAL(Nvert)
          yb(iel) = SUM(X(T(iel,1:Nvert),2))/REAL(Nvert)
       ENDDO

       SELECT CASE (dtype)
       CASE (DD_SCOTCH)
          ! partition already computed by graph_partition
       CASE (1)
          ! Radial division
          DO iel = 1,Nelems
             thetab(iel) = ATAN(yb(iel),xb(iel))
             IF (thetab(iel) .LT. 0) THEN
                thetab(iel) = thetab(iel) + 2*pi_loc
             ENDIF
          ENDDO
       CASE (2)
          ! X division
          thetab = xb
          step = (MAXVAL(X(:,1)) - MINVAL(X(:,1))) / REAL(nDiv)
          theta = [(MINVAL(X(:,1)) + (i - 1) * step, i = 1, nDiv + 1)]
       CASE (3)
          ! Y division
          thetab = yb
          step = (MAXVAL(X(:,2)) - MINVAL(X(:,2))) / REAL(nDiv)
          theta = [(MINVAL(X(:,2)) + (i - 1) * step, i = 1, nDiv + 1)]
       CASE DEFAULT
          WRITE(*,*) "Option not supported for divType, domain decompositon. STOP."
          STOP
       END SELECT


       !******************** Balancing Division *********************!

       ! !***** Division Loop->ndiv-1 *****!
       !
       IF (dtype .EQ. DD_SCOTCH) THEN
          redo = .FALSE.
       ELSE

       nbal = 0
       DO WHILE(redo)

          nbal = nbal + 1
          redo = .FALSE.
          Nel_tot = 0
          Nel_loc = 0

          ! Balance equation
          DO it = 1, nDiv - 1
             theta0 = theta(it)
             theta1 = theta(it+2)

             DO i = 1, imax
                Nel_loc = 0
                DO j = 1, SIZE(thetab)
                   IF((thetab(j) .GT. theta(it)) .AND. (thetab(j) .LT. theta(it+1))) THEN
                      Nel_loc = Nel_loc + 1
                   ENDIF
                ENDDO

                IF((Nel_loc-Nelems_mean) .GT. 5) THEN
                   theta1 = theta(it+1)
                   theta(it+1) = 0.5*(theta0+theta1)
                ELSEIF((Nel_loc-Nelems_mean) .LT. -5) THEN
                   theta0 = theta(it+1)
                   theta(it+1) = 0.5*(theta0+theta1)
                ELSE
                   ! partitions are well balanced
                   EXIT
                ENDIF
             ENDDO

             IF(i .EQ. imax+1) THEN
                !WRITE(*,*) "redo"
                redo = .TRUE.
             ENDIF
             Nel_tot = Nel_tot + Nel_loc
          ENDDO
          IF (redo .AND. nbal .GE. 10) EXIT   ! give up this angle, rotate below
       ENDDO

       Nel_tot = 0
       DO it = 1, SIZE(theta)-1
          Nel_loc = 0
          DO i = 1, SIZE(thetab)
             IF((thetab(i) .GT. theta(it)) .AND. (thetab(i) .LT. theta(it+1))) THEN
                Nel_loc = Nel_loc + 1
             ENDIF
          ENDDO
          IF(verbose) THEN
             IF(MPIvar%glob_id .EQ. 0) THEN
                WRITE(*,*) "Number of elements per division:", Nelems_mean, Nel_loc
             ENDIF
          ENDIF

          Nel_tot = Nel_tot + Nel_loc
       ENDDO

       IF(Nel_tot .NE. Nelems) THEN
          IF(verbose) THEN
             IF(MPIvar%glob_id .EQ. 0) THEN
                WRITE(*,*) 'Error in numbering of elements: ', Nel_tot, Nelems
             ENDIF
          ENDIF
          redo = .TRUE.
       ENDIF

#ifdef PARALL
       ! Use MPI_Allreduce to check if any process encountered the exit condition
       CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

       IF(redo) THEN
          ! no balanced strip cut: rotate the mesh by 5 deg and retry
          X(:,1) = X(:,1) - xcenter
          X(:,2) = X(:,2) - ycenter
          RotMat(1,1) = COS(5.*pi_loc/180.)
          RotMat(1,2) = -SIN(5.*pi_loc/180.)
          RotMat(2,1) = SIN(5.*pi_loc/180.)
          RotMat(2,2) = COS(5.*pi_loc/180.)
          X = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(X)))
          X(:,1) = X(:,1) + xcenter
          X(:,2) = X(:,2) + ycenter
          angle = angle + 5.
          CYCLE
       ENDIF

       ! strip index of each element (all elements are strictly inside a strip here)
       part = 0
       DO it = 1, nDiv
          WHERE ((thetab .GT. theta(it)) .AND. (thetab .LT. theta(it+1))) part = it
       ENDDO

       ENDIF ! strips

       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifndef PARALL
       DO it = 1, nDiv
#else
          DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif
             IF(verbose) THEN
                WRITE(*,*) "Division: ", it
             ENDIF

             ALLOCATE(loc2glob_el(COUNT(part .EQ. it)))
             loc2glob_el = 0
             Tb_sizes = 0

             ! set the index of the global element in the elemental loc2glob
             ind = .FALSE.

             counter = 1
             DO i = 1, Nelems
                IF(part(i) .EQ. it) THEN
                   loc2glob_el(counter) = i
                   counter = counter + 1
                   ind(i) = .TRUE.
                ENDIF
             ENDDO

             IF(verbose) THEN
                WRITE(*,*) "Add ghost element to the division"
             ENDIF

             IF(nDiv .GT. 1) THEN

                auxgh = 0
                DO i = 1, n_int_faces
                   IF (ANY(intFaces_Glob(i, 1) .EQ. loc2glob_el)) auxgh(i,1) = 1
                   IF (ANY(intFaces_Glob(i, 3) .EQ. loc2glob_el)) auxgh(i,2) = 1
                ENDDO

                indgh = 0
                DO i = 1, Nfaces
                   IF (SUM(auxgh(i, :)) .EQ. 1)  indgh(i) = 1
                ENDDO

                ALLOCATE(ghel_temp(2*SIZE(intFaces_Glob)))
                ghel_temp = -1
                counter = 1
                DO i = 1, SIZE(indgh)
                   IF(indgh(i) .NE. 0) THEN
                      IF(intFaces_Glob(i,1) .NE. intFaces_Glob(i,3)) THEN
                         ghel_temp(counter) = intFaces_Glob(i,1)
                         ghel_temp(counter+1) = intFaces_Glob(i,3)
                         counter = counter + 2
                      ELSE
                         ghel_temp(counter) = intFaces_Glob(i,1)
                         counter = counter + 1
                      ENDIF
                   ENDIF
                ENDDO

                DO i = 1, SIZE(ghel_temp)
                   IF(ghel_temp(i) .EQ. -1) ghel_temp(i) = ghel_temp(1)
                ENDDO

                CALL UNIQUE_1D(ghel_temp, ghel)
                DEALLOCATE(ghel_temp)


                counter = 1
                DO i = 1, SIZE(ghel)
                   IF(part(ghel(i)) .GT. it) THEN
                      counter = counter + 1
                   ENDIF
                ENDDO

                ALLOCATE(innext(counter-1))
                innext = 0
                counter = 1
                DO i = 1, SIZE(ghel)
                   IF(part(ghel(i)) .GT. it) THEN
                      innext(counter) = i
                      counter = counter + 1
                   ENDIF
                ENDDO

                counter = 1
                DO i = 1, SIZE(innext)
                   ind(ghel(innext(i))) = .TRUE.
                ENDDO

                DEALLOCATE(loc2glob_el)
                ALLOCATE(loc2glob_el(COUNT(ind .EQV. .TRUE.)))
                loc2glob_el = 0

                counter = 1
                DO i = 1, SIZE(ind)
                   IF(ind(i)) THEN
                      loc2glob_el(counter) = i
                      counter = counter + 1
                   ENDIF
                ENDDO

                ALLOCATE(M(it)%ghostElems(SIZE(loc2glob_el)))
                M(it)%ghostElems = 0  ! Initialize all elements to false

                DO i = 1, SIZE(loc2glob_el)
                   DO j = 1, SIZE(innext)
                      IF (loc2glob_el(i) .EQ. ghel(innext(j))) THEN
                         M(it)%GhostElems(i) = 1  ! Set to true if loc2glob_el(i) is in ghel(innext)
                         EXIT
                      ENDIF
                   ENDDO
                ENDDO

                DEALLOCATE(ghel)
             ELSE
                ALLOCATE(M(it)%ghostElems(Nelems))
                M(it)%ghostElems = 0
             ENDIF

             IF(verbose) THEN
                IF(MPIvar%glob_id .EQ. 0) THEN
                   WRITE(*,*) "Generate X and T for the current division"
                ENDIF
             ENDIF

             ALLOCATE(Tprov(COUNT(ind .EQV. .TRUE.),Mesh%Nnodesperelem))
             ALLOCATE(Tfix(SIZE(Tprov,1), Mesh%Nnodesperelem))
             Tprov = 0
             Tfix = 0

             counter = 1
             DO i = 1, SIZE(ind)
                IF(ind(i)) THEN
                   Tprov(counter,:) = T(i,:)
                   counter = counter + 1
                ENDIF
             ENDDO

             CALL UNIQUE_1D(RESHAPE(Tprov, [SIZE(Tprov,1)*SIZE(Tprov,2)]), nodes)

             ALLOCATE(Xprov(SIZE(nodes),Mesh%ndim))
             Xprov = X(nodes,:)

             ALLOCATE(mapTfix(SIZE(nodes)))
             mapTfix = nodes

             ALLOCATE(invmap(mapTfix(SIZE(mapTfix))))
             invmap = 0
             DO i = 1, SIZE(mapTfix)
                invmap(mapTfix(i)) = i
             ENDDO

             DO i = 1, SIZE(Tprov, 1)
                DO j = 1, SIZE(Tprov, 2)
                   Tfix(i, j) = invmap(Tprov(i, j))
                ENDDO
             ENDDO

             ! nodesfix is allocated in here
             CALL UNIQUE_1D(RESHAPE(Tfix, [SIZE(Tfix,1)*SIZE(Tfix,2)]),nodesfix)

             ! check that everything is ok so far...
             IF (ALL(nodesfix .NE. (/(j, j=1,MAXVAL(Tfix))/))) THEN
                IF(verbose) THEN
                   WRITE(*,*) 'Problem in fixed connectivity. Retrying.'
                ENDIF
                redo = .TRUE.
                !EXIT
             ENDIF

#ifdef PARALL
             ! Use MPI_Allreduce to check if any process encountered the exit condition
             CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
             IF(redo) THEN
                EXIT
             ENDIF

             IF(verbose) THEN
                IF(MPIvar%glob_id .EQ. 0) THEN
                   WRITE(*,*) "Create loc2glob_fa: interior faces"
                ENDIF
             ENDIF

             ! Get faces takes ages!!
             CALL GetFaces_mod(Tfix(:,1:Nvert), intfaces, extFaces)

             ALLOCATE(M(it)%intfaces_loc(SIZE(intfaces,1),SIZE(intfaces,2)))
             ALLOCATE(M(it)%extfaces_loc(SIZE(extfaces,1),SIZE(extfaces,2)))

             M(it)%intfaces_loc = intfaces
             M(it)%extfaces_loc = extfaces
             M(it)%Nfaces = SIZE(intfaces,1) + SIZE(extfaces,1)
             M(it)%Nintfaces = SIZE(intfaces,1)
             M(it)%Nextfaces = SIZE(extfaces,1)
             M(it)%Nnodesperface = refElPol%nDeg + 1
             M(it)%elemType = Mesh%elemType

             ALLOCATE(loc2glob_fa(M(it)%Nfaces))
             ALLOCATE(reptFaces(M(it)%Nfaces))
             ALLOCATE(faceglo(5))

             loc2glob_fa = 0
             reptFaces = 0
             faceglo = 0

             DO j = 1, SIZE(intfaces,1)

                faceglo(1) = loc2glob_el(intfaces(j,1))
                faceglo(2) = intfaces(j,2)
                faceglo(3) = loc2glob_el(intfaces(j,3))
                faceglo(4) = intfaces(j,4)
                faceglo(5) = intfaces(j,5)

                DO jj = 1, SIZE(intFaces_Glob,1)
                   IF(SUM(ABS(intFaces_Glob(jj,:)-faceglo)) .EQ. 0) THEN
                      loc2glob_fa(j) = jj
                   ENDIF
                ENDDO

                IF(loc2glob_fa(j) .EQ. 0) THEN
                   WRITE(*,*) "Problem in the creation of loc2glob_fa. STOP"
                   STOP
                ENDIF
             ENDDO
             IF(verbose) THEN
                IF(MPIvar%glob_id .EQ. 0) THEN
                   WRITE(*,*) "Generate boundary connectivity for the current division."
                ENDIF
             ENDIF

             offset = 0
             ALLOCATE(Tb_loc(SIZE(TTb,1),SIZE(TTb,2)))
             ALLOCATE(boundaryFlag_loc(SIZE(TTb,1)))
             boundaryFlag_loc = 0
             Tb_sizes = 0
             Tb_loc = 0
             clg = 1


             DO ib = 1, n_boundaries
                n = COUNT(boundaryFlag .EQ. unique_boundaryFlag(ib))
                ALLOCATE(Tb_bound(n, Mesh%Nnodesperface))
                ALLOCATE(tbidx(n))
                Tb_bound = 0
                tbidx = 0

                counter = 1
                DO i = 1, SIZE(boundaryFlag)
                   IF(boundaryFlag(i) .EQ. unique_boundaryFlag(ib)) THEN
                      Tb_bound(counter,:) = TTb(i, :)
                      tbidx(counter) = i          ! global boundary face (row of Tb / extfaces)
                      counter = counter + 1
                   ENDIF
                ENDDO

                ALLOCATE(Tbprov_temp(n, Mesh%Nnodesperface))
                Tbprov_temp = 0

                counter = 1

                DO ifa = 1, n
                   ! the boundary face belongs to this division if its element does
                   ! (rows of Tb and extfaces are aligned, see preprocess)
                   IF(ind(extfaces_glob(tbidx(ifa),1))) THEN
                      Tbprov_temp(counter,:) = Tb_bound(ifa,:)
                      counter = counter + 1
                      loc2glob_fa(clg+M(it)%Nintfaces) = tbidx(ifa)+n_int_faces
                      clg = clg + 1
                   ENDIF
                ENDDO
                DEALLOCATE(tbidx)

                Tb_sizes(ib) = counter - 1

                IF(counter .NE. 1) THEN
                   ALLOCATE(Tb_prov(counter-1,Mesh%Nnodesperface))
                   Tb_prov = Tbprov_temp(1:counter-1,:)

                   ! Tb_prov_unique is allocated here
                   CALL UNIQUE_1D(RESHAPE(Tb_prov, [SIZE(Tb_prov,1)*SIZE(Tb_prov,2)]), Tb_prov_unique)

                   DO i = 1, SIZE(Tb_prov_unique)
                      DO j = 1, SIZE(nodes)
                         IF(nodes(j) .EQ. Tb_prov_unique(i)) THEN
                            DO ii = 1, SIZE(Tb_prov,1)
                               DO jj = 1, SIZE(Tb_prov,2)
                                  IF(Tb_prov(ii, jj) .EQ. Tb_prov_unique(i)) THEN
                                     Tb_prov(ii, jj) = j
                                  ENDIF
                               ENDDO
                            ENDDO
                            EXIT
                         ENDIF
                      ENDDO

                      IF(j .EQ. SIZE(nodes)+1) THEN
                         IF(verbose) THEN
                            WRITE(*,*) "Line was cut too close to boundary. Retrying."
                         ENDIF
                         redo = .TRUE.
                         EXIT
                      ENDIF
                   ENDDO

                   DEALLOCATE(Tb_prov_unique)

                   IF(redo) THEN
                      EXIT
                   ENDIF

                   Tb_loc(offset+1:offset+SIZE(Tb_prov,1),:) = Tb_prov
                   boundaryFlag_loc(offset+1:offset+SIZE(Tb_prov,1)) = unique_boundaryFlag(ib)
                   offset = offset + SIZE(Tb_prov,1)
                   Tb_sizes(ib) = SIZE(Tb_prov,1)
                   DEALLOCATE(Tb_prov)

                ENDIF
                DEALLOCATE(Tb_bound)
                DEALLOCATE(Tbprov_temp)

             ENDDO

#ifdef PARALL
             ! Use MPI_Allreduce to check if any process encountered the exit condition
             CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
             ! if redo rotate the mesh
             IF(redo) THEN
                X(:,1) = X(:,1) - xcenter
                X(:,2) = X(:,2) - ycenter
                RotMat(1,1) = COS(5.*pi_loc/180.)
                RotMat(1,2) = -SIN(5.*pi_loc/180.)
                RotMat(2,1) = SIN(5.*pi_loc/180.)
                RotMat(2,2) = COS(5.*pi_loc/180.)
                X = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(X)))
                X(:,1) = X(:,1) + xcenter
                X(:,2) = X(:,2) + ycenter
                angle = angle + 5.
             ENDIF

             ! if redo rotate deallocate all allocated structures
             IF(redo) THEN

                DEALLOCATE(loc2glob_fa)
                DEALLOCATE(intfaces)
                DEALLOCATE(extfaces)
                DEALLOCATE(reptFaces)
                DEALLOCATE(Tfix)
                DEALLOCATE(Tprov)
                DEALLOCATE(Xprov)
                DEALLOCATE(nodes)
                DEALLOCATE(loc2glob_el)
                DEALLOCATE(innext)
                DEALLOCATE(mapTfix)
                DEALLOCATE(invmap)
                DEALLOCATE(nodesfix)
                DEALLOCATE(faceglo)

#ifndef PARALL
                DO i = 1, it
                   DEALLOCATE(M(i)%ghostElems)
                   DEALLOCATE(M(i)%intfaces_loc)
                   DEALLOCATE(M(i)%extfaces_loc)
                ENDDO
#else

                DEALLOCATE(M(it)%ghostElems)
                DEALLOCATE(M(it)%intfaces_loc)
                DEALLOCATE(M(it)%extfaces_loc)
#endif


#ifndef PARALL
                DO i = 1, it-1
#else
                   DO i = it, it
#endif
                      IF(ALLOCATED(M(i)%loc2glob_fa))  DEALLOCATE(M(i)%loc2glob_fa)
                      IF(ALLOCATED(M(i)%reptFaces))    DEALLOCATE(M(i)%reptFaces)
                      IF(ALLOCATED(M(i)%TT))           DEALLOCATE(M(i)%TT)
                      IF(ALLOCATED(M(i)%XX))           DEALLOCATE(M(i)%XX)
                      IF(ALLOCATED(M(i)%loc2glob_no))  DEALLOCATE(M(i)%loc2glob_no)
                      IF(ALLOCATED(M(i)%loc2glob_el))  DEALLOCATE(M(i)%loc2glob_el)
                      IF(ALLOCATED(M(i)%Tb))           DEALLOCATE(M(i)%Tb)
                      IF(ALLOCATED(M(i)%boundaryFlag)) DEALLOCATE(M(i)%boundaryFlag)
                      IF(ALLOCATED(M(i)%Tb_sizes))     DEALLOCATE(M(i)%Tb_sizes)
                   ENDDO


                   IF(ALLOCATED(Tb_bound)) DEALLOCATE(Tb_bound)
                   IF(ALLOCATED(Tbprov_temp)) DEALLOCATE(Tbprov_temp)
                   IF(ALLOCATED(Tb_prov)) DEALLOCATE(Tb_prov)
                   IF(ALLOCATED(Tb_loc)) DEALLOCATE(Tb_loc)
                   IF(ALLOCATED(boundaryFlag_loc)) DEALLOCATE(boundaryFlag_loc)
                   EXIT
                ENDIF

                ! store subdivision
                ALLOCATE(M(it)%loc2glob_fa(SIZE(loc2glob_fa)))
                ALLOCATE(M(it)%reptFaces(SIZE(reptFaces)))
                ALLOCATE(M(it)%TT(SIZE(Tfix,1),SIZE(Tfix,2)))
                ALLOCATE(M(it)%XX(SIZE(Xprov,1),SIZE(Xprov,2)))
                ALLOCATE(M(it)%loc2glob_no(SIZE(nodes)))
                ALLOCATE(M(it)%loc2glob_el(SIZE(loc2glob_el)))
                ALLOCATE(M(it)%Tb(SUM(Tb_sizes),SIZE(Tb_loc,2)))
                ALLOCATE(M(it)%boundaryFlag(SUM(Tb_sizes)))
                ALLOCATE(M(it)%Tb_sizes(SIZE(Tb_sizes)))

                M(it)%Loc2Glob_fa   = loc2glob_fa
                M(it)%ReptFaces     = reptFaces
                M(it)%TT            = Tfix
                M(it)%XX            = Xprov
                M(it)%loc2glob_no   = nodes
                M(it)%loc2glob_el   = loc2glob_el
                M(it)%intfaces_loc  = intfaces
                M(it)%extfaces_loc  = extfaces
                M(it)%Tb            = Tb_loc(1:offset,:)
                M(it)%boundaryFlag  = boundaryFlag_loc(1:offset)
                M(it)%Tb_sizes      = Tb_sizes
                M(it)%Nfaces        = SIZE(intfaces,1) + SIZE(extfaces,1)
                M(it)%Nintfaces     = SIZE(intfaces,1)
                M(it)%Ndim          = SIZE(M(it)%XX,2)
                M(it)%Nnodes        = SIZE(M(it)%XX,1)
                M(it)%Nelems        = SIZE(M(it)%TT,1)
                M(it)%Nnodesperelem = SIZE(M(it)%TT,2)

                IF(ALLOCATED(loc2glob_fa)) DEALLOCATE(loc2glob_fa)
                IF(ALLOCATED(Tb_loc)) DEALLOCATE(Tb_loc)
                IF(ALLOCATED(intfaces)) DEALLOCATE(intfaces)
                IF(ALLOCATED(extfaces)) DEALLOCATE(extfaces)
                IF(ALLOCATED(reptFaces)) DEALLOCATE(reptFaces)
                IF(ALLOCATED(Tfix)) DEALLOCATE(Tfix)
                IF(ALLOCATED(Tprov)) DEALLOCATE(Tprov)
                IF(ALLOCATED(Xprov)) DEALLOCATE(Xprov)
                IF(ALLOCATED(nodes)) DEALLOCATE(nodes)
                IF(ALLOCATED(loc2glob_el)) DEALLOCATE(loc2glob_el)
                IF(ALLOCATED(innext)) DEALLOCATE(innext)
                IF(ALLOCATED(mapTfix)) DEALLOCATE(mapTfix)
                IF(ALLOCATED(invmap)) DEALLOCATE(invmap)
                IF(ALLOCATED(nodesfix)) DEALLOCATE(nodesfix)
                IF(ALLOCATED(faceglo)) DEALLOCATE(faceglo)
                IF(ALLOCATED(boundaryFlag_loc)) DEALLOCATE(boundaryFlag_loc)

             ENDDO
          ENDDO

          ! Boundary connectivity for the overlapping faces

          SELECT CASE (refElPol%elemType)
          CASE (0) !triangle
             Efaces = RESHAPE((/1, 2, 3, 2, 3, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
          CASE (1) !quadrilaterals
             Efaces = RESHAPE((/1, 2, 3, 4, 2, 3, 4, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
          CASE (2) !thetrahedra
             !            Efaces = computeNodesFacesTetra(1);
          CASE (3) !hexahedra
             !        Efaces = computeNodesFacesHexa(1);
          END SELECT

#ifndef PARALL
          DO it = 1, nDiv
#else
             DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif
                IF(verbose) THEN
                   WRITE(*,*) "Division: ", it
                ENDIF

                ALLOCATE(Tp(SIZE(M(it)%TT,1),SIZE(M(it)%TT,2)))
                ALLOCATE(loc2glob_fa(SIZE(M(it)%loc2glob_fa)))
                ALLOCATE(loc2glob_el(SIZE(M(it)%loc2glob_el)))
                ALLOCATE(GhostFaces(M(it)%Nfaces))
                ALLOCATE(extfaces(SIZE(M(it)%extfaces_loc,1),SIZE(M(it)%extfaces_loc,2)))

                Tp          = M(it)%TT
                loc2glob_fa = M(it)%loc2glob_fa
                loc2glob_el = M(it)%loc2glob_el
                ghostFaces  = 0
                extfaces    = M(it)%extfaces_loc

                ! one row per border face of the local mesh: bounded by the external faces,
                ! not by the elements (small partitions have more border faces than elements)
                ALLOCATE(T_con(SIZE(extFaces,1),SIZE(refElPol%Face_nodes,2)))
                ALLOCATE(flipFace(SIZE(extFaces,1)))
                ALLOCATE(mapextf(SIZE(extFaces,1)))

                T_con = 0
                flipFace = 0
                mapextf = 0
                clg = 0

                DO j = 1, SIZE(loc2glob_fa)
                   IF(loc2glob_fa(j) .EQ. 0) THEN
                      clg = j
                      EXIT
                   ENDIF
                ENDDO


                counter = 1
                DO ifa = 1, SIZE(extFaces,1)
                   iel = extFaces(ifa,1);
                   iface = extFaces(ifa,2);
                   n1 = Tp(iel,Efaces(iface,1));
                   n2 = Tp(iel,Efaces(iface,2));
                   check = .TRUE.
                   s = 0

                   DO ib = 1, n_boundaries

                      IF(M(it)%Tb_sizes(ib) .NE. 0) THEN
                         ALLOCATE(Tb_bound(M(it)%Tb_sizes(ib), Mesh%Nnodesperface))
                         Tb_bound = M(it)%Tb(s+1:s+M(it)%Tb_sizes(ib),:)

                         IF(ANY(Tb_bound .EQ. n1) .AND. ANY(Tb_bound .EQ. n2)) THEN

                            DO j = 1, SIZE(Tb_bound,1)
                               IF((n1 .EQ. Tb_bound(j,1)) .AND. (n2 .EQ. Tb_bound(j,SIZE(Tb_bound,2)))) THEN
                                  mapextf(ifa) = j+s
                                  check = .FALSE.
                                  EXIT
                               ELSEIF((n2 .EQ. Tb_bound(j,1)) .AND. (n1 .EQ. Tb_bound(j,SIZE(Tb_bound,2)))) THEN
                                  mapextf(ifa) = j+s
                                  check = .FALSE.
                                  EXIT
                               ENDIF

                               IF(j .EQ. (SIZE(Tb_bound,1)+1)) THEN
                                  WRITE(*,*) "Strange."
                                  STOP
                               ENDIF
                            ENDDO
                         ENDIF
                         s = s + M(it)%Tb_sizes(ib)
                         DEALLOCATE(Tb_bound)
                      ENDIF
                   ENDDO

                   IF (check) THEN
                      T_CON(counter,:) = Tp(iel,refElPol%Face_nodes(iface,:))
                      mapextf(ifa) = counter+s
                      ghostFaces(counter+s+M(it)%Nintfaces) = 1

                      b1 = 0
                      b2 = 0
                      DO j = 1, SIZE(intFaces_Glob,1)
                         IF ((loc2glob_el(iel) .EQ. intFaces_Glob(j,1)) .AND. (iface .EQ. intFaces_Glob(j,2))) THEN
                            b1 = j
                         ENDIF
                         IF ((loc2glob_el(iel) .EQ. intFaces_Glob(j,3)) .AND. (iface .EQ. intFaces_Glob(j,4))) THEN
                            b2 = j
                         ENDIF
                      ENDDO

                      IF (b1 .NE. 0) THEN
                         loc2glob_fa(clg) = b1
                         IF (b2 .NE. 0) THEN
                            WRITE(*,*) "That is strange. STOP."
                            STOP
                         ENDIF
                         IF (loc2glob_el(iel) .GT. intFaces_Glob(b1,3)) THEN
                            flipFace(counter) = 1
                         ENDIF
                      ELSEIF (b2 .NE. 0) THEN
                         loc2glob_fa(clg) = b2
                         IF (loc2glob_el(iel) .GT. intFaces_Glob(b2,1)) THEN
                            flipFace(counter)=1
                         ENDIF
                      ELSE
                         WRITE(*,*) "That is very strange. STOP."
                         STOP
                      ENDIF

                      counter = counter + 1
                      clg = clg + 1
                   ENDIF
                ENDDO

                ! little check
                IF (ANY(mapextf .EQ. 0)) THEN
                   WRITE(*,*) "Error in mapextf. STOP."
                   STOP
                ENDIF

                ! watch out
                CALL UNIQUE_1D(mapextf, mapextf_unique)
                IF(SIZE(mapextf_unique) .NE. MAXVAL(mapextf)) THEN
                   WRITE(*,*) "Error in mapextf. STOP."
                   STOP
                ENDIF
                DO j = 1, SIZE(mapextf_unique)
                   IF(mapextf_unique(j) .NE. j) THEN
                      WRITE(*,*) "Error in mapextf. STOP."
                      STOP
                   ENDIF
                ENDDO

                ALLOCATE(M(it)%TTb_CON(counter-1,Mesh%Nnodesperface))
                ALLOCATE(M(it)%ghostFlp(counter-1))
                ALLOCATE(M(it)%ghostFaces(M(it)%Nfaces))
                ALLOCATE(M(it)%mapextfa_loc(SIZE(mapextf)))

                M(it)%TTb_CON         = T_CON(1:counter-1,:)
                M(it)%ghostFlp        = flipFace(1:counter-1)
                M(it)%loc2glob_fa     = loc2glob_fa
                M(it)%mapextfa_loc    = mapextf
                M(it)%ghostFaces      = ghostFaces

                DEALLOCATE(T_CON)
                DEALLOCATE(flipFace)
                DEALLOCATE(mapextf)
                DEALLOCATE(ghostFaces)
                DEALLOCATE(loc2glob_fa)
                DEALLOCATE(loc2glob_el)
                DEALLOCATE(mapextf_unique)
                DEALLOCATE(Tp)
                DEALLOCATE(extfaces)

             ENDDO

             !  Find multiple occurrence of the same face across various processes, this part needs to be run un serial
             IF(verbose) THEN
                IF(MPIvar%glob_id .EQ. 0) THEN
                   WRITE(*,*) "Find multiple occurrence of the same face across various processes"
                ENDIF
             ENDIF

             ! Face ownership. An interior face is assembled by the division that
             ! owns the lower-numbered of its two elements, a boundary face by the
             ! division that owns its element. That division always holds every
             ! element adjacent to the face (its own element, plus the neighbour as a
             ! ghost since ghosts are taken from higher-numbered divisions). All the
             ! other occurrences of the face are ghost faces. This replaces the former
             ! sequential "first division that contains the face" rule, which could
             ! leave a face unassembled where three or more divisions meet.
#ifndef PARALL
             DO it = 1, ndiv
#else
                DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif
                   IF(verbose) THEN
                      WRITE(*,*) "Division: ", it
                   ENDIF

                   DO i = 1, M(it)%Nfaces
                      gfa = M(it)%loc2glob_fa(i)
                      IF (gfa .LE. n_int_faces) THEN
                         owner = MIN(part(intfaces_glob(gfa,1)), part(intfaces_glob(gfa,3)))
                      ELSE
                         owner = part(extfaces_glob(gfa-n_int_faces,1))
                      ENDIF
                      ! a face on the border of the local mesh can never be assembled here
                      IF (M(it)%ghostFaces(i) .EQ. 1 .AND. owner .EQ. it) THEN
                         WRITE(*,*) "Domain decomposition: inconsistent ownership of face ", gfa, &
                              " in division ", it, ". STOP."
                         STOP
                      ENDIF
                      IF (owner .NE. it) THEN
                         M(it)%reptFaces(i)  = 1
                         M(it)%ghostFaces(i) = 1
                      ELSE
                         M(it)%reptFaces(i)  = 0
                         M(it)%ghostFaces(i) = 0
                      ENDIF
                   ENDDO

                   ! count the number of ghost faces
                   counter = COUNT(M(it)%ghostFaces .EQ. 1)

                   ! ghostFlp: flips of the ghost faces in local face order. Only the
                   ! faces on the border of the local mesh (stored last) can be flipped.
                   ALLOCATE(auxflp(counter))
                   auxflp = 0
                   s = 1
                   DO i = counter - SIZE(M(it)%ghostFlp) + 1, counter
                      auxflp(i) = M(it)%ghostFlp(s)
                      s = s + 1
                   ENDDO
                   M(it)%ghostFlp = auxflp
                   DEALLOCATE(auxflp)
                ENDDO


#ifdef PARALL
                ! before sending make sure all processes have finished
                CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
                ! Share Loc2Glob_el and Loc2Glob_fa with the other processes.
                DO it = 1, nDiv
                   IF (MPIvar%glob_id+1 .EQ. it) THEN
                      ! If we are the root process, send our data to everyone
                      DO i = 0, MPIvar%glob_size-1
                         IF (i .NE. MPIvar%glob_id) THEN
                            n = SIZE(M(it)%Loc2Glob_el)
                            CALL MPI_Send(n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, ierr)
                            CALL MPI_Send(M(it)%Loc2Glob_el, SIZE(M(it)%Loc2Glob_el), MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
                            n = SIZE(M(it)%Loc2Glob_fa)
                            CALL MPI_Send(n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, ierr)
                            CALL MPI_Send(M(it)%Loc2Glob_fa, SIZE(M(it)%Loc2Glob_fa), MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
                         ENDIF
                      ENDDO
                   ELSE
                      ! If we are a receiver process, receive the data from the root (it-th process)
                      CALL MPI_Recv(n, 1, MPI_INT, it-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                      ALLOCATE(M(it)%Loc2Glob_el(n))
                      CALL MPI_Recv(M(it)%Loc2Glob_el, SIZE(M(it)%Loc2Glob_el), MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                      CALL MPI_Recv(n, 1, MPI_INT, it-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                      ALLOCATE(M(it)%Loc2Glob_fa(n))
                      CALL MPI_Recv(M(it)%Loc2Glob_fa, SIZE(M(it)%Loc2Glob_fa), MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                   ENDIF
                ENDDO
                CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

                ! Create structure for the communications
                IF(verbose) THEN
                   IF(MPIvar%glob_id .EQ. 0) THEN
                      WRITE(*,*) "Create structure for the communications"
                   ENDIF
                ENDIF

#ifndef PARALL
                DO it = 1, nDiv
#else
                   DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif
                      IF(verbose) THEN
                         WRITE(*,*) "Division: ", it
                      ENDIF

                      ! Count number of ghostfaces
                      n = COUNT(M(it)%ghostFaces .EQ. 1)

                      IF(n .EQ. 0) THEN
                         ALLOCATE(M(it)%GhostPro(1))
                         ALLOCATE(M(it)%GhostLoc(1))
                         M(it)%GhostPro = -1
                         M(it)%GhostLoc = -1
                      ELSE
                         ALLOCATE(indices(n))
                         ALLOCATE(indices2(n))
                         ALLOCATE(M(it)%GhostPro(n))
                         ALLOCATE(M(it)%GhostLoc(n))

                         M(it)%GhostPro = 0
                         M(it)%GhostLoc = 0
                         indices = 0
                         indices2 = 0

                         counter = 1
                         DO i = 1, M(it)%Nfaces
                            IF(M(it)%ghostFaces(i) .EQ. 1) THEN
                               indices(counter) = i
                               counter = counter + 1
                            ENDIF
                         ENDDO

                         ! ghost face -> division that assembles it, and its local index there
                         DO j = 1, n
                            gfa = M(it)%Loc2Glob_fa(indices(j))
                            IF (gfa .LE. n_int_faces) THEN
                               owner = MIN(part(intfaces_glob(gfa,1)), part(intfaces_glob(gfa,3)))
                            ELSE
                               owner = part(extfaces_glob(gfa-n_int_faces,1))
                            ENDIF
                            M(it)%GhostPro(j) = owner
                            DO jj = 1, SIZE(M(owner)%Loc2Glob_fa)
                               IF(M(owner)%Loc2Glob_fa(jj) .EQ. gfa) THEN
                                  M(it)%GhostLoc(j) = jj
                                  EXIT
                               ENDIF
                            ENDDO
                            IF (M(it)%GhostLoc(j) .EQ. 0) THEN
                               WRITE(*,*) "Domain decomposition: ghost face ", gfa, " not found in division ", owner, ". STOP."
                               STOP
                            ENDIF
                         ENDDO
                         DEALLOCATE(indices2)
                         DEALLOCATE(indices)
                      ENDIF

                      ! Elements
                      n = COUNT(M(it)%ghostElems .EQ. 1)

                      IF(n .EQ. 0) THEN
                         ALLOCATE(M(it)%GhelsPro(1))
                         ALLOCATE(M(it)%GhelsLoc(1))
                         M(it)%GhelsPro = -1
                         M(it)%GhelsLoc = -1
                      ELSE
                         ALLOCATE(indices(n))
                         ALLOCATE(indices2(n))
                         ALLOCATE(M(it)%GhelsPro(n))
                         ALLOCATE(M(it)%GhelsLoc(n))
                         M(it)%GhelsPro = 0
                         M(it)%GhelsLoc = 0
                         indices = 0
                         indices2 = 0

                         counter = 1
                         DO i = 1, M(it)%Nelems
                            IF(M(it)%ghostElems(i) .EQ. 1) THEN
                               indices(counter) = i
                               counter = counter + 1
                            ENDIF
                         ENDDO

                         ! ghost element -> division that owns it, and its local index there
                         DO j = 1, n
                            iel = M(it)%Loc2Glob_el(indices(j))
                            owner = part(iel)
                            M(it)%GhelsPro(j) = owner
                            DO jj = 1, SIZE(M(owner)%Loc2Glob_el)
                               IF(M(owner)%Loc2Glob_el(jj) .EQ. iel) THEN
                                  M(it)%GhelsLoc(j) = jj
                                  EXIT
                               ENDIF
                            ENDDO
                            IF (M(it)%GhelsLoc(j) .EQ. 0) THEN
                               WRITE(*,*) "Domain decomposition: ghost element ", iel, " not found in division ", owner, ". STOP."
                               STOP
                            ENDIF
                         ENDDO
                         DEALLOCATE(indices2)
                         DEALLOCATE(indices)
                      ENDIF
                   ENDDO

#ifndef PARALL
                   DO it = 1, nDiv
#else
                      DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif

                         ! Merge Tb_CON with Tb
                         ALLOCATE(Tb_prov(SIZE(M(it)%Tb,1),SIZE(M(it)%Tb,2)))
                         ALLOCATE(boundaryFlag_loc(SIZE(M(it)%Tb,1)))

                         Tb_prov = M(it)%Tb
                         boundaryFlag_loc = M(it)%boundaryFlag

                         DEALLOCATE(M(it)%Tb)
                         DEALLOCATE(M(it)%boundaryFlag)

                         ALLOCATE(M(it)%Tb(SIZE(M(it)%TTb_CON,1)+SIZE(Tb_prov,1),SIZE(Tb_prov,2)))
                         ALLOCATE(M(it)%boundaryFlag(SIZE(M(it)%TTb_CON,1)+SIZE(Tb_prov,1)))

                         M(it)%Tb(1:SIZE(Tb_prov,1),:) = Tb_prov
                         M(it)%Tb(SIZE(Tb_prov,1)+1:SIZE(M(it)%Tb,1),:) = M(it)%TTb_CON
                         M(it)%boundaryFlag = 0
                         M(it)%boundaryFlag(1:SIZE(Tb_prov,1)) = boundaryFlag_loc

                         DEALLOCATE(Tb_prov)
                         DEALLOCATE(boundaryFlag_loc)
                      ENDDO


                      ! Rotate back the mesh
#ifndef PARALL
                      DO it = 1, nDiv
#else
                         DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif

                            IF(angle .GT. 1.e-12) THEN
                               M(it)%XX(:,1) = M(it)%XX(:,1) - xcenter
                               M(it)%XX(:,2) = M(it)%XX(:,2) - ycenter
                               RotMat(1,1) = COS((-angle)*pi_loc/180.)
                               RotMat(1,2) = -SIN((-angle)*pi_loc/180.)
                               RotMat(2,1) = SIN((-angle)*pi_loc/180.)
                               RotMat(2,2) = COS((-angle)*pi_loc/180.)
                               M(it)%XX = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(M(it)%XX)))
                               M(it)%XX(:,1) = M(it)%XX(:,1) + xcenter
                               M(it)%XX(:,2) = M(it)%XX(:,2) + ycenter
                            ENDIF

                         ENDDO

#ifdef PARALL
                         CALL send_to_processes
#endif

#ifndef PARALL
                         ! free the structures
                         DO it = 1, nDiv
#else
                            DO it = MPIvar%glob_id+1, MPIvar%glob_id+1
#endif
                               DEALLOCATE(M(it)%TT)
                               DEALLOCATE(M(it)%XX)
                               DEALLOCATE(M(it)%intfaces_loc)
                               DEALLOCATE(M(it)%extfaces_loc)
                               DEALLOCATE(M(it)%TTb_CON)
                               DEALLOCATE(M(it)%Tb)
                               DEALLOCATE(M(it)%reptFaces)
                               DEALLOCATE(M(it)%ghostFaces)
                               DEALLOCATE(M(it)%ghostElems)
                               DEALLOCATE(M(it)%loc2glob_no)
                               DEALLOCATE(M(it)%ghostFlp)
                               DEALLOCATE(M(it)%mapextfa_loc)
                               DEALLOCATE(M(it)%Tb_sizes)
                               DEALLOCATE(M(it)%GhelsPro)
                               DEALLOCATE(M(it)%GhelsLoc)
                               DEALLOCATE(M(it)%GhostPro)
                               DEALLOCATE(M(it)%GhostLoc)
                               DEALLOCATE(M(it)%boundaryFlag)

                            ENDDO

                            ! loc2glob_el and loc2glob_fa are shared by all processes
                            DO it = 1, nDiv
                               DEALLOCATE(M(it)%loc2glob_fa)
                               DEALLOCATE(M(it)%loc2glob_el)
                            ENDDO

                            DEALLOCATE(part)
                            NULLIFY(intfaces_glob, extfaces_glob, T, TTb, boundaryFlag, X)


#ifdef PARALL
                          CONTAINS

                            SUBROUTINE send_to_processes
                              USE printutils

                              IMPLICIT NONE

                              INTEGER                          :: ierr
                              REAL*8                           :: xmin
                              CHARACTER(10)                    :: str
                              REAL*8                           :: tol = 1e-6
#ifdef TOR3D
                              INTEGER                          :: i
#endif
                              CALL free_mesh


                              Mesh%Nelems        = M(MPIvar%glob_id+1)%Nelems
                              Mesh%Nnodesperelem = M(MPIvar%glob_id+1)%Nnodesperelem
                              Mesh%Nnodesperface = M(MPIvar%glob_id+1)%Nnodesperface
                              Mesh%Nnodes        = M(MPIvar%glob_id+1)%Nnodes
                              Mesh%ndim          = M(MPIvar%glob_id+1)%ndim
                              Mesh%Nintfaces     = M(MPIvar%glob_id+1)%Nintfaces
                              Mesh%Nextfaces     = M(MPIvar%glob_id+1)%Nextfaces
                              Mesh%elemType      = M(MPIvar%glob_id+1)%elemType
                              Mesh%Nfaces        = M(MPIvar%glob_id+1)%Nfaces

                              ALLOCATE (Mesh%T(Mesh%Nelems, Mesh%Nnodesperelem))
                              ALLOCATE (Mesh%X(Mesh%Nnodes, Mesh%ndim ))
                              ALLOCATE (Mesh%Tb(Mesh%Nextfaces, Mesh%Nnodesperface))
                              ALLOCATE (Mesh%boundaryFlag(Mesh%Nextfaces))
                              ALLOCATE (Mesh%ghostFaces(Mesh%Nfaces))
                              ALLOCATE (Mesh%loc2glob_fa(Mesh%Nfaces))
                              ALLOCATE (Mesh%loc2glob_el(Mesh%Nelems))
                              ALLOCATE (Mesh%loc2glob_nodes(Mesh%Nnodes))
                              ALLOCATE (Mesh%ghostElems(Mesh%Nelems))
                              !ALLOCATE (Mesh%intFaces(Mesh%Nintfaces,Mesh%Nnodesperface))
                              !ALLOCATE (Mesh%extFaces(Mesh%Nextfaces,Mesh%Nnodesperface))


                              Mesh%T                = M(MPIvar%glob_id+1)%TT
                              Mesh%X                = M(MPIvar%glob_id+1)%XX
                              Mesh%Tb               = M(MPIvar%glob_id+1)%Tb
                              Mesh%boundaryFlag     = M(MPIvar%glob_id+1)%boundaryFlag
                              Mesh%ghostFaces       = M(MPIvar%glob_id+1)%ghostFaces
                              Mesh%ghostElems       = M(MPIvar%glob_id+1)%ghostElems
                              Mesh%loc2glob_fa      = M(MPIvar%glob_id+1)%loc2glob_fa
                              Mesh%loc2glob_el      = M(MPIvar%glob_id+1)%loc2glob_el
                              Mesh%loc2glob_nodes   = M(MPIvar%glob_id+1)%loc2glob_no
                              !Mesh%intFaces        = M(MPIvar%glob_id+1)%intfaces_loc
                              !Mesh%intFaces        = M(MPIvar%glob_id+1)%extfaces_loc

                              ! Find the number of ghost faces
                              Mesh%nghostfaces = SUM(Mesh%ghostFaces)
                              ! Find the number of ghost elements
                              Mesh%nghostElems = SUM(Mesh%ghostElems)

                              IF(Mesh%nghostfaces .EQ. 0) THEN
                                 ALLOCATE (Mesh%ghostflp(1))
                                 ALLOCATE (Mesh%ghostpro(1))
                                 ALLOCATE (Mesh%ghostloc(1))
                              ELSE
                                 ALLOCATE (Mesh%ghostflp(Mesh%nghostfaces))
                                 ALLOCATE (Mesh%ghostpro(Mesh%nghostfaces))
                                 ALLOCATE (Mesh%ghostloc(Mesh%nghostfaces))
                              ENDIF

                              IF(Mesh%nghostElems .EQ. 0) THEN
                                 ALLOCATE (Mesh%ghelspro(1))
                                 ALLOCATE (Mesh%ghelsloc(1))
                              ELSE
                                 ALLOCATE (Mesh%ghelspro(Mesh%nghostElems))
                                 ALLOCATE (Mesh%ghelsloc(Mesh%nghostElems))
                              ENDIF


                              Mesh%ghostflp = M(MPIvar%glob_id+1)%ghostflp
                              Mesh%ghostpro = M(MPIvar%glob_id+1)%ghostpro
                              Mesh%ghostloc = M(MPIvar%glob_id+1)%ghostloc
                              Mesh%ghelspro = M(MPIvar%glob_id+1)%ghelspro
                              Mesh%ghelsloc = M(MPIvar%glob_id+1)%ghelsloc

#ifdef TOR3D
                              IF (MPIvar%ntor > 1) THEN
                                 DO i = 1, SIZE(Mesh%ghostPro)
                                    IF (Mesh%ghostPro(i) .GT. -1) THEN
                                       Mesh%ghostPro(i) = Mesh%ghostPro(i) + (MPIvar%itor - 1)*MPIvar%npol
                                    ENDIF
                                 END DO
                              ENDIF

                              IF (MPIvar%ntor .GT. 1) THEN
                                 DO i = 1, SIZE(Mesh%ghelspro)
                                    IF (Mesh%ghelspro(i) .GT. -1) THEN
                                       Mesh%ghelsPro(i) = Mesh%ghelsPro(i) + (MPIvar%itor - 1)*MPIvar%npol
                                    END IF
                                 END DO
                              END IF

#endif
                              Mesh%Nextfaces_nogho = COUNT((Mesh%boundaryFlag .NE. 0) .AND. (Mesh%ghostFaces(Mesh%Nintfaces+1:SIZE(Mesh%ghostFaces)) .EQ. 0))
                              Mesh%Nintfaces_nogho = COUNT(Mesh%ghostFaces(1:Mesh%Nintfaces) .EQ. 0)
                              CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_el), Mesh%Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_fa), Mesh%Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_nodes), Mesh%Nno_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(Mesh%Nextfaces_nogho, Mesh%Nextfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(Mesh%Nintfaces_nogho, Mesh%Nintfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                              !CALL MPI_ALLREDUCE(Mesh%ndir, Mesh%Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(Mesh%nghostfaces, Mesh%Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

                              xmin = MINVAL(Mesh%X(:,1))

                              CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

                              ! Apply shift if axisymmetric case
                              IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. MINVAL(Mesh%X(:,1)) < tol)) THEN
                                 IF (MPIvar%glob_id .EQ. 0) THEN
                                    WRITE (6, *) "*** Applying translation in axisymmetric case!"
                                 ENDIF
                                 Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
                              END IF

                              ! Apply length scale
                              Mesh%X = Mesh%X/phys%lscale

                              Mesh%xmax = MAXVAL(Mesh%X(:, 1))
                              Mesh%xmin = MINVAL(Mesh%X(:, 1))
                              Mesh%ymax = MAXVAL(Mesh%X(:, 2))
                              Mesh%ymin = MINVAL(Mesh%X(:, 2))

                              CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
                              CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

                              IF (utils%printint > 0) THEN
                                 IF (MPIvar%glob_id .EQ. 0) THEN
                                    IF (Mesh%elemType == 0) THEN
                                       WRITE (str, '(A)') 'triangles'
                                    ELSEIF (Mesh%elemType == 1) THEN
                                       WRITE (str, '(A)') 'quads'
                                    ELSEIF (Mesh%elemType == 2) THEN
                                       WRITE (str, '(A)') 'thetra'
                                    ELSEIF (Mesh%elemType == 3) THEN
                                       WRITE (str, '(A)') 'hexa'
                                    END IF
                                    WRITE (6, *) '*************************************************'
                                    WRITE (6, *) '*                    MESH                       *'
                                    WRITE (6, *) '*************************************************'
                                    WRITE (6, '(A,I18)') ' Number of dimensions:         ', Mesh%ndim
                                    WRITE (6, '(A,A34)') ' Element type: ', TRIM(str)
                                    WRITE (6, '(A,I18)') ' Number of elements:           ', Mesh%Nelems
                                    WRITE (6, '(A,I18)') ' Number of nodes:              ', Mesh%Nnodes
                                    WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Mesh%Nnodesperelem
                                    WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Mesh%Nnodesperface
                                    WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Mesh%Nextfaces
                                    WRITE (6, *) ' '
                                    WRITE (6, *) ' '

                                    IF (utils%printint > 1) THEN
                                       WRITE (6, *) "Connectivity matrix T:"
                                       CALL displayMatrixInt(Mesh%T)
                                       WRITE (6, *) "Boundary connectivity matrix Tb:"
                                       CALL displayMatrixInt(Mesh%Tb)
                                    END IF
                                 ENDIF
                              ENDIF

                            ENDSUBROUTINE send_to_processes
#endif

                          ENDSUBROUTINE split_mesh

  !> Partition the elements into npart parts by graph partitioning of the
  !> element connectivity (dual graph: one vertex per element, one edge per
  !> interior face) with SCOTCH. Computed on rank 0 and broadcast.
  !> part(iel) in 1..npart on success; ok = .FALSE. if SCOTCH is not compiled
  !> in (-DWITH_SCOTCH), fails, or leaves a part empty.
  SUBROUTINE graph_partition(npart, nelems, intfaces, part, ok)
    INTEGER, INTENT(IN)    :: npart, nelems
    INTEGER, INTENT(IN)    :: intfaces(:,:)
    INTEGER, INTENT(OUT)   :: part(nelems)
    LOGICAL, INTENT(OUT)   :: ok
#ifdef WITH_SCOTCH
    INCLUDE 'scotchf.h'
    DOUBLE PRECISION                     :: grafdat(SCOTCH_GRAPHDIM), stradat(SCOTCH_STRATDIM)
    INTEGER(SCOTCH_NUMSIZE), ALLOCATABLE :: verttab(:), edgetab(:), parttab(:)
    INTEGER(SCOTCH_NUMSIZE)              :: baseval, vertnbr, edgenbr, partnbr
    INTEGER, ALLOCATABLE                 :: deg(:), pos(:), cnt(:)
    INTEGER                              :: i, e1, e2, ierr, ncut
    INTEGER                              :: mpierr
#else
    INTEGER                              :: mpierr
#endif

    part = 0
    ok = .FALSE.
#ifdef WITH_SCOTCH
    IF (MPIvar%glob_id .EQ. 0) THEN
       ! dual graph in compressed (CSR) form, base 1, no self loops
       ALLOCATE(deg(nelems), pos(nelems))
       deg = 0
       DO i = 1, SIZE(intfaces,1)
          e1 = intfaces(i,1)
          e2 = intfaces(i,3)
          IF (e1 .EQ. e2) CYCLE
          deg(e1) = deg(e1) + 1
          deg(e2) = deg(e2) + 1
       ENDDO
       ALLOCATE(verttab(nelems+1))
       verttab(1) = 1
       DO i = 1, nelems
          verttab(i+1) = verttab(i) + deg(i)
       ENDDO
       edgenbr = verttab(nelems+1) - 1
       ALLOCATE(edgetab(MAX(edgenbr,1)))
       pos = INT(verttab(1:nelems))
       DO i = 1, SIZE(intfaces,1)
          e1 = intfaces(i,1)
          e2 = intfaces(i,3)
          IF (e1 .EQ. e2) CYCLE
          edgetab(pos(e1)) = e2
          pos(e1) = pos(e1) + 1
          edgetab(pos(e2)) = e1
          pos(e2) = pos(e2) + 1
       ENDDO

       baseval = 1
       vertnbr = nelems
       partnbr = npart
       ALLOCATE(parttab(nelems))
       parttab = -1

       CALL SCOTCHFstratInit(stradat, ierr)
       IF (ierr .EQ. 0) CALL SCOTCHFgraphInit(grafdat, ierr)
       IF (ierr .EQ. 0) THEN
          ! absent optional arrays are passed as verttab (vendtab: compact graph,
          ! velotab, vlbltab) and edgetab (edlotab), see library_graph.c
          CALL SCOTCHFgraphBuild(grafdat, baseval, vertnbr, verttab(1), verttab(2), &
               verttab(1), verttab(1), edgenbr, edgetab(1), edgetab(1), ierr)
          IF (ierr .EQ. 0) CALL SCOTCHFgraphCheck(grafdat, ierr)
          IF (ierr .EQ. 0) CALL SCOTCHFgraphPart(grafdat, partnbr, stradat, parttab(1), ierr)
          CALL SCOTCHFgraphExit(grafdat)
       ENDIF
       CALL SCOTCHFstratExit(stradat)

       IF (ierr .EQ. 0) THEN
          part = INT(parttab) + 1             ! SCOTCH parts are numbered from 0
          ALLOCATE(cnt(npart))
          cnt = 0
          ok = ALL(part .GE. 1 .AND. part .LE. npart)
          IF (ok) THEN
             DO i = 1, nelems
                cnt(part(i)) = cnt(part(i)) + 1
             ENDDO
             ok = ALL(cnt .GT. 0)
          ENDIF
          IF (ok) THEN
             ncut = 0
             DO i = 1, SIZE(intfaces,1)
                IF (part(intfaces(i,1)) .NE. part(intfaces(i,3))) ncut = ncut + 1
             ENDDO
             WRITE(*,'(A,I6,A)') " SCOTCH graph partitioning into ", npart, " parts"
             WRITE(*,'(A,I8,A,I8,A,F8.1)') "   elements per part: min ", MINVAL(cnt), &
                  ", max ", MAXVAL(cnt), ", mean ", REAL(nelems)/REAL(npart)
             WRITE(*,'(A,I8,A,I8,A)') "   interface faces: ", ncut, " (of ", SIZE(intfaces,1), " interior faces)"
          ELSE
             WRITE(*,*) "SCOTCH returned an invalid partition (empty part)."
          ENDIF
          DEALLOCATE(cnt)
       ELSE
          WRITE(*,*) "SCOTCH graph partitioning failed, error code ", ierr
       ENDIF
       DEALLOCATE(deg, pos, verttab, edgetab, parttab)
    ENDIF
#endif

#ifdef PARALL
    CALL MPI_Bcast(ok, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
    IF (ok) CALL MPI_Bcast(part, nelems, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
    IF (.NOT. ok) part = 0
  END SUBROUTINE graph_partition

                        ENDMODULE domain_decomposition_module
