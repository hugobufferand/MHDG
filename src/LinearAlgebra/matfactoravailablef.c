/*
 * Fortran-callable wrapper of MatGetFactorAvailable() for PETSc installations
 * that do not provide its Fortran binding.
 *
 * Fortran usage (implicit interface, see solve_petsc.f90):
 *
 *   EXTERNAL :: matgetfactoravailable_wrapper
 *   INTEGER  :: avail
 *   CALL matgetfactoravailable_wrapper(A, 'mumps', MAT_FACTOR_LU, avail, ierr)
 *
 *   A      Mat            matrix (Fortran handle, passed by reference)
 *   solver CHARACTER(*)   factorization package; trailing blanks are ignored
 *   ftype  MatFactorType  MAT_FACTOR_LU, MAT_FACTOR_CHOLESKY, ...
 *   avail  INTEGER        1 if the factorization is available, 0 otherwise
 *   ierr   PetscErrorCode
 *
 * The flag is a plain int rather than a PetscBool/LOGICAL, whose size and
 * value of .TRUE. depend on the Fortran compiler. 'len' is the hidden length
 * of the CHARACTER argument, appended by the Fortran compiler.
 */
#include <petscmat.h>

#define SOLVER_NAME_MAX 256

void matgetfactoravailable_wrapper_(Mat *mat, const char *solver, int *ftype,
                                    int *avail, PetscErrorCode *ierr,
                                    PETSC_FORTRAN_CHARLEN_T len)
{
  char      name[SOLVER_NAME_MAX];
  size_t    n = (size_t)len;
  size_t    i;
  PetscBool flg = PETSC_FALSE;

  while (n > 0 && solver[n - 1] == ' ') n--;          /* strip Fortran blank padding */
  if (n >= SOLVER_NAME_MAX) n = SOLVER_NAME_MAX - 1;
  for (i = 0; i < n; i++) name[i] = solver[i];
  name[n] = '\0';

  *ierr  = MatGetFactorAvailable(*mat, name, (MatFactorType)(*ftype), &flg);
  *avail = (flg == PETSC_TRUE) ? 1 : 0;
}
