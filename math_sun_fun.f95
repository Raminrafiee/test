! the function & subroutine I use
module math1
contains
subroutine cramers(a,b,x,n)

  implicit none
  integer,intent(in) :: n
  integer :: i
  real *8, dimension(n),intent(in) :: b
  real *8, dimension(n,1),intent(out) :: x
  real *8, dimension(n,n,n+1),intent(inout) :: a
  ! real *8, external :: detf

  do i=1,n
    a(:,:,i+1)= a(:,:,1)
    a(:,i,i+1)=b
    x(i,1)=detf(a(:,:,i+1),n)/detf(a(:,:,1),n)
    print *, "the result is: ",x(i,1)
    end do

end subroutine cramers

recursive real*8 function detf(mat,n) result(det2)

    implicit none
    integer , intent(in) :: n
    integer ::i
    real *8, intent(in), dimension(n,n) :: mat
    real *8,dimension(n-1,n-1) :: slice

  if (n==2) then
      det2 = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
    return
    else if (n==1) then
      det2=mat(1,1)
    return
    else

      do i=1,3
        call clip(mat,slice,n,i)
        det2=det2+(-1.0)**(i+1)*(mat(1,i))*detf(slice,n-1)
      end do
    return
  end if

end function detf


subroutine clip(mat,slice,n,colum)

 implicit none
  integer , intent(in) :: n
  integer :: row,colum
  real *8, intent(in), dimension(n,n) :: mat
  real *8, dimension(n-1,n-1),intent(out):: slice
  logical, dimension(n,n)::mask
  row=1
  mask=.true.
  mask(row,:)=.false.
  mask(:,colum)=.false.

  slice=reshape(pack(mat,mask),(/n-1,n-1/))

end subroutine clip

subroutine print_matrics(a,row,colum)
    implicit none
  integer,intent(in) :: row,colum
  integer :: i,j
  real(kind=8), dimension(row,colum),intent(in) :: a

  5 format(f8.3,t4)
  do i=1,row
    write (*,'("|")',advance='no')
      do j=1,colum
          write (*,fmt=5,advance='no') a(i,j)
      end do
    write (*,'("|")',advance='no')
    write (*,*)
  end do

end subroutine print_matrics

subroutine print_array(a,colum)
  implicit none
  integer,intent(in) :: colum
  integer :: i
  real(kind=8), dimension(colum),intent(in) :: a

  do i=1,colum
      write (*,'("|",f8.3,"|")',advance='no') a(i)
    end do

end subroutine print_array





end module
