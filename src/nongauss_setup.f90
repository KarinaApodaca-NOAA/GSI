      subroutine diag_1stguess(tsim,nchannl,diag_tsim,diag_tsim_inv)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    diag_1stguess   
!
!   prgmmr: Karina Apodaca (karina.apodaca@noaa.gov)
!
! abstract: 
! abstract: 
!       Calculate the diagonal, inverse, and iverse transpose of the
!       first guess 
!       [Wo^-T=(diag{Hi(xb)})^-T and Wo^1=(diag{Hi(xb)})^-1] for a
!       non-Gaussian DA application aimed at addressing background errors. 
!       This work is based on Fletcher (2017). 
! 
! program history log:
!       input argument list: tsim, nchanl
!       output argument list: diag_tsim_inv, diag_tsim_invt
!   language: f90 and above
!   machine:  
!   
!$$$
!--------
     use kinds, only: r_kind,i_kind
     use constants, only: zero
      implicit none
! Declare local variables
     real(r_kind),intetnt(in) :: tsim(:)
     real(r_kind)             :: diag_tsim(:)
      integer(i_kind)          :: ii,nchannl
     real(r_kind),intent(out) :: diag_tsim_inv(:),diag_tsim_invt(:)
! Get the diagonal, transpose, and inverse 
      diag_tsim=zero
      diag_tsim_inv=zero
      do ii=1,nchannl
            diag_tsim(ii)=tsim(ii,ii)
            diag_tsim_inv(ii)=-1._r_kind/diag_tsim(ii)
            diag_tsim_t(ii)=transpose(diag_tsim(ii))
            diag_tsim_invt=-1._r_kind/diag_tsim_t(ii)
      end do

      end subroutine diag_1stguess
