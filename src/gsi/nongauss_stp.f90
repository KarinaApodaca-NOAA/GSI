      subroutine diag_cvincr(ii,diag_stv_invt,diag_sq_invt,&
                             diag_scw_invt,diag_soz_invt,diag_sst_invt,&
                             diag_su_invt,diag_sv_invt,diag_sqg_invt,&
                             diag_sqh_invt,diag_sqi_invt,diag_sql_invt,&
                             diag_sqr_invt,diag_sqs_invt)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    diag_cvincr   
!
!   prgmmr: Karina Apodaca (karina.apodaca@noaa.gov)
!
! abstract: 
!       Calculate the diagonal, inverse, and iverse transpose of the 
!       control variable increment [delta_Wb=(diag{delta_Xi})^-T] 
!       for a non-Gaussian DA application aimed at addressing background
!       errors. 
!       This work is based on:
!       Fletcher, S. J.(2017). Non-Gaussian Variational Data Assimilation. In
!       S.J. Fletcher (Ed.). Data Assimilation for the Geosciences (pp.
!       823-868). Netherlands, United Kingdom, United States: Elsevier. 
!       Refered to as Fletcher (2017), hereafter. 
! 
! program history log:
!
!       input argument list:
!       stv,sq,scw,soz,su,sv,sqg,sqh,sqi,sql,sqr,sqs
!
!       output argument list:
!       diag_stv_invt,diag_sq_invt,diag_scw_invt,diag_soz_invt
!       diag_su_invt,diag_sv_invt,diag_sqg_invt,diag_sqh_invt
!       diag_sqi_invt,diag_sql_invt,diag_sqr_invt,diag_sqs_invt
!
!   language: f90 and above
!   machine:  
!   
!$$$
!--------
      use kinds, only: r_kind,i_kind
      use constants, only: zero
      use gridmod, only: nsig
      use gsi_4dvar, only: nobs_bins
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      implicit none
! Declare local variables
      integer(i_kind) ii
      real(r_kind),pointer,dimension(:) :: st,sq,scw,soz,su,sv,sqg
      real(r_kind),pointer,dimension(:) :: sqh,sqi,sql,sqr,sqs
      real(r_kind),pointer,dimension(:) :: sst
      real(r_kind),intent(in) :: stv(:),sq(:),scw(:),soz(:),sst(:),su(:)
      real(r_kind),intent(in) :: sv(:),sqg(:),sqh(:),sqi(:)
      real(r_kind),intent(in) :: sql(:),sqr(:),sqs(:)
      real(r_kind),intent(out):: diag_tv_invt(:),diag_sq_invt(:)
      real(r_kind),intent(out):: diag_scw_invt(:),diag_soz_invt(:)
      real(r_kind),intent(out):: diag_sst_invt(:),diag_su_invt(:)
      real(r_kind),intent(out):: diag_sv_invt(:),diag_sqg_invt(:)
      real(r_kind),intent(out):: diag_sqh_invt(:),diag_sqi_invt(:)
      real(r_kind),intent(out):: diag_sql_invt(:),diag_sqr_invt(:)
      real(r_kind)            :: diag_stv(:),diag_sqs_invt(:)
      real(r_kind)            :: diag_sq(:),diag_scw(:),diag_soz(:)
      real(r_kind)            :: diag_sst(:),diag_su(:),diag_sv(:)
      real(r_kind)            :: diag_sqg(:),diag_sqh(:),diag_sqi(:)
      real(r_kind)            :: diag_sql(:),diag_sqr(:),diag_sqs(:)
      type(gsi_bundle), intent(in   ) :: sval

! Retrieve pointers; return when not found (except in case of
! non-essentials)
      if (nong_solver_l .or. nong_solver_m) then

         call gsi_bundlegetpointer(sval,'u',  su, istatus)
         call gsi_bundlegetpointer(sval,'v',  sv, istatus)
         call gsi_bundlegetpointer(sval,'tv' ,st, istatus)
         call gsi_bundlegetpointer(sval,'q',  sq, istatus)
         call gsi_bundlegetpointer(sval,'cw' ,scw,istatus)
         call gsi_bundlegetpointer(sval,'oz' ,soz,istatus)
         call gsi_bundlegetpointer(sval,'sst',sst,istatus)
         call gsi_bundlegetpointer(sval,'qg' ,sqg,istatus)
         call gsi_bundlegetpointer(sval,'qh' ,sqh,istatus)
         call gsi_bundlegetpointer(sval,'qi' ,sqi,istatus)
         call gsi_bundlegetpointer(sval,'ql' ,sql,istatus)
         call gsi_bundlegetpointer(sval,'qr' ,sqr,istatus)
         call gsi_bundlegetpointer(sval,'qs' ,sqs,istatus)

! Initialize arrays
         diag_stv_invt=zero; diag_sq_invt=zero;  diag_scw_invt=zero;
         diag_soz_invt=zero; diag_sst_invt=zero; diag_su_invt=zero;
         diag_sv_invt=zero;  diag_sqg_invt=zero; diag_sqh_invt=zero;
         diag_sqi_invt=zero; diag_sql_invt=zero; diag_sqr_invt=zero;
         diag_sqs_invt=zero
         diag_tv=zero; diag_sq=zero; diag_scw=zero; diag_soz=zero;
         diag_sst=zero; diag_su=zero; diag_sv=zero;
         diag_sqg=zero; diag_sqh=zero; diag_sqi=zero;
         diag_sql=zero; diag_sqr=zero; diag_sqs=zero;

! Get the diagonal, transpose, and/or inverse 

      do ii=i,nobs_bins

         diag_stv(ii)=stv(ii,ii)
         diag_stv_invt(ii)=1._r_kind/(transpose(diag_stv(ii)))

         diag_sq(ii)=sq(ii,ii)
         diag_sq_invt(ii)=1._r_kind/(transpose(diag_sq(ii)))

         diag_scw(ii)=scw(ii,ii)
         diag_scw_invt(ii)=1._r_kind/(transpose(diag_scw(ii)))

         diag_soz(ii)=soz(ii,ii)
         diag_soz_invt(ii)=1._r_kind/(transpose(diag_soz(ii)))

         diag_sst(ii)=sst(ii,ii)
         diag_sst_invt(ii)=1._r_kind/(transpose(diag_sst(ii)))

         diag_su(ii)=su(ii,ii)
         diag_su_invt(ii)=1._r_kind/(transpose(diag_su(ii)))

         diag_sv(ii)=sv(ii,ii)
         diag_sv_invt(ii)=1._r_kind/(transpose(diag_sv(ii)))

         diag_sqg(ii)=sqg(ii,ii)
         diag_sqg_invt(ii)=1._r_kind/(transpose(diag_sqg(ii)))

         diag_sqh(ii)=sqh(ii,ii)
         diag_sqh_invt(ii)=1._r_kind/(transpose(diag_sqh(ii)))

         diag_sqi(ii)=s(ii,ii)
         diag_sqi_invt(ii)=1._r_kind/(transpose(diag_sqi(ii)))

         diag_sql(ii)=sql(ii,ii)
         diag_sql_invt(ii)=1._r_kind/(transpose(diag_sql(ii)))

         diag_sqr(ii)=sqr(ii,ii)
         diag_sqr_invt(ii)=1._r_kind/(transpose(diag_sqr(ii)))

         diag_sqs(ii)=sqs(ii,ii)
         diag_sqs_invt(ii)=1._r_kind/(transpose(diag_sqs(ii)))

      end do

      end if

      end subroutine diag_cvincr
                                                                 
