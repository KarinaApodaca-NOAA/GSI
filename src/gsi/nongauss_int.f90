      subroutine diag_bstatev(i,j,k,diag_tv,diag_qv,diag_oz,diag_cw,&
                              diag_u,diag_v,diag_sst,diag_qg,&
                              diag_qh,diag_qi,diag_ql,diag_qr,diag_qs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    diag_bstatev   
!
!   prgmmr: Karina Apodaca (karina.apodaca@noaa.gov)
!
! abstract: Calculate the transpose and diagonal of the background
! state vector/guess fields [Wb^T=(diag{Xb})^T] for a non-Gaussian DA
! application aimed at addressing background errors. This work is based on:
! Fletcher, S. J.(2017). Non-Gaussian Variational Data Assimilation. In
! S.J. Fletcher (Ed.). Data Assimilation for the Geosciences (pp. 823-868).
! Netherlands, United Kingdom, United States: Elsevier. 
! Refered to as Fletcher (2017), hereafter. 
!
! program history log:
!
!   input argument list:
!   Guess fields:
!   tvges_itisig, qvges_itsig, ozges_itsigp,  ges_cwmr_itsig,
!   uges_itsig, vges_itsig, sstges_itsig, qgges_itsig, qhges_itsig
!   qiges_itsig, qlges_itsig, qrges_itsig, qhges_itsig
!
!   output argument list:
!   Diagonal of the traspossed guess fields:
!   diag_tv, diag_qv, diag_oz, diag_cw, diag_u, diag_v, diag_sst,
!   diag_qg, diag_qh, diag_qi, diag_ql, diag_qr, diag_qs
!
!   language: f90 and above
!   machine:  
!   
!$$$
!--------
      use kinds, only: r_kind,i_kind
      use constants, only: zero
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_metguess_mod,  only: gsi_metguess_bundle
      implicit none
! Declare input variables
     real(r_kind),pointer,dimension(:,:  )::psges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:  )::psges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::uges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::uges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::vges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::vges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::qges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::qges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::ozges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::ozges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::tgasges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::tgasges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::aeroges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::aeroges_itsigp=>NULL()
     real(r_kind),pointer,dimension(:,:,:)::qgges_itsig =>NULL()
     real(r_kind),pointer,dimension(:,:,:)::qhges_itsig =>NULL()!KA 
     real(r_kind),pointer,dimension(:,:,:)::qiges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::qlges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::qrges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::qsges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::qvges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::tvges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:,:)::cwges_itsig =>NULL()!KA  
     real(r_kind),pointer,dimension(:,:  )::sst_itsig =>NULL()!KA 

! Declare local variables
      integer(i_kind) :: i,j,k,ierr,istatus,itsig,itsigp
      integer(i_kind)::ierr,istatus,imax,jmax
      integer(i_kind)::kmax_q
      real(r_kind) :: tvges_t,qvges_t,ozges_t,cwges_t,uges_t,vges_t,sstges_t
      real(r_kind) :: qgges_t,qhges_t,qiges_t,qlges_t,qrges_t,qsges_t
      logical :: ngauss_solver_l,ngauss_solver_m

! Declare passed variables
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::uges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::vges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::ozges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qgges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qhges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qiges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qlges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qrges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qsges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::qvges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::tvges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(in)::cwges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1),intent(in)::sstges_itsig
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(out)::diag_tv,diag_qv,diag_oz
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(out)::diag_cw,diag_u,diag_v
     real(r_kind),dimension,(1:imax,1:jmax,1),intent(out)::diag_sst
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(out)::diag_qi,diag_ql,diag_qr
     real(r_kind),dimension(1:imax,1:jmax,1:kmax_q),intent(out)::diag_qs,diag_qg,diag_qh

! Read background fields
      if (ngauss_solver_l .or. ngauss_solver_m) then
      
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'ps',psges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'ps',psges_itsigp,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'u' ,uges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'u' ,uges_itsigp,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'v' ,vges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'v' ,vges_itsigp,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'oz',ozges_itsig,iozs)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'oz',ozges_itsigp,iozs)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'q',qges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'q',qges_itsigp,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'sst',sstges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'tv',tvges_itisig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'cw',ges_cwmr_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qi',qiges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'ql',qlges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qv',qvges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qr',qrges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qg',qgges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qh',qhges_itsig,istatus)
     call gsi_bundlegetpointer(gsi_metguess_bundle(itsig),'qs',qsges_itsig,istatus)

! Initialize arrays
         sstges_t=zero; sstges_itsig=zero; tvges_t=zero; tvges_itisig=zero; 
         qvges_t=zero; qvges_itisig=zero; ozges_t=zero; ozges_itisig=zero
         cwges_t=zero; uges_t=zero; vges_t=zero; qgges_t=zero; qgges_itisig=zero; 
         qhges_t=zero; qhges_itisig=zero; qiges_t=zero; qiges_itisig=zero;
         qlges_t=zero; qrges_itisig=zero; qrges_t=zero; qrges_itisig=zero; 
         qsges_t=zero; qsges_itisig=zero

! Transpose background/guess fields
! 2D-fields
      do i=1,imax
         do j=1,jmax
            sstges_t=sstges_itsig(i,j)
            sstges_itsig(i,j)=sstges_itsig(j,i)
            sstges_itsig(j,i)=sstges_t
         end do
      end do
! 3D-fields
      do k=1,kmax_q
         do i=1,imax
            do j=1,jmax
              tvges_t=tvges_itisig(i,j,k)
              tvges_itisig(i,j,k)=tvges_itisig(j,i,k)
              tvges_itisig(j,i,k)=tvges_t
              qvges_t=qvges_itisig(i,j,k)
              qvges_itisig(i,j,k)=qvges_itisig(j,i,k)
              qvges_itisig(j,i,k)=qvges_t
              ozges_t=ozges_itisig(i,j,k)
              ozges_itisig(i,j,k)=ozges_itisig(j,i,k)
              ozges_itisig(j,i,k)=ozges_t
              cwges_t=cwges_itisig(i,j,k)
              cwges_itisig(i,j,k)=cwges_itisig(j,i,k)
              cwges_itisig(j,i,k)=cwges_t
              uges_t=uges_itisig(i,j,k)
              uges_itisig(i,j,k)=uges_itisig(j,i,k)
              uges_itisig(j,i,k)=uges_t
              vges_t=vges_itisig(i,j,k)
              vges_itisig(i,j,k)=vges_itisig(j,i,k)
              vges_itisig(j,i,k)=vges_t
              qgges_t=qgges_itisig(i,j,k)
              qgges_itisig(i,j,k)=qgges_itisig(j,i,k)
              qgges_itisig(j,i,k)=qgges_t
              qhges_t=qhges_itisig(i,j,k)
              qhges_itisig(i,j,k)=qhges_itisig(j,i,k)
              qhges_itisig(j,i,k)=qhges_t
              qiges_t=qiges_itisig(i,j,k)
              qiges_itisig(i,j,k)=qiges_itisig(j,i,k)
              qiges_itisig(j,i,k)=qiges_t
              qlges_t=qrges_itisig(i,j,k)
              qsges_itisig(i,j,k)=qsges_itisig(j,i,k)
              qsges_itisig(j,i,k)=qsges_t
            end do
         end do
      end do

! Get the diagonal of the transpossed fields
! 2D-fields
! 2D-fields
      do i=1,imax
         do j=1,jmax
            if (i.eq.j) then
                diag_sst(i,j)=sstges_itsig(j,i)
            else
                diag_sst(i,j)=zero
            end if
         end do
      end do
! 3D-fields   
      do k=1,kmax_q
         do i=1,imax
            do j=1,jmax
               if (i.eq.j) then
                   diag_tv(i,j,k)=tvges_itsig(j,i,k)
                   diag_qv(i,j,k)=qvges_itsig(j,i,k)
                   diag_oz(i,j,k)=ozges_itsig(j,i,k)
                   diag_cw(i,j,k)=cwges_itsig(j,i,k)
                   diag_u(i,j,k)=uges_itsig(j,i,k)
                   diag_v(i,j,k)=vges_itsig(j,i,k)
                   diag_qg(i,j,k)=qgges_itsig(j,i,k)
                   diag_qh(i,j,k)=qhges_itsig(j,i,k)
                   diag_qi(i,j,k)=qiges_itsig(j,i,k)
                   diag_ql(i,j,k)=qlges_itsig(j,i,k)
                   diag_qr(i,j,k)=qrges_itsig(j,i,k)
                   diag_qs(i,j,k)=qsges_itsig(j,i,k)
               else
                   diag_tv(i,j,k)=zero
                   diag_qv(i,j,k)=zero
                   diag_cw(i,j,k)=zero
                   diag_u(i,j,k)=zero
                   diag_v(i,j,k)=zero
                   diag_qg(i,j,k)=zero
                   diag_qh(i,j,k)=zero
                   diag_qi(i,j,k)=zero
                   diag_ql(i,j,k)=zero
                   diag_qr(i,j,k)=zero
                   diag_qs(i,j,k)=zero
                   diag_tv(i,j,k)=zero
               end if
            end do
         end do
      end do
      end if
      end subroutine diag_bstatev
                                                                                                                             228,1         Bot
