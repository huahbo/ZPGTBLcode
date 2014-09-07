
      subroutine SetupDomain
      use mesh
      use mesh_common
      implicit none
      call RectangularShockTube
      return
      end
c
c----------------------------------------------------------------------
      subroutine RectangularShockTube
      use mesh
      use mesh_common
      implicit none
        integer:: i,j, k
       integer:: ig, jg, kg
c
c        ilo=1
c        ihi=Nxlsize
c
c     Uncomment the following 4 lines only for the sinus case.
c        xl=-1.0
c        xr=8.0*atan(1.0)+1.0
c        yl=xl-3
c        yr=xr-3
       dx=(xbound(2)-xbound(1))/(nx)
       dy=(ybound(2)-ybound(1))/(ny)
       dz=(zbound(2)-zbound(1))/nz
       dxi=1.D0/dx
       dyi=1.D0/dy
       dzi=1.D0/dz
       dxi3=dxi/3.D0
       dyi3=dyi/3.D0
       dzi3=dzi/3.D0
c       dy=dx
        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
           write(8,*) 'Index Bounds',ixlo,ixhi,iylo,iyhi,izlo,izhi
           write(8,*) 'xl=',xbound(1), 'xr=',xbound(2)
           write(8,*) 'yl=',ybound(1), 'yr=',ybound(2)
           write(8,*) 'zl=',ybound(1), 'zr=',zbound(2)
           write(8,*) 'dx=', dx
           write(8,*) 'dy=', dy
           write(8,*) 'dz=', dz
        endif
c       Get the global indices
        call iglobal(iprocx,1,ig)
        call jglobal(iprocy,1,jg)
        call kglobal(iprocz,1,kg)
c
      do k=IZLO, IZHI,1
      do j=IYLO, IYHI,1
      do i=IXLO, IXHI,1
         xc(i)=xbound(1)+(ig+i-1)*dx-0.5*dx
         yc(j)=ybound(1)+(jg+j-1)*dy-0.5*dy
         zc(k)=zbound(1)+(kg+k-1)*dz-0.5*dz
      enddo
      enddo
      enddo
c
      return
      end
c     

c     







