#ifdef DYNAMIC
       include 'mesh_uparms.h'
       include 'mesh_parms.h'
#else
       include 'mesh_uparms_static.h'
       include 'mesh_parms_static.h'
#endif
       include 'mesh_common.h'
       include 'mesh.h'
       include 'properties.h'

#ifdef PARALLEL
       include 'mpistuff.h'
#endif
c
        module time
        save
        double precision:: dt, ttot, tstart
        end module time
c
c
       module iounits
        integer, parameter:: istderr=0
        integer, parameter:: istdin=5
        integer, parameter:: istdout=6
        integer, parameter:: ihis=8
        integer, parameter:: imesh=16
        integer, parameter:: ifluid=21
        integer, parameter:: idump=71
        integer, parameter:: ibinoutput=41
        integer, parameter:: itecoutput=61
        integer, parameter:: itime=81
        integer, parameter:: idiag=91
       integer, parameter:: icons=92
       integer, parameter:: ieng=94
       integer, parameter:: iflux=97
       integer, parameter:: itflux=98
       integer, parameter:: ipflux=99
       integer, parameter:: imass=95
       integer, parameter:: ival=88
       integer, parameter:: ivelProf=89
       end module iounits
c
       module OptionFlags
       save
       integer::  tmp
       end module OptionFlags
c
       module RKCoefficients
       save
       double precision::  alpha(3),beta(3),gamma(3),zeta(2)
       end module RKCoefficients
c
       module PressBC
       save
       double precision,allocatable::  pBCx(:,:,:),
     &      pBCy(:,:,:),pBCz(:,:,:)
       end module PressBC
c
       module InflowVelBC
       save
       double precision,allocatable::  Uin(:,:,:,:)
       end module InflowVelBC
c
