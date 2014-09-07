
       module mpistuff
           save

           include "mpif.h"

           integer:: comm3d, master
           integer:: status(MPI_STATUS_SIZE), ierr
           integer:: my_id
           integer:: left, right, top, bottom, behind, forward

           integer:: comm1dre, comm1d, comm2dyz

           integer,parameter::ERROR_CARTCOORDS=1
           integer,parameter::ERROR_CARTSHIFT=2
           integer,parameter::ERROR_WAIT=3
           integer,parameter::ERROR_SEND=4
           integer,parameter::ERROR_RECV=5
           integer,parameter::ERROR_ALLREDUCE=6
           integer,parameter::ERROR_SCATTER=7

        end module mpistuff

