        module mesh_common
        save
C                iproc_idx is its process number.

        integer:: iprocx, iprocy, iprocz, iproc_idx
        integer,allocatable:: my_iproc(:)
c        common /procids/ iprocx, iprocy, iprocz, iproc_idx
        end module mesh_common
