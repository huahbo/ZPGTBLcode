#include <stdio.h>
#include <stdlib.h>
/* #include <mpi.h> */
#include <hdf5.h>   /* mpi.h already included in hdf5.h if this version
                       of the library supports parallel input/output */

static hid_t h5_file_id[11];  /* HDF5 file identifier. We allow 10 files to
                                 be opened at once. */
static hid_t h5_grp_id[11];   /* Group identifier */

void open_parhdf5_file_rdonly(int *iun, int *communicator, char *filename)
{
     MPI_Comm  comm;
  /* This function opens an HDF5 file collectively so that all the MPI
     processes know how to access it. */
     hid_t  plist_id;   /* Property list */

     comm = (MPI_Comm)*communicator;
  /* Set up file access property list with parallel I/O access */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

  /* Create a new file collectively and release property list identifier.*/
     h5_file_id[*iun] = H5Fopen(filename,
                            H5F_ACC_RDONLY, /* read only */
                            plist_id    /* use mpio property list */
                            );
     H5Pclose(plist_id);

}

void open_parhdf5_file_rdonly_(int *iun, int *communicator, char *filename)
{
  /* This is just an Fortran interface to the "open_hdf5_file" function. */
     open_parhdf5_file_rdonly(iun, communicator, filename);
}

void close_hdf5_file_rdonly(int *iun)
{
  /* Close group id and  HDF5 file */

     printf("*** Closing file: iunit = %d  h5_file_id = %d\n", *iun, h5_file_id[*iun]);
     fflush(stdout);
     H5Fclose(h5_file_id[*iun]);
}

void close_hdf5_file_rdonly_(int *iun)
{
  /* Fortran interface to the function "close_hdf5_file_rdonly".
     We add an underscore at the end of the function name. */

     close_hdf5_file_rdonly(iun);
}

/* Read dataset from hdf5 files:
 * data: pointer to the memory with data stored,
 * data_type: data type  H5T_NATIVE_{INT, FLOAT, DOUBLE},
 * rank: 0 is scalar and >=1 is array
 * dims, countc, and offsetc are all given as parameters
 * data memory must have been allocated outside this function
 */
void read_parhdf5(int *iun, void *data, int *data_type, int *rank,
       int *dims, int *countc, int *offsetc, int *communicator, char *name)
{
     hid_t aspace_id;   /* Dataspace identifier */
     hid_t dset_id;     /* Dataset identifier */
     hid_t plist_id;    /* property list identifier */
     hid_t memspace;
     hid_t filespace;
     hid_t type_id;     /* Data type: H5T_NATIVE_{INT,FLOAT,DOUBLE} */
     hsize_t  dimsf[10], count[10], offset[10];
     int    hdferr, i;
     herr_t  status;
     MPI_Comm comm;

     comm = (MPI_Comm)*communicator;

  /* Open dataset */
     dset_id = H5Dopen(h5_file_id[*iun], name);
                         

     printf("*** Read the data ***\n");
       fflush(stdout);

  /* Get type of the data*/
     type_id = H5Dget_type(dset_id);

     printf("Data type = ");
     if(H5Tequal(type_id, H5T_NATIVE_INT)){
       printf("INT \n");
       *data_type = 0;
       fflush(stdout);
     } else if (H5Tequal(type_id, H5T_NATIVE_FLOAT)){
       printf("FLOAT \n");
       *data_type = 1;
       fflush(stdout);
     } else if (H5Tequal(type_id, H5T_NATIVE_DOUBLE)){
       printf("DOUBLE \n");
       *data_type = 2;
       fflush(stdout);
     } else {
       printf("UNRECOGNIZED \n Aborting...\n");
       fflush(stdout);
       H5Dclose(dset_id);
       return;
     }

  /* Get filespace */
     filespace = H5Dget_space(dset_id);
  
  /* Get dimensions */
     *rank = H5Sget_simple_extent_ndims(filespace);
/*   status = H5Sget_simple_extent_dims(filespace, dims, NULL); */

     printf("Rank: %d\n", *rank);
/*     printf("Dimensions: ");
     for(i=0;i<rank;i++){
       printf("%d ",dims[i]);
     }
     printf("\n"); */
        fflush(stdout);

     for (i=0; i<*rank; i++) {
       dimsf[i]=(hsize_t) dims[i];
       count[i]=(hsize_t) countc[i];
       offset[i]=(hsize_t) offsetc[i];
     }

     if ( rank == 0 ) {
       memspace = H5Screate(H5S_SCALAR);
     }
     else {
       memspace = H5Screate_simple(*rank, count, NULL);
     }

    /*
     * Select hyperslab in the file.
     */
    aspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(aspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Create property list for collective dataset read.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    /*    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */

  /* Read dataset from file */ 
     status = H5Dread(dset_id, type_id, memspace, aspace_id,
		      plist_id, data);

     printf("*** Done reading the data ***\n\n");
       fflush(stdout);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(aspace_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

}

/* Fortran interface to the C function write_parhdf5 */
void read_parhdf5_f(int *iun, void *data, int *data_type, int *rank_f,
       int *dims_f, int *count_f, int *offset_f, int *communicator, char *name)
{
  int dims[10], count[10], offset[10], i;

  int *rank =  rank_f;

  /* We need to invert the order of the array dimensions from Fortran
     column-wise to C row-wise. */
  for (i=0; i<*rank; i++) {
    dims[i] = dims_f[*rank-i-1];
    count[i] = count_f[*rank-i-1];
    offset[i] = offset_f[*rank-i-1];
  }
  read_parhdf5(iun, data, data_type, rank, dims, count, offset,
                 communicator, name);
}

/* Fortran interface to the C function write_parhdf5 */
void read_parhdf5_f_(int *iun, void *data, int *data_type, int *rank_f,
       int *dims_f, int *count_f, int *offset_f, int *communicator, char *name)
{
  int dims[10], count[10], offset[10], i;
  
  int *rank =  rank_f;

  /* We need to invert the order of the array dimensions from Fortran
     column-wise to C row-wise. */
  for (i=0; i<*rank; i++) {
    dims[i] = dims_f[*rank-i-1];
    count[i] = count_f[*rank-i-1];
    offset[i] = offset_f[*rank-i-1];
  }
  read_parhdf5(iun, data, data_type, rank, dims, count, offset,
                communicator, name);
}


void open_hdf5_file_rdonly(int *iun, char *filename)
{
  /* This function opens an HDF5 file serially using a single processor */
     hid_t  plist_id;   /* Property list */

     plist_id = H5Pcreate(H5P_FILE_ACCESS);

  /* Create a new file collectively and release property list identifier.*/
     h5_file_id[*iun] = H5Fopen(filename,
                            H5F_ACC_RDONLY, /* read only */
                            plist_id    /* use mpio property list */
                            );
     H5Pclose(plist_id);
}

void open_hdf5_file_rdonly_(int *iun, char *filename)
{
  /* This is just an Fortran interface to the "open_hdf5_file" function. */
     open_hdf5_file_rdonly(iun, filename);
}

/* Read dataset from hdf5 files:
 * data: pointer to the memory with data stored,
 * data_type: data type  H5T_NATIVE_{INT, FLOAT, DOUBLE},
 * rank: 0 is scalar or 1 is 1D array
 * dims are given as parameters
 * data memory must have been allocated outside this function
 */
void read_hdf5(int *iun, void *data, int *data_type, int *rank,
                int *dims, int *communicator, char *name)
{
     hid_t aspace_id;   /* Dataspace identifier */
     hid_t dset_id;     /* Dataset identifier */
     hid_t type_id;     /* Data type: H5T_NATIVE_{INT,FLOAT,DOUBLE} */
     hid_t filespace;
     hsize_t  dimsf[10];
     int    hdferr,i;
     herr_t  status;

     int mpi_size, mype;
     MPI_Comm comm;
     MPI_Info info  = MPI_INFO_NULL;

     comm = (MPI_Comm)*communicator;
     MPI_Comm_rank(comm, &mype);


  /* Open dataset */
     dset_id = H5Dopen(h5_file_id[*iun], name);

     printf("*** Read the data ***\n");
       fflush(stdout);

  /* Get type of the data*/
     type_id = H5Dget_type(dset_id);

     printf("Data type = ");
     if(H5Tequal(type_id, H5T_NATIVE_INT)){
       printf("INT \n");
       *data_type = 0;
       fflush(stdout);
     } else if (H5Tequal(type_id, H5T_NATIVE_FLOAT)){
       printf("FLOAT \n");
       *data_type = 1;
       fflush(stdout);
     } else if (H5Tequal(type_id, H5T_NATIVE_DOUBLE)){
       printf("DOUBLE \n");
       *data_type = 2;
       fflush(stdout);
     } else {
       printf("UNRECOGNIZED \n Aborting...\n");
       fflush(stdout);
       H5Dclose(dset_id);
       return;
     }

  /* Get filespace */
     filespace = H5Dget_space(dset_id);
  
  /* Get dimensions */
     *rank = H5Sget_simple_extent_ndims(filespace);
/*   status = H5Sget_simple_extent_dims(filespace, dims, NULL); */

     printf("Rank: %d\n", *rank);
/*     printf("Dimensions: ");
     for(i=0;i<*rank;i++){
       printf("%d ",dims[i]);
     }
     printf("\n"); */
        fflush(stdout);

     for (i=0; i<*rank; i++) {
     dimsf[i]=(hsize_t) dims[i];
     }

  /* Read dataset from file */ 
     status = H5Dread(dset_id, type_id, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, data);

     printf("*** Done reading the data ***\n\n");
       fflush(stdout);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
}

/*
 * Fortran interface to the C function write_hdf5.
 * We need to invert the order of the array dimensions when going from
 * Fortran to C. This needs to be done only for multi-dimensional arrays.
 */
void read_hdf5_f(int *iun, void *data, int *data_type, int *rank_f,
                int *dims_f, int *communicator, char *name)
{
  int dims[10], i;

  int *rank =  rank_f;
  for (i=0; i<*rank; i++)
    dims[i] = dims_f[*rank-i-1];

  read_hdf5(iun, data, data_type, rank, dims, communicator, name);
}

/*
 * Same Fortran interface but with an underscore added to the routine name.
 */
void read_hdf5_f_(int *iun, void *data, int *data_type, int *rank_f,
                int *dims_f, int *communicator, char *name)
{
  int dims[10], i;

  int *rank =  rank_f;
  for (i=0; i<*rank; i++)
    dims[i] = dims_f[*rank-i-1];

  read_hdf5(iun, data, data_type, rank, dims, communicator, name);
}
