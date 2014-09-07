#include <stdio.h>
#include <stdlib.h>
/* #include <mpi.h> */
#include <hdf5.h>   /* mpi.h already included in hdf5.h if this version
                       of the library supports parallel input/output */

static hid_t h5_file_id[11];  /* HDF5 file identifier. We allow 10 files to
                                 be opened at once. */
static hid_t h5_grp_id[11];   /* Group identifier */

void open_parhdf5_file(int *iun, int *communicator, char *filename)
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
     h5_file_id[*iun] = H5Fcreate(filename,
                            H5F_ACC_TRUNC, /* overwrite existing file */
                            H5P_DEFAULT,   /* default file creation */
                            plist_id    /* use mpio property list */
                            );
     H5Pclose(plist_id);

}

void open_parhdf5_file_(int *iun, int *communicator, char *filename)
{
  /* This is just an Fortran interface to the "open_hdf5_file" function. */
     open_parhdf5_file(iun, communicator, filename);
}

void close_hdf5_file(int *iun)
{
  /* Close group id and  HDF5 file */

//     printf("*** Closing file: iunit = %d  h5_file_id = %d\n", *iun, h5_file_id[*iun]);
     fflush(stdout);
     H5Fclose(h5_file_id[*iun]);
}

void close_hdf5_file_(int *iun)
{
  /* Fortran interface to the function "close_hdf5_file".
     We add an underscore at the end of the function name. */

     close_hdf5_file(iun);
}

void write_parhdf5(int *iun, void *data, char *data_type, int rank,
       int *dims, int *countc, int *offsetc, int *communicator, char *name)
{
  /* Write a scalar or array quantity using multiple processors.
   * iun is an integer number representing the output file.
   * data is the actual data to be written
   * d_type is the data type (I=integer, R or F=float, D=double)
   * rank = 0 for a scalar or >= 1 for an array
   * dims[] are the dimensions of the array
   * count[] is the actual number of elements written for each dimension
   * offset[] gives the position in the file where the processor starts writing
   * name is the name of the variable as it will appear in the HDF5 file
   */
     hid_t aspace_id;   /* Dataspace identifier */
     hid_t dset_id;     /* Dataset identifier */
     hid_t plist_id;    /* property list identifier */
     hid_t memspace;
     hid_t type_id;     /* Data type: H5T_NATIVE_{INT,FLOAT,DOUBLE} */
     hsize_t  dimsf[10], count[10], offset[10];
     int    hdferr, i;
     herr_t  status;
     MPI_Comm comm;

     comm = (MPI_Comm)*communicator;

     for (i=0; i<rank; i++) {
       dimsf[i]=(hsize_t) dims[i];
       count[i]=(hsize_t) countc[i];
       offset[i]=(hsize_t) offsetc[i];
     }

     //printf("data name = %s\n",name);
     //printf("data type = %d    %c\n",*data_type, *data_type);

     switch (*data_type) {
     case 'I':
       type_id = H5T_NATIVE_INT;
       break;
     case 'F':
     case 'R':
       type_id = H5T_NATIVE_FLOAT;
       break;
     case 'D':
       type_id = H5T_NATIVE_DOUBLE;
     }

  /* Create data space for the dataset */
     if ( rank == 0 ) {
//      printf("Creating SCALAR  rank=%d  dimsf[0]=%ld\n\n",rank,dimsf[0]);
       aspace_id = H5Screate(H5S_SCALAR);
     }
     else {
//       printf("rank=%d  dimsf[1]=%lli  dimsf[0]=%lli\n\n",rank,dimsf[1],dimsf[0]);
       aspace_id = H5Screate_simple(rank, dimsf, NULL);
     }

  /* Create dataset */
     dset_id = H5Dcreate(h5_file_id[*iun], name, type_id, aspace_id,
                         H5P_DEFAULT);
     H5Sclose(aspace_id);

     //printf("*** Writing the data ***\n");
       fflush(stdout);

    memspace = H5Screate_simple(rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    aspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(aspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    /*    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */
    
    status = H5Dwrite(dset_id, type_id, memspace, aspace_id,
		      plist_id, data);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(aspace_id);
    H5Sclose(memspace);
    H5Pclose(plist_id);


     //printf("*** Done writing the data ***\n\n");
       fflush(stdout);
}

/* Fortran interface to the C function write_parhdf5 */
void write_parhdf5_f(int *iun, void *data, char *data_type, int *rank_f,
       int *dims_f, int *count_f, int *offset_f, int *communicator, char *name)
{
  int dims[10], count[10], offset[10], rank, i;

  rank =  *rank_f;

  /* We need to invert the order of the array dimensions from Fortran
     column-wise to C row-wise. */
  for (i=0; i<rank; i++) {
    dims[i] = dims_f[rank-i-1];
    count[i] = count_f[rank-i-1];
    offset[i] = offset_f[rank-i-1];
  }
  write_parhdf5(iun, data, data_type, rank, dims, count, offset,
                 communicator, name);
}

/* Fortran interface to the C function write_parhdf5 */
void write_parhdf5_f_(int *iun, void *data, char *data_type, int *rank_f,
       int *dims_f, int *count_f, int *offset_f, int *communicator, char *name)
{
  int dims[10], count[10], offset[10], rank, i;

  rank =  *rank_f;

  /* We need to invert the order of the array dimensions from Fortran
     column-wise to C row-wise. */
  for (i=0; i<rank; i++) {
    dims[i] = dims_f[rank-i-1];
    count[i] = count_f[rank-i-1];
    offset[i] = offset_f[rank-i-1];
  }
  write_parhdf5(iun, data, data_type, rank, dims, count, offset,
                communicator, name);
}


void open_hdf5_file(int *iun, char *filename)
{
  /* This function opens an HDF5 file serially using a single processor */
     hid_t  plist_id;   /* Property list */

     plist_id = H5Pcreate(H5P_FILE_ACCESS);

  /* Create a new file collectively and release property list identifier.*/
     h5_file_id[*iun] = H5Fcreate(filename,
                            H5F_ACC_TRUNC, /* overwrite existing file */
                            H5P_DEFAULT,   /* default file creation */
                            plist_id    /* use mpio property list */
                            );
     H5Pclose(plist_id);
}

void open_hdf5_file_(int *iun, char *filename)
{
  /* This is just an Fortran interface to the "open_hdf5_file" function. */
     open_hdf5_file(iun, filename);
}

void write_hdf5(int *iun, void *data, char *data_type, int rank,
                int *dims, int *communicator, char *name)
{
  /* Write an integer scalar or 1D array using a single processor.
   * name is the name of the variable as it will appear in the HDF5 file
   * data is the actual data to be written
   * d_type is the data type (I=integer, R or F=float, D=double)
   * rank = 0 for a scalar or 1 for a 1D array
   * n is the number of values to be written
   */
     hid_t aspace_id;   /* Dataspace identifier */
     hid_t dset_id;     /* Dataset identifier */
     hid_t type_id;     /* Data type: H5T_NATIVE_{INT,FLOAT,DOUBLE} */
     hsize_t  dimsf[10];
     int    hdferr,i;
     herr_t  status;

     int mpi_size, mype;
     MPI_Comm comm;
     MPI_Info info  = MPI_INFO_NULL;

     comm = (MPI_Comm)*communicator;
     MPI_Comm_rank(comm, &mype);

     for (i=0; i<rank; i++) {
     dimsf[i]=(hsize_t) dims[i];
     }

     //printf("data name = %s\n",name);
     //printf("data type = %d    %c\n",*data_type, *data_type);

     switch (*data_type) {
     case 'I':
       type_id = H5T_NATIVE_INT;
       break;
     case 'F':
     case 'R':
       type_id = H5T_NATIVE_FLOAT;
       break;
     case 'D':
       type_id = H5T_NATIVE_DOUBLE;
     }

  /* Create data space for the dataset */
     if ( rank == 0 ) {
      //printf("Creating SCALAR  rank=%d  dimsf[0]=%ld\n\n",rank,dimsf[0]);
       aspace_id = H5Screate(H5S_SCALAR);
     }
     else {
       //printf("rank=%d  dimsf[1]=%lli  dimsf[0]=%lli\n\n",rank,dimsf[1],dimsf[0]);
       fflush(stdout);
       aspace_id = H5Screate_simple(rank, dimsf, NULL);
     }

  /* Create dataset */
     dset_id = H5Dcreate(h5_file_id[*iun], name, type_id, aspace_id,
                         H5P_DEFAULT);

     //printf("*** Writing the data ***\n");
       fflush(stdout);

  /* Write dataset to file */ 
     if (mype == 0)
     status = H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, 
                       H5P_DEFAULT, data);

     //printf("*** Done writing the data ***\n\n");
       fflush(stdout);
  /* Close handles */
     H5Dclose(dset_id);
     H5Sclose(aspace_id);
}

/*
 * Fortran interface to the C function write_hdf5.
 * We need to invert the order of the array dimensions when going from
 * Fortran to C. This needs to be done only for multi-dimensional arrays.
 */
void write_hdf5_f(int *iun, void *data, char *data_type, int *rank_f,
                int *dims_f, int *communicator, char *name)
{
  int dims[10], rank, i;

  rank =  *rank_f;
  for (i=0; i<rank; i++)
    dims[i] = dims_f[rank-i-1];

  write_hdf5(iun, data, data_type, rank, dims, communicator, name);
}

/*
 * Same Fortran interface but with an underscore added to the routine name.
 */
void write_hdf5_f_(int *iun, void *data, char *data_type, int *rank_f,
                int *dims_f, int *communicator, char *name)
{
  int dims[10], rank, i;

  rank =  *rank_f;
  for (i=0; i<rank; i++)
    dims[i] = dims_f[rank-i-1];

  write_hdf5(iun, data, data_type, rank, dims, communicator, name);
}
