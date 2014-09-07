#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
/* #include <mpi.h> */

static hid_t h5_file_id[11];  /* HDF5 file identifier. We allow 10 files to
                                 be opened at once. */
static hid_t h5_grp_id[11];   /* Group identifier */

void open_hdf5_file(int *iun, char *filename)
{
  /* This function opens an HDF5 file collectively so that all the MPI
     processes know how to access it. */
     hid_t  plist_id;   /* Property list */

  /* Set up file access property list with parallel I/O access */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     /*  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);*/

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

void close_hdf5_file(int *iun)
{
  /* Close group id and  HDF5 file */

     H5Fclose(h5_file_id[*iun]);
}

void close_hdf5_file_(int *iun)
{
  /* Fortran interface to the function "close_hdf5_file".
     We add an underscore at the end of the function name. */

     close_hdf5_file(iun);
}

void write_hdf5(int *iun, void *data, char *data_type, int *rank,
                 int *dims, char *name)
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

     for (i=0; i<*rank; i++) {
     dimsf[i]=(hsize_t) dims[i];
     }

     printf("data name = %s\n",name);
     printf("data type = %d    %c\n",*data_type, *data_type);

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
     if ( *rank == 0 ) {
      printf("Creating SCALAR  rank=%d\n\n",*rank);
       aspace_id = H5Screate(H5S_SCALAR);
     }
     else {
       printf("rank=%d  dimsf[0]=%lli  dimsf[1]=%lli\n\n",*rank,dimsf[0],dimsf[1]);
       fflush(stdout);
       aspace_id = H5Screate_simple(*rank, dimsf, NULL);
     }

  /* Create dataset */
     dset_id = H5Dcreate(h5_file_id[*iun], name, type_id, aspace_id,
                         H5P_DEFAULT);

     printf("*** Writing the data ***\n");
       fflush(stdout);

  /* Write dataset to file */ 
     status = H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, 
                       H5P_DEFAULT, data);

     printf("*** Done writing the data ***\n\n");
       fflush(stdout);
  /* Close handles */
     H5Dclose(dset_id);
     H5Sclose(aspace_id);
}

/* ===============================================================
   Fortran interfaces to the C function "write_hdf5()".
   We have 2 possibilities since the Fortran compiler usually
   adds an underscore to the function name although not always. */

void write_hdf5_f(int *iun, void *data, char *data_type, int *rank,
                int *dimsf, char *name)
{
  int i, dims_c[10];

  /* We need to reverse the order of the array dimensions when going from
     Fortran to C.  */
     for (i=0; i<*rank; i++) {
     dims_c[i]=(hsize_t) dimsf[*rank-i-1];
     }

  /* Call C function */
     write_hdf5(iun, data, data_type, rank, dims_c, name);
}

void write_hdf5_f_(int *iun, void *data, char *data_type, int *rank,
                int *dimsf, char *name)
{
  int i, dims_c[10];

  /* We need to reverse the order of the array dimensions when going from
     Fortran to C.  */
     for (i=0; i<*rank; i++) {
     dims_c[i]=(hsize_t) dimsf[*rank-i-1];
     }

  /* Call C function */
     write_hdf5(iun, data, data_type, rank, dims_c, name);
}
