#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <hdf5.h>


#define ITOT 1024
#define JTOT 1024

#define DELTAX (1./ITOT)
#define DELTAY (1./JTOT)

#define ALLOCATE(ptr, size){ \
  if(((ptr) = calloc(1, (size))) == NULL){ \
    printf("memory allocation error: %s\n", #ptr); \
    exit(EXIT_FAILURE); \
  } \
}

#define DEALLOCATE(ptr){ \
  free((ptr)); \
  (ptr) = NULL; \
}

#define FOPEN(fp, fname, mode){ \
  if(((fp) = fopen((fname), (mode))) == NULL) { \
    printf("file open error %s\n", (fname)); \
    exit(EXIT_FAILURE); \
  } \
}

#define FCLOSE(fp) { \
  fclose(fp); \
  fp = NULL; \
}

static int read(const hid_t file_id, const char dataset_name[], const hid_t type_id, void *data){
  // Read the dataset "dataset_name" from the given HDF5 file "file_id"
  hid_t dataspace_id, dataset_id;
  const hsize_t dims[2] = {JTOT, ITOT};
  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  H5Dread(dataset_id, type_id, H5S_ALL, dataspace_id, H5P_DEFAULT, data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  return 0;
}

static int transpose(double *data){
  // u2 is aligned in 1 direction, which is not suitable to compute the longitudinal autocorrelation.
  // This subroutine aims to convert the row-major order to column-major order to improve the efficiency.
  double *buf = NULL;
  int i, j;
  ALLOCATE(buf, ITOT*JTOT*sizeof(double));
  for(j=0; j<JTOT; j++){
    for(i=0; i<ITOT; i++){
      buf[i*JTOT+j] = data[j*ITOT+i];
    }
  }
  memcpy(data, buf, ITOT*JTOT*sizeof(double));
  DEALLOCATE(buf);
  return 0;
}

static int subtract_average_1d(double *u1, double *u2){
  // Make the average of each row (u1) / column (u2) 0, so that the energy of zero frequency is forced to be 0
  double u1ave_1d, u2ave_1d;
  int i, j;
  for(j=0; j<JTOT; j++){
    u1ave_1d = 0.;
    for(i=0; i<ITOT; i++){
      u1ave_1d += u1[j*ITOT+i];
    }
    u1ave_1d /= ITOT;
    for(i=0; i<ITOT; i++){
      u1[j*ITOT+i] -= u1ave_1d;
    }
  }
  for(i=0; i<ITOT; i++){
    u2ave_1d = 0.;
    for(j=0; j<JTOT; j++){
      u2ave_1d += u2[i*JTOT+j];
    }
    u2ave_1d /= JTOT;
    for(j=0; j<JTOT; j++){
      u2[i*JTOT+j] -= u2ave_1d;
    }
  }
  return 0;
}

static int compute_two_point_correlations(double *R11, double *R22, const double *u1, const double *u2){
  // Compute the longitudinal two-point correlations R11 and R22
  // based on the equation (6.41), "Turbulent Flows" (S. B. Pope, 2001)
  double *u1_1d = NULL;
  double *u2_1d = NULL;
  printf("proceeding %s...", __func__);
  ALLOCATE(u1_1d, ITOT*sizeof(double));
  ALLOCATE(u2_1d, JTOT*sizeof(double));
  int i, j, shift;
  // R11
  for(j=0; j<JTOT; j++){
    for(shift=0; shift<ITOT; shift++){
      memcpy(u1_1d,            u1+j*ITOT+shift, (ITOT-shift)*sizeof(double));
      memcpy(u1_1d+ITOT-shift, u1+j*ITOT+0,           shift *sizeof(double));
      // compute autocorrelation
      for(i=0; i<ITOT; i++){
        R11[i] += u1_1d[i]*u1_1d[0];
      }
    }
  }
  for(i=0; i<ITOT; i++){
    R11[i] /= JTOT*ITOT;
  }
  // R22
  for(i=0; i<ITOT; i++){
    for(shift=0; shift<JTOT; shift++){
      memcpy(u2_1d,            u2+i*JTOT+shift, (JTOT-shift)*sizeof(double));
      memcpy(u2_1d+JTOT-shift, u2+i*JTOT+0,           shift *sizeof(double));
      // compute autocorrelation
      for(j=0; j<JTOT; j++){
        R22[j] += u2_1d[j]*u2_1d[0];
      }
    }
  }
  for(j=0; j<JTOT; j++){
    R22[j] /= ITOT*JTOT;
  }
  DEALLOCATE(u1_1d);
  DEALLOCATE(u2_1d);
  printf("done\n");
  return 0;
}

static int compute_one_dimensional_spectra(double *E11, double *E22, double *R11, double *R22){
  // Compute the one-dimensional spectra E_{ii}, which is "defined to be twice the one-dimensional Fourier transform of R_{ii}"
  // based on the equation (6.206), "Turbulent Flows" (S. B. Pope, 2001)
  fftw_plan plan = NULL;
  fftw_complex *buf = NULL;
  int i, j;
  //
  ALLOCATE(buf, (ITOT/2+1)*sizeof(fftw_complex));
  plan = fftw_plan_dft_r2c_1d(ITOT, R11, buf, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  // multiply pre-factors
  // 2 comes from "twice" the one-dimensional Fourier transform of R_{ii}
  // ITOT appears since FFTW3 does not normalize the transform
  // (ref: http://www.fftw.org/fftw3_doc/What-FFTW-Really-Computes.html)
  for(i=0; i<ITOT/2+1; i++){
    E11[i] = 2.*cabs(buf[i])/ITOT;
  }
  DEALLOCATE(buf);
  //
  ALLOCATE(buf, (JTOT/2+1)*sizeof(fftw_complex));
  plan = fftw_plan_dft_r2c_1d(JTOT, R22, buf, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  for(j=0; j<JTOT/2+1; j++){
    E22[j] = 2.*cabs(buf[j])/JTOT;
  }
  DEALLOCATE(buf);
  //
  return 0;
}

static int check_energy_consistency(const double *E11, const double *E22, const double *u1, const double *u2){
  // check whether the spectrum is computed correctly through
  // equation (6.209), "Turbulent Flows" (S. B. Pope, 2001)
  // volume-averaged energy in real space: <u^2>
  // volume-averaged energy in k    space: \int_0^{\infty} E_{ii}(k) dk
  double energy_u1, energy_u2;
  double energy_E11, energy_E22;
  int i, j;
  // in real space
  energy_u1 = 0.;
  for(j=0; j<JTOT; j++){
    for(i=0; i<ITOT; i++){
      energy_u1 += u1[j*ITOT+i]*u1[j*ITOT+i];
    }
  }
  energy_u1 /= ITOT*JTOT;
  energy_u2 = 0.;
  for(i=0; i<ITOT; i++){
    for(j=0; j<JTOT; j++){
      energy_u2 += u2[i*JTOT+j]*u2[i*JTOT+j];
    }
  }
  energy_u2 /= ITOT*JTOT;
  // in k space
  energy_E11 = 0.;
  for(i=0; i<ITOT/2+1; i++){
    energy_E11 += E11[i];
  }
  energy_E22 = 0.;
  for(j=0; j<JTOT/2+1; j++){
    energy_E22 += E22[j];
  }
  printf("energy computed from one-dimensional spectra\n");
  printf("\t%.7e %.7e\n", energy_E11, energy_E22);
  printf("energy computed from velocity\n");
  printf("\t%.7e %.7e\n", energy_u1, energy_u2);
  printf("rate\n");
  printf("\t%.7e %.7e\n", energy_E11/energy_u1, energy_E22/energy_u2);
  return 0;
}

int main(void){
  FILE *fp = NULL;
  hid_t file_id;
  double *u1 = NULL;
  double *u2 = NULL;
  double *R11 = NULL; // correlation in 1 dir, 1 vel
  double *R22 = NULL; // correlation in 2 dir, 2 vel
  double *E11 = NULL;
  double *E22 = NULL;
  int i, j;
  ALLOCATE(u1,   ITOT*JTOT*sizeof(double));
  ALLOCATE(u2,   ITOT*JTOT*sizeof(double));
  ALLOCATE(R11,       ITOT*sizeof(double));
  ALLOCATE(R22,       JTOT*sizeof(double));
  ALLOCATE(E11, (ITOT/2+1)*sizeof(double));
  ALLOCATE(E22, (JTOT/2+1)*sizeof(double));
  //
  file_id = H5Fopen("sample.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  read(file_id, "u", H5T_NATIVE_DOUBLE, u1);
  read(file_id, "v", H5T_NATIVE_DOUBLE, u2);
  H5Fclose(file_id);
  transpose(u2);
  //
  subtract_average_1d(u1, u2);
  compute_two_point_correlations(R11, R22, u1, u2);
  compute_one_dimensional_spectra(E11, E22, R11, R22);
  check_energy_consistency(E11, E22, u1, u2);
  //
  FOPEN(fp, "output/R11.dat", "w");
  for(i=0; i<ITOT; i++){
    fprintf(fp, "%.7f %.7e\n", i*DELTAX, R11[i]);
  }
  FCLOSE(fp);
  FOPEN(fp, "output/R22.dat", "w");
  for(j=0; j<JTOT; j++){
    fprintf(fp, "%.7f %.7e\n", j*DELTAY, R22[j]);
  }
  FCLOSE(fp);
  FOPEN(fp, "output/E11.dat", "w");
  for(i=0; i<ITOT/2+1; i++){
    fprintf(fp, "%4d %.7e\n", i, E11[i]);
  }
  FCLOSE(fp);
  FOPEN(fp, "output/E22.dat", "w");
  for(j=0; j<JTOT/2+1; j++){
    fprintf(fp, "%4d %.7e\n", j, E22[j]);
  }
  FCLOSE(fp);
  //
  DEALLOCATE(u1);
  DEALLOCATE(u2);
  DEALLOCATE(R11);
  DEALLOCATE(R22);
  DEALLOCATE(E11);
  DEALLOCATE(E22);
  return 0;
}

