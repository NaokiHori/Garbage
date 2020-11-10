#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>


static int write_data(const int itot, const int jtot){
  const char fname[] = {"test.h5"};
  const int rank = 2;
  const hsize_t dims[2] = {jtot, itot};
  hid_t file_id, dataspace_id, dataset_id;
  herr_t status;
  double *data = NULL;
  int i, j;
  data = (double*)calloc(itot*jtot, sizeof(double));
  for(j=0; j<jtot; j++){
    for(i=0; i<itot; i++){
      data[j*itot+i] = (double)(10*j+i);
    }
  }
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(rank, dims, dims);
  dataset_id = H5Dcreate(file_id, "/dset", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5S_ALL, data);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
  free(data);
  return 0;
}

static int write_coordinate(const int itot, const int jtot){
  const char fname[] = {"coordinate.h5"};
  const int rank = 2;
  const hsize_t dims[2] = {jtot, itot};
  hid_t file_id, dataspace_id, dataset_id;
  herr_t status;
  double *xs= NULL;
  double *ys= NULL;
  int i, j;
  xs = (double*)calloc(itot*jtot, sizeof(double));
  ys = (double*)calloc(itot*jtot, sizeof(double));
  for(j=0; j<jtot; j++){
    for(i=0; i<itot; i++){
      xs[j*itot+i] = 1.*i;
    }
  }
  for(j=0; j<jtot; j++){
    for(i=0; i<itot; i++){
      ys[j*itot+i] = 1.*j;
    }
  }
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id = H5Screate_simple(rank, dims, dims);
  dataset_id = H5Dcreate(file_id, "/x", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5S_ALL, xs);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  dataspace_id = H5Screate_simple(rank, dims, dims);
  dataset_id = H5Dcreate(file_id, "/y", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5S_ALL, ys);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
  free(xs);
  free(ys);
  return 0;
}

static int make_xdmf(const int itot, const int jtot){
  FILE *fp = NULL;
  const char fname[] = {"test.xmf"};
  if((fp=fopen(fname, "w"))==NULL){
    printf("file open error\n");
    exit(1);
  }
  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"2.0\">\n");
  fprintf(fp, "<Domain>\n");
  fprintf(fp, "<Grid Name=\"GridName\" GridType=\"Uniform\">\n");
  fprintf(fp, "<Topology TopologyType=\"2DSMesh\" NumberOfElements= \"%d %d\"/>\n", jtot, itot);
  fprintf(fp, "<Geometry GeometryType=\"X_Y\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", jtot, itot);
  fprintf(fp, "coordinate.h5:/x\n");
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", jtot, itot);
  fprintf(fp, "coordinate.h5:/y\n");
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Geometry>\n");
  fprintf(fp, "<Attribute Name=\"data\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fp, "<DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", jtot, itot);
  fprintf(fp, "test.h5:/dset\n");
  fprintf(fp, "</DataItem>\n");
  fprintf(fp, "</Attribute>\n");
  fprintf(fp, "</Grid>\n");
  fprintf(fp, "</Domain>\n");
  fprintf(fp, "</Xdmf>\n");
  fclose(fp);
  return 0;
}

int main(void){
  const int itot = 6;
  const int jtot = 4;
  write_data(itot, jtot);
  write_coordinate(itot, jtot);
  make_xdmf(itot, jtot);
  return 0;
}

