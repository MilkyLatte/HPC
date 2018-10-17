
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const int nx, const int ny, float * restrict image, float * restrict  tmp_image);
void init_image(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float * restrict image);
double wtime(void);

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // Allocate the image
  float *image = malloc(sizeof(float)*nx*ny);
  float *tmp_image = malloc(sizeof(float)*nx*ny);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  // Call the stencil kernel
  double tic = wtime();
  #pragma omp simd
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }
  double toc = wtime();


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
}

void stencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image) {
   float x = 0.6f;
   float y = 0.1f;  
for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {
      float variable = image[j+i*ny] * x;
      variable += image[j  +(i-1)*ny] * y;
      variable += image[j  +(i+1)*ny] * y;
      variable += image[j-1+i*ny] * y;
      variable += image[j+1+i*ny] * y;
      tmp_image[j+i*ny] = variable;
    }
  }
// when i = 0; 
 for (int j = 1; j < ny-1; ++j) {
    float variable = image[j] * x;
    variable += image[j + nx] * y;
    variable += image[j-1] * y;
    variable += image[j+1] * y;
    tmp_image[j] = variable; 
 }
// when i == nx - 1
 for (int j = 1; j < ny-1; ++j) {
    float variable = image[j+(nx-1)*ny] * x;
    variable += image[j +(nx-2)*ny] * y;
    variable += image[j-1+(nx-1)*ny] * y;
    variable += image[j+1+(nx-1)*ny] * y;
    tmp_image[j+(nx-1)*ny] = variable;
} 
// when j == 0
 for (int i = 1; i < nx-1; ++i){
   float variable = image[i*ny] * x;
   variable += image[(i-1)*ny] * y;
   variable += image[(i+1)*ny] * y;
   variable += image[1+i*ny] * y;
   tmp_image[i*ny] = variable;
}
// when j == nx -1
  for (int i = 1; i < nx-1; ++i) {
    float variable = image[(ny-1)+i*ny]*x;
    variable += image[(ny-1)+(i-1)*ny] * y;
    variable += image[(ny-1)+(i+1)*ny] * y;
    variable += image[(ny-1)-1+i*ny] * y;
    tmp_image[(ny-1)+i*ny] = variable;
}
    float leftBot = image[0] * x;
    leftBot += image[ny] * x;
    leftBot += image[1] * x;
    tmp_image[0] = leftBot;
    float rightBot = image[ny-1] * x;
    rightBot += image[(ny-1)+ny] * y;
    rightBot += image[(ny-1)-1] * y;
    tmp_image[ny-1] = rightBot; 
    float leftTop = image[(nx-1) * ny] * x;
    leftTop += image[((nx-1)-1) * ny] * y;
    leftTop += image[1+(nx-1)*ny] * y;
    tmp_image[(nx-1)*ny] = leftTop;
    float rightTop = image[(ny-1)+(nx-1) * ny] * x;
    rightTop += image[(ny-1) + (nx-1-1) * ny] * y;
    rightTop += image[((ny-1) - 1) + (nx-1) * ny] * y; 
    tmp_image[(ny-1)+(nx-1) * ny] = rightTop;

    }

// Create the input image
void init_image(const int nx, const int ny, float * restrict image, float * restrict  tmp_image) {
  // Zero everything
  #pragma omp simd
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int jj = j*ny*0.125f; jj < (j+1)*ny*0.125f; ++jj) {
        for (int ii = i*nx*0.125f; ii < (i+1)*nx*0.125f; ++ii) {
          if ((i+j)%2)
          
	image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float * restrict image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      fputc((char)(255*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
