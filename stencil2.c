#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(const int nx, const int ny, float * restrict  image, float *  tmp_image);
void init_image(const int nx, const int ny, float * restrict image, float *  tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *  image);
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
// #pragma omp simd
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

void stencil(const int nx, const int ny, float * restrict image, float * tmp_image) {
   float x = 0.6f;
   float y = 0.1f;  
for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {
      tmp_image[j+i*ny] = image[j+i*ny] * x;
      tmp_image[j+i*ny] += image[j  +(i-1)*ny] * y;
      tmp_image[j+i*ny] += image[j  +(i+1)*ny] * y;
      tmp_image[j+i*ny] += image[j-1+i*ny] * y;
      tmp_image[j+i*ny] += image[j+1+i*ny] * y;
    }
  }
// when i = 0; 
 for (int j = 1; j < ny-1; ++j) {
    tmp_image[j] = image[j] * x;
    tmp_image[j] += image[j + nx] * y;
    tmp_image[j] += image[j-1] * y;
    tmp_image[j] += image[j+1] * y; 
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
void init_image(const int nx, const int ny, float * restrict image, float * tmp_image) {
  
// Zero everything
  float zero = 0.0f;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      image[j+i*ny] = zero;
      tmp_image[j+i*ny] = zero;
    }
  }
  float val = 100.0f;
  // Checkerboard
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < 8; ++i) {
      for (int jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (int ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
          if ((i+j)%2)
          image[jj+ii*ny] = val;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float * image) {

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
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }
 float val = 255.0f;
  // Output image, converting to numbers 0-255
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      fputc((char)(val*image[j+i*ny]/maximum), fp);
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
