#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <inttypes.h>
#include <arm_neon.h>

#define NUM_OF_CHANNELS 3

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

typedef struct
{
  unsigned int maxcolor, width, height;
  unsigned char *image_data;
} image;

void avg_to_gray(image *source, image *grayscale)
{

  if (grayscale->image_data == NULL)
  {
    perror("Memory allocation failed for grayscale image data");
    exit(EXIT_FAILURE);
  }
  unsigned int sum;
  for (int i = 0; i < grayscale->height * grayscale->width; i++)
  {
    sum = 0;
    for (int j = 0; j < NUM_OF_CHANNELS; j++)
    {
      sum += source->image_data[i * NUM_OF_CHANNELS + j];
    }
    grayscale->image_data[i] = sum / NUM_OF_CHANNELS;
  }
}
// Grayscale = 0.299R + 0.587G + 0.114B
void weighted_to_gray(image *source)
{
  image grayscale;
  grayscale.height = source->height;
  grayscale.width = source->width;
  grayscale.maxcolor = source->maxcolor;
  grayscale.image_data =
      malloc(grayscale.height * grayscale.width * sizeof(unsigned char));
  if (grayscale.image_data == NULL)
  {
    perror("Memory allocation failed for grayscale image data");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < grayscale.height * grayscale.width; i++)
  {
    grayscale.image_data[i] = 0.299 * source->image_data[i * 3] +
                              0.587 * source->image_data[i * 3 + 1] +
                              0.114 * source->image_data[i * 3 + 2];
  }

  memcpy(source->image_data, grayscale.image_data, grayscale.height * grayscale.width * sizeof(unsigned char));
  free(grayscale.image_data);
}

image read_image_from_file(FILE *fd)
{
  char buffor[256];

  image input_image;
  int i = 0;
  while (i < 3)
  {
    do
    {
      memset(buffor, 0, sizeof(buffor));
      fgets(buffor, 10, fd);
    } while (buffor[0] == '#');
    switch (i)
    {
    case 0:
      if (buffor[1] != '6' || buffor[0] != 'P')
      {
        perror("Supports only P6 format");
        exit(EXIT_FAILURE);
      }
      i++;
      break;
    case 1:
      if (2 !=
          sscanf(buffor, "%d %d", &input_image.width, &input_image.height))
      {
        perror("Width and height format error");
        exit(EXIT_FAILURE);
      }
      i++;
      break;
    case 2:
      if (1 != sscanf(buffor, "%d", &input_image.maxcolor))
      {
        perror("Maxval format error");
        exit(EXIT_FAILURE);
      }
      i++;
      break;
    }
  }
  input_image.image_data = malloc(input_image.height * input_image.width * 3 *
                                  sizeof(unsigned char));
  if (input_image.image_data == NULL)
  {
    perror("Could not allocate memory for image");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < input_image.height; i++)
  {
    size_t bytes_read =
        fread(&input_image.image_data[i * input_image.width * 3], 1,
              input_image.width * 3, fd);
    if (bytes_read != input_image.width * 3)
    {
      perror("Error reading from input file");
      exit(EXIT_FAILURE);
    }
  }
  return input_image;
}
int write_grayscale_image_to_file(FILE *fd, image *grayscale)
{

  fprintf(fd, "P5\n%d %d\n%d\n", grayscale->width, grayscale->height,
          grayscale->maxcolor);
  for (int i = 0; i < grayscale->height; i++)
  {
    size_t bytes_written = fwrite(&grayscale->image_data[i * grayscale->width],
                                  1, grayscale->width, fd);
    if (bytes_written != grayscale->width)
    {
      perror("Error writing to output file");
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}
void histogram_equalization(image *source)
{
  unsigned int histogram[source->maxcolor + 1],
      histogram_c[source->maxcolor + 1], g_min = source->maxcolor;
  memset(histogram, 0, sizeof(unsigned int) * (source->maxcolor + 1));
  memset(histogram_c, 0, sizeof(unsigned int) * (source->maxcolor + 1));

  for (int i = 0; i < source->height * source->width; i++)
  {
    histogram[source->image_data[i]]++;
  }

  for (int i = 0; i <= source->maxcolor; i++)
  {
    if (histogram[i] > 0)
    {
      g_min = i;
      break;
    }
  }

  histogram_c[0] = histogram[0];
  for (int i = 1; i <= source->maxcolor; i++)
  {
    histogram_c[i] = histogram_c[i - 1] + histogram[i];
  }

  unsigned int h_min = histogram_c[g_min];

  unsigned char t[source->maxcolor + 1];
  memset(t, 0, sizeof(unsigned char) * (source->maxcolor + 1));

  for (int i = 0; i <= source->maxcolor; i++)
  {
    double licz = histogram_c[i] - h_min;
    double mian = source->height * source->width - h_min;
    if (mian != 0)
    {
      t[i] = (licz / mian) * source->maxcolor;
    }
    else
    {
      t[i] = source->maxcolor;
    }
  }

  for (int i = 0; i < source->height * source->width; i++)
  {
    source->image_data[i] = t[source->image_data[i]];
  }
}

void otsu(image *source)
{
  unsigned int histogram[source->maxcolor + 1];
  float histogram_norm[source->maxcolor + 1];
  int mint, minsigma = INT_MAX;
  int sigma;
  memset(histogram, 0, sizeof(unsigned int) * (source->maxcolor + 1));

  for (int i = 0; i < source->height * source->width; i++)
  {
    histogram[source->image_data[i]]++;
  }

  for (int i = 0; i <= source->maxcolor; i++)
  {
    histogram_norm[i] = (float)histogram[i] / (source->height * source->width);
  }
  for (int t = 0; t <= source->maxcolor - 1; t++)
  {
    float omegaB = 0.0f, omegaF = 0.0f, sigmaF = 0.0f, sigmaB = 0.0f, expectedB = 0.0f, expectedF = 0.0f;
    for (int i = 0; i <= t; i++)
    {
      omegaB += histogram_norm[i];
      expectedB += histogram_norm[i] * i;
    }
    for (int i = 0; i <= t; i++)
    {
      sigmaB += histogram_norm[i] * pow(histogram[i] - expectedB, 2);
    }

    for (int i = t + 1; i <= source->maxcolor; i++)
    {
      omegaF += histogram_norm[i];
      expectedF += histogram_norm[i] * i;
    }
    for (int i = t + 1; i <= source->maxcolor; i++)
    {
      sigmaF += histogram_norm[i] * pow(histogram[i] - expectedF, 2);
    }
    sigma = sigmaF * omegaF + sigmaB * omegaB;
    if (sigma < minsigma)
    {
      mint = t;
      minsigma = sigma;
    }
  }
  for (int i = 0; i < source->height * source->width; i++)
  {
    source->image_data[i] = source->image_data[i] < mint ? 0 : source->maxcolor;
  }
}

void gamma_correction(image *source, float gamma)
{
  for (int i = 0; i < source->width * source->height; i++)
  {
    source->image_data[i] = (unsigned char)(pow((float)source->image_data[i] / (float)source->maxcolor, 1.0f / gamma) * (float)source->maxcolor);
  }
}

void mexican_hat_filter(image *source)
{
  int h[5][5] = {{0, 0, -1, 0, 0},
                 {0, -1, -2, -1, 0},
                 {-1, -2, 16, -2, -1},
                 {0, -1, -2, -1, 0},
                 {0, 0, -1, 0, 0}};
  image mexico;

  mexico.height = source->height;
  mexico.width = source->width;
  mexico.maxcolor = source->maxcolor;
  mexico.image_data =
      malloc(mexico.height * mexico.width * sizeof(unsigned char));

  memset(mexico.image_data, 0, mexico.height * mexico.width * sizeof(unsigned char));
  for (int i = 0; i < source->height - 5; i++)
  {
    for (int j = 0; j < source->width - 5; j++)
    {
      int sum = 0;
      for (int k = 0; k < 5; k++)
      {
        for (int l = 0; l < 5; l++)
        {
          sum += source->image_data[(i + k) * source->width + (j + l)] * h[k][l];
        }
      }
      mexico.image_data[(i + 2) * source->width + (j + 2)] = sum / 16;
    }
  }
  memcpy(source->image_data, mexico.image_data, mexico.height * mexico.width * sizeof(unsigned char));
  free(mexico.image_data);
}
void normal_sum(uint8_t src1[8], uint8_t src2[8], uint16_t *dest)
{
  for (int i = 0; i < 8; i++)
    dest[i] = src1[i] + src2[i];
}
void normal_mul(uint8_t src1[8], uint8_t src2[8], uint16_t *dest)
{
  for (int i = 0; i < 8; i++)
    dest[i] = src1[i] * src2[i];
}
void neon_sum(uint8x8_t src1, uint8x8_t src2, uint16_t *dest)
{
  uint16x8_t result = vmovl_u8(vadd_u8(src1, src2)); // Convert the result to uint16x8_t
  vst1q_u16(dest, result);                           // Use vst1q_u16 to store the result in dest
}

void neon_mul(uint8x8_t src1, uint8x8_t src2, uint16_t *dest)
{
  uint16x8_t result = vmull_u8(src1, src2); // Convert the result to uint16x8_t
  vst1q_u16(dest, result);                  // Use vst1q_u16 to store the result in dest
}
void neon_grayscale(image *source, image *grayscale)
{
  uint8x8_t r_ratio = vdup_n_u8(77);
  uint8x8_t g_ratio = vdup_n_u8(151);
  uint8x8_t b_ratio = vdup_n_u8(28);
  for (int i = 0; i < source->height * source->width; i += 8)
  {
    uint16x8_t temp;
    uint8x8x3_t rgb = vld3_u8(&source->image_data[i * 3]);
    uint8x8_t result;
    temp = vmull_u8(rgb.val[0], r_ratio);
    temp = vmlal_u8(temp, rgb.val[1], g_ratio);
    temp = vmlal_u8(temp, rgb.val[2], b_ratio);
    result = vshrn_n_u16(temp, 8);
    vst1_u8(&grayscale->image_data[i], result);
  }
}

void neon_mexican_hat_filter(image *source)
{
  unsigned char h[8][8] = {
      {0, 0, 1, 1, 1, 0, 0, 0},  // 3
      {0, 1, 3, 3, 3, 1, 0, 0},  // 11
      {1, 3, 0, 7, 0, 3, 1, 0},  // 15
      {1, 3, 7, 24, 7, 3, 1, 0}, // 46
      {1, 3, 0, 7, 0, 3, 1, 0},  // 15
      {0, 1, 3, 3, 3, 1, 0, 0},  // 11
      {0, 0, 1, 1, 1, 0, 0, 0}}; // 3
  // sum=29+46+29=104
  // Create an output image structure

  // unsigned char h[8][8] = {
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1},
  //     {1, 1, 1, 1, 1, 1, 1, 1}};

  // unsigned char h[8][8] = {
  //   {0,0,0,0,0,0,0,0},
  //   {0,0,0,0,0,0,0,0},
  //   {0,0,0,0,0,0,0,0},
  //   {0,0,0,0,1,0,0,0},
  //   {0,0,0,0,0,0,0,0},
  //   {0,0,0,0,0,0,0,0},
  //   {0,0,0,0,0,0,0,0}};
  image mexico;
  mexico.height = source->height;
  mexico.width = source->width;
  mexico.maxcolor = source->maxcolor;
  mexico.image_data = malloc(mexico.height * mexico.width * sizeof(unsigned char));

  uint8x8_t v_kernel[8];
  for (int i = 0; i < 8; i++)
  {
    int8x8_t temp = vld1_s8(h[i]);
    v_kernel[i] = vreinterpret_u8_s8(temp);
  }

  // Initialize the output image data to zero
  memset(mexico.image_data, 0, mexico.height * mexico.width * sizeof(unsigned char));

  for (int i = 0; i < source->height - 8; i++)
  {
    for (int j = 0; j < source->width - 8; j++)
    {
      int32_t suma = 0;
      uint16x8_t dest[8];
      for (int k = 0; k < 8; k++)
      {

        uint8x8_t v_src = vld1_u8(&source->image_data[(i + k) * source->width + j]);
        dest[k] = vmull_u8(v_kernel[k], v_src);
      }
      for (int s = 0; s < 8; s++)
      {
        suma += vgetq_lane_u16(dest[s], 0) + vgetq_lane_u16(dest[s], 1) + vgetq_lane_u16(dest[s], 2) + vgetq_lane_u16(dest[s], 3) +
                vgetq_lane_u16(dest[s], 4) + vgetq_lane_u16(dest[s], 5) + vgetq_lane_u16(dest[s], 6) + vgetq_lane_u16(dest[s], 7);
      }
      mexico.image_data[i * source->width + j] = (unsigned char)(suma / 104);
    }
  }

  memcpy(source->image_data, mexico.image_data, mexico.height * mexico.width * sizeof(unsigned char));

  free(mexico.image_data);
}

void erosion(image *source, int strength)
{
  image result;

  result.height = source->height;
  result.width = source->width;
  result.maxcolor = source->maxcolor;
  result.image_data = malloc(result.height * result.width * sizeof(unsigned char));

  memset(result.image_data, 0, result.height * result.width * sizeof(unsigned char));
  for (int i = 0; i < source->height - strength; i++)
  {
    for (int j = 0; j < source->width - strength; j++)
    {
      int f = INT_MAX;
      for (int k = 0; k < strength; k++)
      {
        for (int l = 0; l < strength; l++)
        {
          f = MIN(source->image_data[(i + k) * source->width + (j + l)], f);
        }
      }
      result.image_data[(i + 2) * source->width + (j + 2)] = f;
    }
  }
  memcpy(source->image_data, result.image_data, result.height * result.width * sizeof(unsigned char));
  free(result.image_data);
}
void dilation(image *source, int strength)
{
  image result;

  result.height = source->height;
  result.width = source->width;
  result.maxcolor = source->maxcolor;
  result.image_data = malloc(result.height * result.width * sizeof(unsigned char));

  memset(result.image_data, 0, result.height * result.width * sizeof(unsigned char));
  for (int i = 0; i < source->height - strength; i++)
  {
    for (int j = 0; j < source->width - strength; j++)
    {
      int f = INT_MIN;
      for (int k = 0; k < strength; k++)
      {
        for (int l = 0; l < strength; l++)
        {
          f = MAX(source->image_data[(i + k) * source->width + (j + l)], f);
        }
      }
      result.image_data[(i + 2) * source->width + (j + 2)] = f;
    }
  }
  memcpy(source->image_data, result.image_data, result.height * result.width * sizeof(unsigned char));
  free(result.image_data);
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    perror("Give input and output filenames as arguments");
    exit(EXIT_FAILURE);
  }
  FILE *fd = fopen(argv[argc - 2], "r");
  if (!fd)
  {
    perror("File could not be opened");
    exit(EXIT_FAILURE);
  }
  image input_image = read_image_from_file(fd);
  fclose(fd);

  printf("%dx%d %d\n", input_image.width, input_image.height,
         input_image.maxcolor);

  // Convert to grayscale
  image grayscale;
  grayscale.height = input_image.height;
  grayscale.width = input_image.width;
  grayscale.maxcolor = input_image.maxcolor;
  grayscale.image_data = malloc(grayscale.height * grayscale.width * sizeof(unsigned char));
  // avg_to_gray(&input_image, &grayscale);
  neon_grayscale(&input_image, &grayscale);
  // image grayscale = weighted_to_gray(&input_image);

  histogram_equalization(&grayscale);

  FILE *fd2 = fopen(argv[argc - 1], "wb");
  if (!fd2)
  {
    perror("Output file could not be opened");
    exit(EXIT_FAILURE);
  }
  free(input_image.image_data);
  // mexican_hat_filter(&grayscale);
  // neon_mexican_hat_filter(&grayscale);
  // otsu(&grayscale);

  printf("%d\n", argc);
  if (argc > 4)
  {
    float strength = atof(argv[2]);
    if (argv[1][0] == 'd')
      dilation(&grayscale, strength);
    else if (argv[1][0] == 'e')
      erosion(&grayscale, strength);
    else if (argv[1][0] == 'g')
      gamma_correction(&grayscale, strength);
  }
  else if (argc > 3)
  {
    if (argv[1][0] == 'o')
      otsu(&grayscale);
    else if (argv[1][0] == 'm')
      neon_mexican_hat_filter(&grayscale);
  }
  write_grayscale_image_to_file(fd2, &grayscale);
  // free the
  fclose(fd2);
  free(grayscale.image_data);

  return 0;
}