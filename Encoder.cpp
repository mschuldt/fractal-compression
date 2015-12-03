/*
 * Fractal Image Compression. Copyright 2004 Alex Kennberg.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <omp.h>
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "Encoder.h"

Encoder::Encoder(){
  temp_ints = new unsigned int* [N_THREADS];
  temp_floats = new float* [N_THREADS];
  for (int i = 0; i < N_THREADS; i++){
    temp_ints[i] = new unsigned int[4];
    temp_floats[i] = new float[4];
  }
}

double Encoder::GetScaleFactor(
                               PixelValue* domainData, int domainWidth, int domainX, int domainY, int domainAvg,
                               PixelValue* rangeData, int rangeWidth, int rangeX, int rangeY, int rangeAvg,
                               int size)
{
  int top = 0;
  int bottom = 0;

  for (int y = 0; y < size; y++)
    {
      for (int x = 0; x < size; x++)
        {
          int domain = (domainData[(domainY + y) * domainWidth + (domainX + x)] - domainAvg);
          int range = (rangeData[(rangeY + y) * rangeWidth + (rangeX + x)] - rangeAvg);

          // According to the formula we want (R*D) / (D*D)
          top += range * domain;
          bottom += domain * domain;

          if (bottom < 0)
            {
              printf("Error: Overflow occured during scaling %d %d %d %d\n",
                     y, domainWidth, bottom, top);
              exit(-1);
            }
        }
    }

  if (bottom == 0)
    {
      top = 0;
      bottom = 1;
    }

  return ((double)top) / ((double)bottom);
}

void print128i(__m128i var)
{
    unsigned int *val = (unsigned int*) &var;
    printf("__m128i: %d %d %d %d\n", val[0], val[1], val[2], val[3]);
}

void print128(__m128 var)
{
    float *val = (float*) &var;
    printf("__m128: %f %f %f %f\n", val[0], val[1], val[2], val[3]);
}

double Encoder::GetError(
                         PixelValue* domainData, int domainWidth, int domainX, int domainY, int domainAvg,
                         PixelValue* rangeData, int rangeWidth, int rangeX, int rangeY, int rangeAvg,
                         int size, double scale)
  {
    // printf("domainAvg = %d\nrangeAvg = %d\nscale = %f\n", domainAvg, rangeAvg, scale);
    float bottom = (float)(size * size);
    PixelValue * domain_ptr;
    PixelValue * range_ptr;

    if (size == 2){
      int top;
      int domain, range, diff;
      domain_ptr = domainData + (domainY) * domainWidth + domainX;
      range_ptr = rangeData + (rangeY) * rangeWidth + rangeX;          
      diff = (int)(scale * (double)(*domain_ptr - domainAvg))
        - (*range_ptr - rangeAvg);
      top = (diff * diff);
      diff = (int)(scale * (double)(*(domain_ptr + 1) - domainAvg))
        - (*range_ptr - rangeAvg);
      top += (diff * diff);          
      domain_ptr = domainData + (domainY) * domainWidth + domainX;
      range_ptr = rangeData + (rangeY) * rangeWidth + rangeX;          
      diff = (int)(scale * (double)(*domain_ptr - domainAvg))
        - (*(range_ptr + 1) - rangeAvg);
      top += (diff * diff);
      diff = (int)(scale * (double)(*(domain_ptr + 1) - domainAvg))
        - (*(range_ptr + 1) - rangeAvg);
      top += (diff * diff);          
      return top/bottom;

    } else {

      __m128i top = _mm_setzero_si128();

      unsigned int* temp_i = temp_ints[omp_get_thread_num()];
      float* temp_f = temp_floats[omp_get_thread_num()];

      temp_i[0] = domainAvg;
      temp_i[1] = domainAvg;
      temp_i[2] = domainAvg;
      temp_i[3] = domainAvg;
      __m128i domainAvgVec = _mm_load_si128((__m128i const *)temp_i);

      temp_i[0] = rangeAvg;
      temp_i[1] = rangeAvg;
      temp_i[2] = rangeAvg;
      temp_i[3] = rangeAvg;
      __m128i rangeAvgVec = _mm_load_si128((__m128i const *)temp_i);

      temp_f[0] = scale;
      temp_f[1] = scale;
      temp_f[2] = scale;
      temp_f[3] = scale;
      __m128 scaleVec = _mm_load_ps(temp_f);

      switch (size){
      case 16:
      case 8:
      case 4:
        for (int y = 0; y < size; y++)
          {
            domain_ptr = domainData + (domainY + y) * domainWidth + domainX;
            range_ptr = rangeData + (rangeY + y) * rangeWidth + rangeX;
            for (int x = 0; x < size; x+=4)
              {
                __m128i domain = _mm_load_si128((__m128i const *)(domain_ptr + x));
                domain = _mm_sub_epi32(domain, domainAvgVec);
                __m128i range = _mm_load_si128((__m128i const *)(range_ptr + x));
                range = _mm_sub_epi32(range, rangeAvgVec);
                __m128 domain_ps = _mm_cvtepi32_ps(domain);
                __m128i scaled = _mm_cvtps_epi32(_mm_mul_ps(scaleVec, domain_ps));
                __m128i diff = _mm_sub_epi32(scaled, range);
                top = _mm_add_epi32(top, _mm_mullo_epi32(diff, diff));
              }
          }
        break;

      default:
        printf("default case, size=%d\n", size);
        exit(1);
      }

      _mm_store_si128((__m128i*)temp_i, top);
      return ((temp_i[0] + temp_i[1] + temp_i[2] + temp_i[3]) / bottom);
    }
  }


int Encoder::GetAveragePixel(PixelValue* domainData, int domainWidth,
                             int domainX, int domainY, int size)
{
  int top = 0;
  int bottom = (size * size);

  // Simple average of all pixels.
  for (int y = domainY; y < domainY + size; y++)
    {
      for (int x = domainX; x < domainX + size; x++)
        {
          top += domainData[y * domainWidth + x];

          if (top < 0)
            {
              printf("Error: Accumulator rolled over averaging pixels.\n");
              exit(-1);
            }
        }
    }

  return (top / bottom);
}
