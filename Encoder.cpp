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

#define use_fast_GetError 1
#define use_fast_GetAveragePixel 1

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include "emmintrin.h"
#include <omp.h>
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "Encoder.h"

Encoder::Encoder(){
  temp_ints = new unsigned int* [N_THREADS];
  for (int i = 0; i < N_THREADS; i++){
    temp_ints[i] = new unsigned int[4];
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

#if use_fast_GetError
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
    diff = (int)(scale * (float)(*domain_ptr - domainAvg))
      - (*range_ptr - rangeAvg);
    top = (diff * diff);
    diff = (int)(scale * (float)(*(domain_ptr + 1) - domainAvg))
      - (*range_ptr - rangeAvg);
    top += (diff * diff);          
    diff = (int)(scale * (float)(*domain_ptr - domainAvg))
      - (*(range_ptr + 1) - rangeAvg);
    top += (diff * diff);
    diff = (int)(scale * (float)(*(domain_ptr + 1) - domainAvg))
      - (*(range_ptr + 1) - rangeAvg);
    top += (diff * diff);          
    return top/bottom;

  } else {

    __m128i top = _mm_setzero_si128();
    __m128i domainAvgVec = _mm_set1_epi32(domainAvg);
    __m128i rangeAvgVec = _mm_set1_epi32(rangeAvg);
    __m128 scaleVec = _mm_set1_ps(scale);

    __m128i domain;
    __m128i range;
    __m128 domain_ps;
    __m128i scaled;
    __m128i diff;
#define body(x) do{                                                     \
        domain = _mm_sub_epi32(_mm_load_si128((__m128i const *)(domain_ptr+x)), domainAvgVec); \
        range = _mm_sub_epi32(_mm_load_si128((__m128i const *)(range_ptr+x)), rangeAvgVec); \
        scaled = _mm_cvtps_epi32(_mm_mul_ps(scaleVec,                   \
                                            _mm_cvtepi32_ps(domain)));  \
        diff = _mm_sub_epi32(scaled, range);                            \
        top = _mm_add_epi32(top, _mm_mullo_epi32(diff, diff));          \
      } while (0);

    switch (size){

    case 16:
      for (int y = 0; y < size; y++){
        domain_ptr = domainData + (domainY + y) * domainWidth + domainX;
        range_ptr = rangeData + (rangeY + y) * rangeWidth + rangeX;            
        body(0);
        body(4);
        body(8);
        body(12);
      }
      break;

    case 8:
      for (int y = 0; y < size; y++){
        domain_ptr = domainData + (domainY + y) * domainWidth + domainX;
        range_ptr = rangeData + (rangeY + y) * rangeWidth + rangeX;
        body(0);
        body(4);
      }
      break;
    case 4:
      for (int y = 0; y < size; y++){
        domain_ptr = domainData + (domainY + y) * domainWidth + domainX;
        range_ptr = rangeData + (rangeY + y) * rangeWidth + rangeX;            
        body(0);
      }
      break;

    default:
      printf("ERROR: (default case) size=%d\n", size);
      exit(1);
    }
    unsigned int* temp_i = temp_ints[omp_get_thread_num()];
    _mm_store_si128((__m128i*)temp_i, top);
    return ((temp_i[0] + temp_i[1] + temp_i[2] + temp_i[3]) / bottom);
  }
}

#else //use_fast_GetError

double Encoder::GetError(
                         PixelValue* domainData, int domainWidth, int domainX, int domainY, int domainAvg,
                         PixelValue* rangeData, int rangeWidth, int rangeX, int rangeY, int rangeAvg,
                         int size, double scale)
{
  double top = 0;
  double bottom = (double)(size * size);

  for (int y = 0; y < size; y++)
    {
      for (int x = 0; x < size; x++)
        {
          int domain = (domainData[(domainY + y) * domainWidth + (domainX + x)] - domainAvg);
          int range = (rangeData[(rangeY + y) * rangeWidth + (rangeX + x)] - rangeAvg);
          int diff = (int)(scale * (double)domain) - range;

          // According to the formula we want (DIFF*DIFF)/(SIZE*SIZE)
          top += (diff * diff);

          if (top < 0)
            {
              printf("Error: Overflow occured during error %lf\n", top);
              exit(-1);
            }
        }
    }

  return (top / bottom);
}
#endif // use_fast_GetError

#if use_fast_GetAveragePixel

int Encoder::GetAveragePixel(PixelValue* domainData, int domainWidth,
                             int domainX, int domainY, int size)
{
    int top = 0;
    int bottom = (size * size);
    if(size == 2){
        top += domainData[domainY * domainWidth + domainX];
        top += domainData[domainY * domainWidth + domainX+1];
        top += domainData[(domainY+1) * domainWidth + domainX];
        top += domainData[(domainY+1) * domainWidth + domainX+1];
        
    }else{
        __m128i v1 = _mm_setzero_si128();
        __m128i v2 = _mm_setzero_si128();
        __m128i v3 = _mm_setzero_si128();
        __m128i v4 = _mm_setzero_si128();
        __m128i vsum = _mm_setzero_si128();
        switch(size){
            case 4:
                v1 = _mm_load_si128 ((__m128i const *)&domainData[domainY * domainWidth + domainX]);
                v2 = _mm_load_si128 ((__m128i const *)&domainData[(domainY+1) * domainWidth + domainX]);
                v3 = _mm_load_si128 ((__m128i const *)&domainData[(domainY+2) * domainWidth + domainX]);
                v4 = _mm_load_si128 ((__m128i const *)&domainData[(domainY+3) * domainWidth + domainX]);
                vsum = _mm_add_epi32(vsum, v1);
                vsum = _mm_add_epi32(vsum, v2);
                vsum = _mm_add_epi32(vsum, v3);
                vsum = _mm_add_epi32(vsum, v4);
                break;
            case 8:
                for(int y = domainY; y<domainY+size; y++){
                    v1 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX]);
                    v2 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX + 4]);
                    vsum = _mm_add_epi32(vsum, v1);
                    vsum = _mm_add_epi32(vsum, v2);
                }
                break;
                
            case 16:
                for(int y = domainY; y<domainY+size; y++){
                    v1 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX]);
                    v2 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX + 4]);
                    v3 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX + 8]);
                    v4 = _mm_load_si128 ((__m128i const *)&domainData[y * domainWidth + domainX + 12]);
                    vsum = _mm_add_epi32(vsum, v1);
                    vsum = _mm_add_epi32(vsum, v2);
                    vsum = _mm_add_epi32(vsum, v3);
                    vsum = _mm_add_epi32(vsum, v4);
                }
                break;
        }
        vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
        vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
        top = _mm_cvtsi128_si32(vsum);
    }
    return (top / bottom);
}

#else // use_fast_GetAveragePixel

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

#endif // use_fast_GetAveragePixel

