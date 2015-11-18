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
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "Encoder.h"
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>


//function to print SIMD Vector
void print128_num(__m128i var)
{
    uint8_t *val = (uint8_t*) &var;
    printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
           val[0], val[1], val[2], val[3], val[4], val[5],
           val[6], val[7], val[8], val[9],val[10],val[11],val[12],val[13],val[14],val[15]);
}

//*****
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
/*****/
int Encoder::GetAveragePixel(PixelValue* domainData, int domainWidth,
                             int domainX, int domainY, int size)
{
    //printf("block size: %d\n", size);
    float top = 0;
    int bottom = (size * size);
    //printf("Buffer Size SIMD: %d\n", size);
    
    // Simple average of all pixels.
    __m128i v;
    __m128i tempSum = _mm_setzero_si128();
    
    
    if(size>=4){
        for( int y = domainY; y< domainY +size; y++){
            for (int x = domainX; x < domainX + size; x+=4){
                
                //Then first unpack those 8-bit values into 16-bit values in the lower 64 bits of the register, interleaving them with 0s:
                v = _mm_cvtsi32_si128(*(const int*)&domainData[y * domainWidth + x]);
                //print128_num(v);
                
                //And again unpack those 16-bit values into 32-bit values:
                v = _mm_unpacklo_epi8(v, _mm_setzero_si128());
                //print128_num(v);
                
                //now have each pixel as 32-bit integer in the respective 4 components of the SSE register.
                v = _mm_unpacklo_epi16(v, _mm_setzero_si128());
                //print128_num(v);
                
                tempSum = _mm_add_epi32(tempSum,v);
                //printf("-----Sum-----\n");
                //print128_num(tempSum);
                //printf("-------------\n");
            }
        }
        
        //Add 4 elements of SIMD vector into one single element.
        tempSum = _mm_hadd_epi32(tempSum,tempSum);
        tempSum = _mm_hadd_epi32(tempSum,tempSum);
        top = _mm_cvtsi128_si32(tempSum);
        
        if (top < 0)
        {
            printf("Error: Accumulator rolled over averaging pixels.\n");
            exit(-1);
        }
        
    }else{
        
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
    }
    return (top / bottom);
}
