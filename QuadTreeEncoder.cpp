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


#include "count_ops.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>
#include <cstring>
#include <cstdint>
#include <stdint.h>
#include <omp.h>
#include "counters.h"
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "Encoder.h"
#include "QuadTreeEncoder.h"

extern int verb;
extern bool useYCbCr;

#define BUFFER_SIZE		(16)

QuadTreeEncoder::QuadTreeEncoder(int threshold, bool symmetry)
{
  this->threshold = threshold;
  this->symmetry = symmetry;
}

QuadTreeEncoder::~QuadTreeEncoder()
{
}


Transforms* QuadTreeEncoder::Encode(Image* source)
{
   //Initialize a hardware counter
    hwCounter_t cl;
    cl.init = false;
    initTicks(cl);
    
  Transforms* transforms = new Transforms;

  img.width = source->GetWidth();
  img.height = source->GetHeight();
  img.channels = source->GetChannels();
  transforms->channels = img.channels;

  omp_set_num_threads(N_THREADS);
  buffers = new PixelValue*[N_THREADS];
  for (int i = 0; i < N_THREADS; i++){
    buffers[i] = new PixelValue[BUFFER_SIZE * BUFFER_SIZE];
  }

  for (int channel = 1; channel <= img.channels; channel++)
    {
      INC_OP(2);

      // Load image into a local copy
      INC_OP(4);
      img.imagedata = new PixelValue[img.width * img.height];
      source->GetChannelData(channel, img.imagedata, img.width * img.height);

      INC_OP(4);
      if (img.width % 32 != 0 || img.height %32 != 0)
        {
          printf("Error: Image must have dimensions that are multiples of 32.\n");
          exit(-1);
        }
        

      // Make second channel the downsampled version of the image.
      //Get time before
        INC_OP(1);
        img.imagedata2 = IFSTransform::DownSample(img.imagedata, img.width, 0, 0, img.width / 2);

      // When using YCbCr we can reduce the quality of colour, because the eye
      // is more sensitive to intensity which is channel 1.
        INC_OP(2);
       if (channel >= 2 && useYCbCr){
         INC_OP(1);
        threshold *= 2;
        }

      // Go through all the range blocks
#pragma omp parallel for schedule(dynamic)
      for (int y = 0; y < img.height; y += BUFFER_SIZE)
        {
          INC_OP(2);
          for (int x = 0; x < img.width; x += BUFFER_SIZE)
            {
              INC_OP(3);
              //printf("****Buffer Size: %d\n", BUFFER_SIZE);
              findMatchesFor(transforms->ch[channel-1], x, y, BUFFER_SIZE);
              printf(".");
            }
          printf("\n");
        }
        
       //elapsed = getTicks(cl) - current_time;
       //printf("Number of Cycles required to take findBestMatch: %lu\n", elapsed);
        
      // Bring the threshold back to original.
      INC_OP(2);
      if (channel >= 2 && useYCbCr){
        INC_OP(1);
        threshold /= 2;
      }

      delete []img.imagedata2;
      img.imagedata2 = NULL;
      delete []img.imagedata;
      img.imagedata = NULL;
      printf("\n");
    }

  return transforms;
}

void QuadTreeEncoder::findMatchesFor(Transform& transforms, int toX, int toY, int blockSize)
{
    //Initialize a hardware counter
    //hwCounter_t cl;
    //cl.init = false;
    //initTicks(cl);
    
  int bestX = 0;
  int bestY = 0;
  int bestOffset = 0;
  IFSTransform::SYM bestSymmetry = IFSTransform::SYM_NONE;
  double bestScale = 0;
  double bestError = 1e9;

  
  //PixelValue* buffer = new PixelValue[blockSize * blockSize];
  
  // Get average pixel for the range block
  
  int rangeAvg = GetAveragePixel(img.imagedata, img.width, toX, toY, blockSize);
    
  
  // Go through all the downsampled domain blocks
    for (int y = 0; y < img.height; y += blockSize * 2)
    {
      INC_OP(3);
      for (int x = 0; x < img.width; x += blockSize * 2)
        {
          INC_OP(3);
          PixelValue* buffer = buffers[omp_get_thread_num()];

          for (int symmetry = 0; symmetry < IFSTransform::SYM_MAX; symmetry++)
            {
              INC_OP(2);
              IFSTransform::SYM symmetryEnum = (IFSTransform::SYM)symmetry;
              IFSTransform* ifs = new IFSTransform(x, y, 0, 0, blockSize, symmetryEnum, 1.0, 0);

              INC_OP(1);
              ifs->Execute(img.imagedata2, img.width / 2, buffer, blockSize, true);
                
              // Get average pixel for the downsampled domain block
              int domainAvg = GetAveragePixel(buffer, blockSize, 0, 0, blockSize);
                
              // Get scale and offset
              double scale = GetScaleFactor(img.imagedata, img.width, toX, toY, domainAvg,
                                            buffer, blockSize, 0, 0, rangeAvg, blockSize);
              INC_OP(2);
              int offset = (int)(rangeAvg - scale * (double)domainAvg);

              // Get error and compare to best error so far
              double error = GetError(buffer, blockSize, 0, 0, domainAvg,
                                      img.imagedata, img.width, toX, toY, rangeAvg, blockSize, scale);
              #pragma omp critical
              {
                INC_OP(1);
              if (error < bestError)
                {
                  bestError = error;
                  bestX = x;
                  bestY = y;
                  bestSymmetry = symmetryEnum;
                  bestScale = scale;
                  bestOffset = offset;
                }
              }

              delete ifs;

              if (!symmetry)
                break;
            }
        }
    }

    INC_OP(3);
  if (blockSize > 2 && bestError >= threshold)
    {
      // Recurse into the four corners of the current block.
      INC_OP(5);
      blockSize /= 2;
      findMatchesFor(transforms, toX, toY, blockSize);
      findMatchesFor(transforms, toX + blockSize, toY, blockSize);
      findMatchesFor(transforms, toX, toY + blockSize, blockSize);
      findMatchesFor(transforms, toX + blockSize, toY + blockSize, blockSize);
    }
  else
    {
      // Use this transformation
      IFSTransform* new_transform = new IFSTransform(
                                                     bestX, bestY,
                                                     toX, toY,
                                                     blockSize,
                                                     bestSymmetry,
                                                     bestScale,
                                                     bestOffset
                                                     );
      #pragma omp critical
      {
        transforms.push_back(new_transform);
      }
      INC_OP(1);
      if (verb >= 1)
        {
          printf("to=(%d, %d)\n", toX, toY);
          printf("from=(%d, %d)\n", bestX, bestY);
          printf("best error=%lf\n", bestError);
          printf("best symmetry=%d\n", (int)bestSymmetry);
          printf("best offset=%d\n", bestOffset);
          printf("best scale=%lf\n", bestScale);
        }
    }
}
