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
#include <string>
#include <vector>
using namespace std;

#include "Image.h"
#include "IFSTransform.h"
#include "count_ops.h"

extern int verb;

/////////////////////////////////////////////////////////////////////
// class Transforms

Transforms::Transforms()
{
  channels = 0;
}

Transforms::~Transforms()
{
  for (int i = 0; i < channels; i++)
    {
      for (int j = 0; j < ch[i].size(); j++)
        delete (ch[i][j]);
      ch[i].clear();
    }
}


/////////////////////////////////////////////////////////////////////
// class IFSTransform

PixelValue* IFSTransform::DownSample(PixelValue* src, int srcWidth,
                                     int startX, int startY, int targetSize)
{
  INC_OP(1);
  PixelValue* dest = new PixelValue[targetSize * targetSize];
  int destX = 0;
  int destY = 0;

  for (int y = startY; y < startY + targetSize * 2; y += 2)
    {
     INC_OP(4);
      for (int x = startX; x < startX + targetSize * 2; x += 2)
        {
         INC_OP(4);
          // Perform simple 2x2 average
          int pixel = 0;
          INC_OP(17);
          pixel += src[y * srcWidth + x];
          pixel += src[y * srcWidth + (x + 1)];
          pixel += src[(y + 1) * srcWidth + x];
          pixel += src[(y + 1) * srcWidth + (x + 1)];
          pixel /= 4;

          INC_OP(3);
          dest[destY * targetSize + destX] = pixel;
          destX++;
        }
      INC_OP(1);
      destY++;
      destX = 0;
    }

  return dest;
}

IFSTransform::IFSTransform(int fromX, int fromY, int toX, int toY, int size,
                           IFSTransform::SYM symmetry, double scale, int offset)
{
  this->fromX = fromX;
  this->fromY = fromY;
  this->toX = toX;
  this->toY = toY;
  this->size = size;
  this->symmetry = symmetry;
  this->scale = scale;
  this->offset = offset;
}

IFSTransform::~IFSTransform()
{
}

void IFSTransform::Execute(PixelValue* src, int srcWidth,
                           PixelValue* dest, int destWidth, bool downsampled)
{
  INC_OP(2);
  int fromX = this->fromX / 2;
  int fromY = this->fromY / 2;
  int dX = 1;
  int dY = 1;
  bool inOrder = isScanlineOrder();

  INC_OP(1);
  if (!downsampled)
    {
      PixelValue* newSrc = DownSample(src, srcWidth, this->fromX, this->fromY, size);
      src = newSrc;
      srcWidth = size;
      fromX = fromY = 0;
    }

  INC_OP(1);
  if (!isPositiveX())
    {
      INC_OP(3); 
      fromX += size - 1;
      dX = -1;
    }

  INC_OP(1);
  if (!isPositiveY())
    {
      INC_OP(3);
      fromY += size - 1;
      dY = -1;
    }

  int startX = fromX;
  int startY = fromY;

  for (int toY = this->toY; toY < (this->toY + size); toY++)
    {
     INC_OP(3);
      for (int toX = this->toX; toX < (this->toX + size); toX++)
        {
          INC_OP(4);
          if (verb >= 4)
            {
              printf("toX=%d\n", toX);
              printf("toY=%d\n", toY);
              printf("fromX=%d\n", fromX);
              printf("fromY=%d\n", fromY);
            }
          INC_OP(4);  
          int pixel = src[fromY * srcWidth + fromX];
          pixel = (int)(scale * pixel) + offset;

          INC_OP(2);
          if (pixel < 0)
            pixel = 0;
          if (pixel > 255)
            pixel = 255;

          INC_OP(1);
          if (verb >= 4)
            printf("pixel=%d\n", pixel);

          INC_OP(2);
          dest[toY * destWidth + toX] = pixel;

         if (inOrder){
            INC_OP(1);
            fromX += dX;
         }else{
            INC_OP(1);
            fromY += dY;
         }
        }

      if (inOrder)
        {
          INC_OP(1);
          fromX = startX;
          fromY += dY;
        }
      else
        {
          INC_OP(1);
          fromY = startY;
          fromX += dX;
        }
    }

  INC_OP(1);
  if (!downsampled)
    {
      delete []src;
      src = NULL;
    }
}

bool IFSTransform::isScanlineOrder()
{
  INC_OP(4);
  return (
          symmetry == SYM_NONE ||
          symmetry == SYM_R180 ||
          symmetry == SYM_HFLIP ||
          symmetry == SYM_VFLIP
          );
}

bool IFSTransform::isPositiveX()
{
  INC_OP(4);
  return (
          symmetry == SYM_NONE ||
          symmetry == SYM_R90 ||
          symmetry == SYM_VFLIP ||
          symmetry == SYM_RDFLIP
          );
}

bool IFSTransform::isPositiveY()
{
  INC_OP(4);
  return (
          symmetry == SYM_NONE ||
          symmetry == SYM_R270 ||
          symmetry == SYM_HFLIP ||
          symmetry == SYM_RDFLIP
          );
}

