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

#ifndef QTE_H
#define QTE_H

class QuadTreeEncoder : public Encoder
{
 public:

  QuadTreeEncoder(int threshold = 100, bool symmetry = true);

  virtual ~QuadTreeEncoder();

  virtual Transforms* Encode(Image* source);

  PixelValue** buffers;

 protected:
  void findMatchesFor(Transform& transforms, int toX, int toY, int blockSize);
  void executeIFS(int blockSize);
  void calculateSummedAreaTable();
  int domainGetAveragePixel(int x, int y, int blockSize);

 protected:
  int threshold;
  bool symmetry;

  PixelValue *executePixels[4];
  PixelValue *averagePixels[4];
  PixelValue *table;
};

#endif // QTE_H
