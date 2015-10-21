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

#ifndef IMAGE_H
#define IMAGE_H

typedef unsigned char PixelValue;

class ImageData
{
 public:
  ImageData();

  ~ImageData();

 public:
  int width;
  int height;
  int channels;
  PixelValue* imagedata;
  PixelValue* imagedata2;
  PixelValue* imagedata3;
};


class Image
{
 public:

  Image(string fileName);

  ~Image();

  void Load();

  void Save();

  void GetChannelData(int channel, PixelValue *buffer, int size);

  void SetChannelData(int channel, PixelValue *buffer, int size);

  int GetWidth();

  int GetHeight();

  int GetChannels();

  int GetOriginalSize();

 private:

  void ConvertToYCbCr(PixelValue& Y, PixelValue& Cb, PixelValue& Cr,
                      PixelValue R, PixelValue G, PixelValue B);

  void ConvertFromYCbCr(PixelValue& R, PixelValue& G, PixelValue& B,
                        PixelValue Y, PixelValue Cb, PixelValue Cr);

 private:

  string fileName;
  ImageData img;
  int originalSize;
};

#endif // IMAGE_H
