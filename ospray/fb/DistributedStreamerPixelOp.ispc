// ======================================================================== //
// Copyright 2009-2018 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "fb/FrameBuffer.ih"
#include "fb/Tile.ih"
#include "math/vec.ih"

inline vec4f soa_to_aos4f(const vec4f s)
{
  vec4f a;
  soa_to_aos4(s.x, s.y, s.z, s.w, (uniform float*uniform)&a);
  return a;
}

// writes the final normalized colors in a tile into a color buffer
#define template_writeTile(name, dst_fmt, cvt)                                  \
export void DSPOP_writeTile_##name(const VaryingTile *uniform tile,             \
                                   void *uniform _color)                        \
{                                                                               \
  const uniform vec2i size = box_size(tile->region);                            \
  const uniform int maxi = size.y*(TILE_SIZE/programCount);                     \
  dst_fmt *uniform color = (dst_fmt*uniform)_color;                             \
  for (uniform int i = 0; i < maxi; i++) {                                      \
    const vec4f c = make_vec4f(tile->r[i], tile->g[i], tile->b[i], tile->a[i]); \
    color[i] = cvt(c);                                                          \
  }                                                                             \
}
template_writeTile(RGBA8, varying uint32, cvt_uint32);
template_writeTile(SRGBA, varying uint32, linear_to_srgba8);
template_writeTile(RGBA32F, varying vec4f, soa_to_aos4f);
#undef template_writeTile
