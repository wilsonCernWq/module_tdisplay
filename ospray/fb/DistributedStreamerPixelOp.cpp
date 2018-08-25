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

#include "DistributedStreamerPixelOp.h"
#include "ospray/fb/FrameBuffer.h"
#include "ospray/common/Data.h"
#include "DistributedStreamerPixelOp_ispc.h"

using namespace ospcommon;

using DInfo = ospray::DistributedStreamerPixelOp::TiledDisplayInfo;

namespace ospray {

//// DistributedStreamerPixelOp ///////////////////////////////////////////////

  DistributedStreamerPixelOp::DistributedStreamerPixelOp()
    : buf_count(0), buf_list(0, nullptr)
  {
  }

  void
  DistributedStreamerPixelOp::commit()
  {
    // parent class
    PixelOp::commit();

    // get MPI window
    win = static_cast<MPI_Win*>(getParam<void*>("windowObject", nullptr));

    // retrieve display information
    const Data *_d = getParamData("displayConfigurations", nullptr);
    dn = getParam<vec2i>("displayDimension", vec2i(0));
    ds = getParam<vec2i>("displayTileSize",  vec2i(0));
    di = reinterpret_cast<const DInfo*>(_d->data);
    df = static_cast<OSPFrameBufferFormat>
      (getParam1i("displayFormat", (int)OSP_FB_RGBA32F));

    // double check for users
    const auto nexp = dn.x * dn.y;
    const auto nget = static_cast<int32_t>(_d->numBytes / sizeof(DInfo));
    if (nexp != nget) {
      throw std::runtime_error("The 'displayConfigurations' array is shorter "
                               "than expected. " + std::to_string(nexp) + " "
                               "expected, found " + std::to_string(nget));
    }
  }

  std::string
  DistributedStreamerPixelOp::toString() const
  {
    return "ospray::DistributedStreamerPixelOp";
  }

  PixelOp::Instance*
  DistributedStreamerPixelOp::createInstance(FrameBuffer* fb,
                                             PixelOp::Instance* prev)
  {
    buf_list = std::vector<uint8_t*>(static_cast<size_t>(fb->getTotalTiles()), nullptr);
    return new DistributedStreamerPixelOp::Instance(fb, this);
  }

////////////////////

//// DistributedStreamerPixelOp::Instance /////////////////////////////////////

  DistributedStreamerPixelOp::Instance::
  Instance(FrameBuffer* fb, const DistributedStreamerPixelOp* pixelop)
    : op(pixelop)
  {
    this->fb = fb;
  }

  /*! gets called once at the beginning of the frame */
  void
  DistributedStreamerPixelOp::Instance::
  beginFrame()
  {
    MPI_CALL(Win_fence(MPI_MODE_NOPRECEDE, *(op->win)));
  }

  /*! gets called once at the end of the frame */
  void
  DistributedStreamerPixelOp::Instance::
  endFrame()
  {
    MPI_CALL(Win_fence(MPI_MODE_NOSTORE, *(op->win)));
    for (auto& b : op->buf_list) delete [] b;
    op->buf_list.clear();
    op->buf_count = 0;
  }

  void
  DistributedStreamerPixelOp::Instance::
  postAccum(Tile& tile)
  {
    // starting & ending index
    const vec2i si = tile.region.lower / op->ds;
    const vec2i ei = min((tile.region.lower + TILE_SIZE) / op->ds, op->dn-1);

    // convert tile into correct data format
    uint8_t* color;
    int32_t  unit;
    switch (op->df) {
    case OSP_FB_RGBA8:
      unit = 4;
      color = new uint8_t[TILE_SIZE * TILE_SIZE * unit];
      ispc::DSPOP_writeTile_RGBA8((ispc::VaryingTile*)&tile, reinterpret_cast<float*>(color));
      break;
    case OSP_FB_SRGBA:
      unit = 4;
      color = new uint8_t[TILE_SIZE * TILE_SIZE * unit];
      ispc::DSPOP_writeTile_SRGBA((ispc::VaryingTile*)&tile, reinterpret_cast<float*>(color));
      break;
    case OSP_FB_RGBA32F:
      unit = sizeof(float) * 4;
      color = new uint8_t[TILE_SIZE * TILE_SIZE * unit];
      ispc::DSPOP_writeTile_RGBA32F((ispc::VaryingTile*)&tile, reinterpret_cast<float*>(color));
      break;
    default:
      throw std::runtime_error("unknown framebuffer color format " + std::to_string(op->df) + ".");
    }
    op->buf_list[op->buf_count.fetch_add(1)] = color;

    // send tile to target
    for (auto dx = si.x; dx <= ei.x; ++dx) {
      for (auto dy = si.y; dy <= ei.y; ++dy) {
        const DInfo& disp = op->di[dx + dy * op->dn.x];
        const vec2i  ll = max( disp.idx      * op->ds, tile.region.lower); //< framebuffer space coord
        const vec2i  ur = min((disp.idx + 1) * op->ds, tile.region.upper); //  ^ same
        const vec2i  dim = ur - ll;                                        //< buffer dimension
        const vec2i  idx_disp = ll % op->ds;            //< index of the LL corner in the display space
        const vec2i  idx_tile = ll - tile.region.lower; //< index of the LL corner in the tile space
        for (auto py = 0; py < dim.y; ++py) {
          const auto off_disp = idx_disp.x + (idx_disp.y + py) * op->ds.x;
          const auto off_tile = idx_tile.x + (py + idx_tile.y) * TILE_SIZE;
          MPI_CALL(Put(&color[unit * off_tile],
                       unit * dim.x, MPI_CHAR, disp.rnk,
                       unit * off_disp + disp.off,
                       unit * dim.x, MPI_CHAR, *(op->win)));
        }
      }
    }
  }

  std::string
  DistributedStreamerPixelOp::Instance::
  toString() const
  {
    return "ospray::DistributedStreamerPixelOp::Instance";
  }

////////////////////

  OSP_REGISTER_PIXEL_OP(DistributedStreamerPixelOp, distributed_streamer);

} // ::ospray
