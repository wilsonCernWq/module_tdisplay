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

#pragma once

// ospray
#include "ospray/fb/PixelOp.h"

// ospcommon
#include "ospcommon/vec.h"

// mpiCommon
#include "mpiCommon/MPICommon.h"

namespace ospray {

  /*! \brief Tool to gather composited images on customized ranks instead of 
   * master. 
   */
  struct OSPRAY_SDK_INTERFACE DistributedStreamerPixelOp : public PixelOp
  {

    struct OSPRAY_SDK_INTERFACE TiledDisplayInfo
    {
      int32_t          rnk; //< receving rank
      size_t           off; //< byte offset in the target data window
      ospcommon::vec2i idx; //< location of the tiled display
    };

    struct OSPRAY_SDK_INTERFACE Instance : public PixelOp::Instance
    {
      Instance(FrameBuffer*, const DistributedStreamerPixelOp*);

      /*! gets called every time the frame buffer got 'commit'ted */
      void  commitNotify() override {}

      /*! gets called once at the beginning of the frame */
      void beginFrame() override;

      /*! gets called once at the end of the frame */
      void endFrame() override;

      /*! called whenever a new tile comes in from a renderer, but
          _before_ the tile gets written/accumulated into the frame
          buffer. this way we can, for example, fill in missing
          samples; however, the tile will _not_ yet contain the
          previous frame's contributions from the accum buffer
          etcpp. In distriubuted mode, it is undefined if this op gets
          executed on the node that _produces_ the tile, or on the
          node that _owns_ the tile (and its accum buffer data)  */
      void preAccum(Tile&) override {}

      /*! called right after the tile got accumulated; i.e., the
          tile's RGBA values already contain the accu-buffer blended
          values (assuming an accubuffer exists), and this function
          defines how these pixels are being processed before written
          into the color buffer */
      void postAccum(Tile&) override;

      //! \brief common function to help printf-debugging
      /*! Every derived class should overrride this! */
      std::string toString() const override;

      // Data members //
      const DistributedStreamerPixelOp* op;
    };

    DistributedStreamerPixelOp();
    
    /*! ospray commit */
    void commit() 
      override;

    //! \brief common function to help printf-debugging
    /*! Every derived class should overrride this! */    
    std::string toString() const 
      override;

    //! \brief create an instance of this pixel op
    PixelOp::Instance* createInstance(FrameBuffer*, PixelOp::Instance*)
      override;

    // Data members //
    const TiledDisplayInfo* di;  // detailed info for each tiled display
    ospcommon::vec2i        dn;  // number (Nx x Ny) of tiled display
    ospcommon::vec2i        ds;  // width and height for each tiled display
    OSPFrameBufferFormat    df;  // color format
    MPI_Win*                win; // MPI window object

    // (Internal) Data members //
    mutable std::atomic<size_t>   buf_count;
    mutable std::vector<uint8_t*> buf_list;
    mutable std::mutex lock;
  };

} // ::ospray
