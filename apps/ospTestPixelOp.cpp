//
// Created by qadwu on 8/6/18.
//

#include "ospcommon/vec.h"
#include "ospcommon/box.h"
#include "ospcommon/utility/SaveImage.h"
#include "ospray/ospray.h"
#include <vector>

#define USE_MPI 1
#define USE_PIXEL_OP (1 & USE_MPI)

#if USE_MPI
#include <mpi.h>
#endif

#define debug(x) if (world_rank == x) std::cout

using namespace ospcommon;

///// Class Definitions ///////////////////////////////////////////////////////

// struct DistributedRegion {
//   box3f bd;
//   int id = 0;
//   DistributedRegion(const box3f &b, const int &i) : bd(b), id(i) {}
// };

#if USE_PIXEL_OP
struct DisplayInfo {
  int32_t rank;
  size_t byte_offset;
  // ^ the byte offset in the receiving buffer, which can non-zero if 
  //   one wants to allow one tiled display receiving multiple images.
  vec2i tile_offset;
  // ^ the location (in terms of index) of this tiled display 
};
#endif

////////////////////

///// Static Variables ////////////////////////////////////////////////////////

/// @name Command Line Arguments
/// @{
static vec2i fb_size{300, 300};
static enum { DATA_RANK, DATA_SPHERE } data_mode;
static int32_t num_of_image_per_display = 3;
static int32_t num_of_displays = 2;
static int32_t num_of_frames = 10;
/// @}

/// @name MPI Information
/// @{
static int world_size = 1; //< Get the number of processes
static int world_rank = 0; //< Get the rank of the process
/// @}

////////////////////

///// Initialize OSPRay ///////////////////////////////////////////////////////

static void OSPContext_ErrorFunc(OSPError, const char *msg) {
  std::cerr << "#osp: (rank " << world_rank << ")" << msg;
}
static void OSPContext_StatusFunc(const char *msg) {
  std::cerr << "#osp: (rank " << world_rank << ")" << msg;
}

void InitOSPRay() {
  OSPError err;
  err = ospLoadModule("ispc");
  if (err != OSP_NO_ERROR) {
    std::cerr << "[Error] can't load ispc module" << std::endl;
  }
  err = ospLoadModule("mpi");
  if (err != OSP_NO_ERROR) {
    std::cerr << "[Error] can't load mpi module" << std::endl;
  }
  err = ospLoadModule("tdisplay");
  if (err != OSP_NO_ERROR) {
    std::cerr << "[Error] can't load tdisplay module" << std::endl;
  }
#if USE_MPI
  OSPDevice device = ospNewDevice("mpi_distributed");
  ospDeviceSet1i(device, "masterRank", 0);
#else
  OSPDevice device = ospNewDevice("default");
#endif
  ospDeviceSetErrorFunc(device, OSPContext_ErrorFunc);
  ospDeviceSetStatusFunc(device, OSPContext_StatusFunc);
  ospDeviceCommit(device);
  ospSetCurrentDevice(device);
}

////////////////////

///// Helper Functions ////////////////////////////////////////////////////////
namespace ospcommon {
namespace utility {
void writePPM(const std::string &fileName,
              const int sizeX, const int sizeY,
              const float *pixel) {
  FILE *file = fopen(fileName.c_str(), "wb");
  if (file == nullptr)
    throw std::runtime_error("Can't open file for writeP[FP]M!");
  fprintf(file, "P6\n%i %i\n255\n", sizeX, sizeY);
  auto out = STACK_BUFFER(uint8_t, 3 * sizeX);
  for (int y = 0; y < sizeY; y++) {
    const auto *in = &pixel[4 * (sizeY - 1 - y) * sizeX];
    for (int x = 0; x < sizeX; x++)
      for (int c = 0; c < 3; c++)
        out[3 * x + c] = (uint8_t) (in[4 * x + c] * 255.f);
    fwrite(out, 3 * sizeX, sizeof(uint8_t), file);
  }
  fprintf(file, "\n");
  fclose(file);
}
};
};
////////////////////

///// Main ////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  ///-- Initialize MPI ----------------------------------------------------////
#if USE_MPI
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  assert(provided == MPI_THREAD_MULTIPLE);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif


  ///-- Initialize OSPRay -------------------------------------------------////
  InitOSPRay();


  ///-- Debug -------------------------------------------------------------////
#if USE_MPI && 0 // Wait for debugger ...
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    std::cout << "please attach a debugger to a process, "
          << "and then type continue ..."
          << std::endl;
    std::string tmp; std::cin >> tmp;
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif


  ///-- Parse Arguments ---------------------------------------------------////
  size_t num_of_bricks_per_rank = 1;
  vec3f vp(-128, -128, 512), vu(0, 1, 0), vi(128, 128, 128);
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-num-bricks") {
      num_of_bricks_per_rank = size_t(std::stoi(argv[++i]));
    } else if (arg == "-vp") {
      vp.x = std::stof(argv[++i]);
      vp.y = std::stof(argv[++i]);
      vp.z = std::stof(argv[++i]);
    } else if (arg == "-vu") {
      vu.x = std::stof(argv[++i]);
      vu.y = std::stof(argv[++i]);
      vu.z = std::stof(argv[++i]);
    } else if (arg == "-vi") {
      vi.x = std::stof(argv[++i]);
      vi.y = std::stof(argv[++i]);
      vi.z = std::stof(argv[++i]);
    } else if (arg == "-fb") {
      fb_size.x = std::stoi(argv[++i]);
      fb_size.y = std::stoi(argv[++i]);
    } else if (arg == "-tile-config") {
      num_of_image_per_display = std::stoi(argv[++i]);
      num_of_displays = std::stoi(argv[++i]);
    } else if (arg == "-data") {
      const std::string str = argv[++i];
      try {
        const int idx = std::stoi(str);
        if (idx == 0) {
          data_mode = DATA_RANK;
        } else if (idx == 1) {
          data_mode = DATA_SPHERE;
        } else {
          data_mode = DATA_RANK;
        }
      }
      catch (const std::invalid_argument &) {
        if (str == "rank") {
          data_mode = DATA_RANK;
        } else if (str == "sphere") {
          data_mode = DATA_SPHERE;
        } else {
          data_mode = DATA_RANK;
        }
      }
    } else if (arg == "-frames") {
      num_of_frames = std::stoi(argv[++i]);
    } else {
      std::cerr << "unknown argument: " << arg << std::endl;
    }
  }
  if (num_of_displays > world_size) {
    throw std::runtime_error("number of tiled display = " +
        std::to_string(num_of_displays) +
        ", which cannot be larger than the number "
        "of MPI ranks. Please rerun the program with "
        "more MPI ranks.");
  }


  ///-- Generate Volume Data ----------------------------------------------////
  const vec_t<size_t, 3> g_dim(256, 256, 256);
  const vec_t<size_t, 3> l_dim(3 + g_dim.x / num_of_bricks_per_rank,
                               3 + g_dim.y / world_size, g_dim.z);
  const size_t l_num_voxels(l_dim.x * l_dim.y * l_dim.z);
  std::vector<std::vector<float>> brickData(num_of_bricks_per_rank,
                                            std::vector<float>(l_num_voxels));
  for (size_t i = 0; i < num_of_bricks_per_rank; ++i) {
    size_t id = num_of_bricks_per_rank * world_rank + i;
    if (data_mode == DATA_RANK) {
      // volume type 0
      std::fill(brickData[i].begin(), brickData[i].end(), float(id));
    } else if (data_mode == DATA_SPHERE) {
      // volume type 1 sphere -> (0, world_size * num_of_bricks_per_rank - 1)
      const vec3f offset((l_dim.x - 3.f) * i - 1.f,
                         (l_dim.y - 3.f) * world_rank - 1.f,
                         0.f);
      for (size_t x = 0; x < l_dim.x; ++x)
        for (size_t y = 0; y < l_dim.y; ++y)
          for (size_t z = 0; z < l_dim.z; ++z) {
            const float d = length(offset + vec3f(x, y, z) - vec3f(128.f))
                / cbrt(255.f * 255.f * 255.f)
                * static_cast<float>(world_size * num_of_bricks_per_rank - 1);
            const size_t j = z * l_dim.x * l_dim.y + y * l_dim.x + x;
            brickData[i][j] = ospcommon::pow(d, 1.2f);
          }
    }
  }


  ///-- Create Renderer ---------------------------------------------------////
#if USE_MPI
  auto ospRen = ospNewRenderer("mpi_raycast");
#else
  auto ospRen = ospNewRenderer("scivis");
#endif


  ///-- Create Transfer Function ------------------------------------------////
  const std::vector<vec3f> colors = {
      vec3f(0, 0, 0.5),
      vec3f(0, 0, 1),
      vec3f(0, 1, 1),
      vec3f(0.5, 1, 0.5),
      vec3f(1, 1, 0),
      vec3f(1, 0, 0),
      vec3f(0.5, 0, 0),
  };
  const std::vector<float> opacities = {
      0.1f,
      0.01f,
      0.1f,
      0.01f,
      0.1f,
      0.01f,
  };
  auto c_data = ospNewData(colors.size(), OSP_FLOAT3, colors.data());
  ospCommit(c_data);
  auto o_data = ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
  ospCommit(o_data);
  auto ospTfn = ospNewTransferFunction("piecewise_linear");
  ospSetData(ospTfn, "colors", c_data);
  ospSetData(ospTfn, "opacities", o_data);
  ospSet2f(ospTfn, "valueRange", 0, world_size * num_of_bricks_per_rank - 1);
  ospCommit(ospTfn);
  ospRelease(c_data);
  ospRelease(o_data);


  ///-- Create Volume & Models --------------------------------------------////
#if USE_MPI
  std::vector<OSPModel> models;
#else
  auto ospMod = ospNewModel();
#endif
  //> loop over all bricks for a rank
  for (size_t i = 0; i < num_of_bricks_per_rank; ++i) {

    //> bounding boxes
    const vec3f l_dim_f(l_dim.x, l_dim.y, l_dim.z);
    const vec3f actual_bd((l_dim.x - 3.f) * i - 1.f,
                          (l_dim.y - 3.f) * world_rank - 1.f,
                          0.f);
    const box3f logical_bd(actual_bd + 1.f, actual_bd + l_dim_f - 2.f);
    const int id = static_cast<int>(world_rank * num_of_bricks_per_rank + i);

    //> debug
    for (int r = 0; r < world_size; ++r) {
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      if (r == world_rank) {
        std::cerr << "rank " << world_rank << " id " << id << " dims "
                  << l_dim << std::endl
                  << " actual  "
                  << actual_bd << " "
                  << " logical "
                  << logical_bd.lower << " "
                  << logical_bd.upper << " "
                  << std::endl;
      }
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    //> volume configuration
    auto ospVol = ospNewVolume("shared_structured_volume");
    auto ospVox = ospNewData(brickData[i].size(), OSP_FLOAT,
                             brickData[i].data(), OSP_DATA_SHARED_BUFFER);
    ospSetString(ospVol, "voxelType", "float");
    ospSet3i(ospVol, "dimensions",
             static_cast<int>(l_dim.x),
             static_cast<int>(l_dim.y),
             static_cast<int>(l_dim.z));
    ospSet3f(ospVol, "gridOrigin", actual_bd.x, actual_bd.y, actual_bd.z);
    ospSet3f(ospVol, "gridSpacing", 1.0f, 1.0f, 1.0f);
    ospSet1f(ospVol, "samplingRate", 0.2f);
    ospSet1i(ospVol, "gradientShadingEnabled", 0);
    ospSet1i(ospVol, "adaptiveSampling", 0);
    ospSet1i(ospVol, "preIntegration", 0);
    ospSet1i(ospVol, "singleShade", 0);
#if !USE_MPI
    ospSetVec3f(ospVol, "volumeClippingBoxLower",
            (const osp::vec3f&)(logical_bd.lower));
    ospSetVec3f(ospVol, "volumeClippingBoxUpper",
            (const osp::vec3f&)(logical_bd.upper));
#endif
    ospSetData(ospVol, "voxelData", ospVox);
    ospSetObject(ospVol, "transferFunction", ospTfn);
    ospCommit(ospVol);
    ospRelease(ospVox);

    //> create model
#if USE_MPI
    auto ospMod = ospNewModel();
#endif

    //> add to model
    ospAddVolume(ospMod, ospVol);

    //> configure model
#if USE_MPI
    ospSet1i(ospMod, "id", id);
    ospSetVec3f(ospMod, "region.lower", (const osp::vec3f &) (logical_bd.lower));
    ospSetVec3f(ospMod, "region.upper", (const osp::vec3f &) (logical_bd.upper));
    ospCommit(ospMod);
    models.push_back(ospMod);
#endif

  }


  ///-- Setup Model -------------------------------------------------------////
#if USE_MPI
  auto mData = ospNewData(models.size(), OSP_OBJECT, models.data(), OSP_DATA_SHARED_BUFFER);
  ospSetData(ospRen, "models", mData);
#else
  ospCommit(ospMod);
  ospSetObject(ospRen, "model", ospMod);
  ospSetObject(ospRen, "world", ospMod);
#endif


  ///-- Create Camera -----------------------------------------------------////
  auto ospCam = ospNewCamera("perspective");
  ospSet3f(ospCam, "dir", vi.x - vp.x, vi.y - vp.y, vi.z - vp.z);
  ospSet3f(ospCam, "pos", vp.x, vp.y, vp.z);
  ospSet3f(ospCam, "up", vu.x, vu.y, vu.z);
  ospSet1f(ospCam, "aspect", static_cast<float>(fb_size.x) / fb_size.y);
  ospCommit(ospCam);
  ospSetObject(ospRen, "camera", ospCam);


  ///-- Configure Renderer ------------------------------------------------////
  const float bgColor = (1 + world_rank) / (float) world_size;
  ospSet3f(ospRen, "bgColor", bgColor, bgColor, bgColor);
  ospSet1i(ospRen, "shadowEnabled", 0);
  ospSet1i(ospRen, "oneSidedLighting", 0);
  ospCommit(ospRen);


  ///-- Testing Pixel Operation -------------------------------------------////
#if USE_PIXEL_OP //@ Setup pixel operation

  //> valid pixel operation color format
  const auto ospPixelOpFormat = OSP_FB_SRGBA;
//const auto ospPixelOpFormat = OSP_FB_RGBA8;
//const auto ospPixelOpFormat = OSP_FB_RGBA32F;

  //> generate display configuration: configuration XxY
  const vec2i display_num{num_of_displays, num_of_image_per_display};
  const vec2i display_dim = fb_size / display_num;

  //> open a MPI window
  MPI_Win the_window;
  const size_t display_nelem = static_cast<size_t>(display_dim.x * display_dim.y) * 4;
  const size_t display_unit = ospPixelOpFormat == OSP_FB_RGBA32F ? sizeof(float) : 1;
  std::vector<uint8_t> display_buf(display_unit * display_nelem * num_of_image_per_display);
  MPI_Win_create(display_buf.data(), display_buf.size(),
                 1 /* make sure disp_unit is always 1 */,
                 MPI_INFO_NULL, MPI_COMM_WORLD, &the_window);

  //> configure the pixel operation
  auto ospPixelOp = ospNewPixelOp("distributed_streamer");
  std::vector<DisplayInfo> display_info(display_num.x * display_num.y);
  for (auto dy = 0; dy < display_num.y; ++dy) {
    for (auto dx = 0; dx < display_num.x; ++dx) {
      const auto d = dx + dy * display_num.x;
      display_info[d].rank = d / num_of_image_per_display;
      display_info[d].byte_offset = display_unit * display_nelem * (d % num_of_image_per_display);
      display_info[d].tile_offset = vec2i(dx, dy);
    }
  }
  ospSetData(ospPixelOp, "displayConfigurations",
             ospNewData(display_info.size() * sizeof(DisplayInfo), OSP_UCHAR, display_info.data()));
  ospSetVec2i(ospPixelOp, "displayDimension", (const osp::vec2i &) display_num);
  ospSetVec2i(ospPixelOp, "displayTileSize",  (const osp::vec2i &) display_dim);
  ospSet1i(ospPixelOp, "displayFormat", ospPixelOpFormat);
  ospSetVoidPtr(ospPixelOp, "windowObject", &the_window);
  ospCommit(ospPixelOp);

#endif


  ///-- Create FrameBuffer ------------------------------------------------////

  //> valid ospray framebuffer format
  const auto ospFbFormat = OSP_FB_NONE;
//const auto ospFbFormat = OSP_FB_SRGBA;
//const auto ospFbFormat = OSP_FB_RGBA8;
//const auto ospFbFormat = OSP_FB_RGBA32F;

  //> valid ospray framebuffer channel
  const auto ospFbChannel = OSP_FB_COLOR | OSP_FB_ACCUM | OSP_FB_VARIANCE;

  //> create framebuffer 
  auto ospFb = ospNewFrameBuffer((const osp::vec2i &) fb_size, ospFbFormat, ospFbChannel);
#if USE_PIXEL_OP
  ospSetPixelOp(ospFb, ospPixelOp);
#endif
  ospCommit(ospFb);
  ospFrameBufferClear(ospFb, ospFbChannel);

  
  ///-- Rendering ---------------------------------------------------------////
  {
    double t1, t2;
    t1 = MPI_Wtime();
    for (int32_t f = 0; f < num_of_frames; ++f) {
      ospRenderFrame(ospFb, ospRen, ospFbChannel);
    }
    t2 = MPI_Wtime();
    printf("Frame rate is %f\n", num_of_frames / (t2 - t1));
  }


  ///-- Save Images -------------------------------------------------------////
  //> normal framebuffer
  if (world_rank == 0 && ospFbFormat != OSP_FB_NONE) {
    const void *mapped = ospMapFrameBuffer(ospFb, OSP_FB_COLOR);
    switch (ospFbFormat) {
    case OSP_FB_SRGBA:
    case OSP_FB_RGBA8: debug(0) << "save image (char)" << std::endl;
      utility::writePPM("framebuffer.ppm", fb_size.x, fb_size.y, (const uint32_t *) mapped);
      break;
    case OSP_FB_RGBA32F: debug(0) << "save image (float)" << std::endl;
      utility::writePPM("framebuffer.ppm", fb_size.x, fb_size.y, (const float *) mapped);
      break;
    default: debug(0) << "uncaught framebuffer format " << ospFbFormat << std::endl;
      break;
    }
    ospUnmapFrameBuffer(mapped, ospFb);
  }    
  //> pixelop tiles
#if USE_PIXEL_OP // Receive images from pixelop
  MPI_Win_free(&the_window); // close mpi window
  if (world_rank < num_of_displays) {
    for (auto i = 0; i < num_of_image_per_display; ++i) {
      const std::string fname = "tile" + std::to_string(world_rank) + "-" + std::to_string(i) + ".ppm";
      switch (ospPixelOpFormat) {
      case OSP_FB_SRGBA:
      case OSP_FB_RGBA8: debug(0) << "save tile (char)" << std::endl;
	utility::writePPM(fname, display_dim.x, display_dim.y, reinterpret_cast<uint32_t*>
			  (display_buf.data() + i * display_nelem * display_unit));
	break;
      case OSP_FB_RGBA32F: debug(0) << "save tile (float)" << std::endl;
	utility::writePPM(fname , display_dim.x, display_dim.y, reinterpret_cast<float *>
			  (display_buf.data() + i * display_nelem * display_unit));
	break;
      default: debug(0) << "uncaught framebuffer format " << ospFbFormat << std::endl;
	break;
      }
    }
  }

#endif


  ///-- Close Program -----------------------------------------------------////
  ospShutdown();
#if USE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

////////////////////
