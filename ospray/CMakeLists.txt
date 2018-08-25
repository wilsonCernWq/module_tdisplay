## ======================================================================== ##
## Copyright 2009-2017 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

OSPRAY_CONFIGURE_MPI()

OSPRAY_BUILD_COMPONENT(mpiCommon)
OSPRAY_BUILD_COMPONENT(mpiMessageLayer)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ospray)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/modules/mpi)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/ospray)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES_ISPC(${CMAKE_CURRENT_BINARY_DIR})

OSPRAY_CREATE_LIBRARY(ospray_module_tdisplay
  fb/DistributedStreamerPixelOp.h
  fb/DistributedStreamerPixelOp.cpp
  fb/DistributedStreamerPixelOp.ispc
  moduleInit.cpp
  moduleInit.h 
  LINK
  ospray
  ospray_common
  ospray_module_ispc
  ospray_mpi_common
  ospray_mpi_maml)
