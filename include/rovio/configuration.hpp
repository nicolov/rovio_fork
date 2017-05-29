#pragma once

#include "lightweight_filtering/State.hpp"
#include "rovio/ImagePyramid.hpp"

namespace rovio {

namespace conf {
  static constexpr int nMax = ROVIO_NMAXFEATURE; // Maximal number of considered features in the filter state.
  static constexpr int nLevels = ROVIO_NLEVELS; // // Total number of pyramid levels considered.
  static constexpr int patchSize = ROVIO_PATCHSIZE; // Edge length of the patches (in pixel). Must be a multiple of 2!
  static constexpr int nCam = ROVIO_NCAM; // Used total number of cameras.
  static constexpr int nPose = ROVIO_NPOSE; // Additional pose states.
}

}

