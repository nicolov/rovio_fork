/*
 * Copyright (c) 2014, Autonomous Systems Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * * Neither the name of the Autonomous Systems Lab, ETH Zurich nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef ROVIO_FILTERSTATES_HPP_
#define ROVIO_FILTERSTATES_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/FilterState.hpp"
#include <map>
#include <unordered_set>
#include "CoordinateTransform/FeatureOutput.hpp"
#include "rovio/FeatureManager.hpp"
#include "rovio/MultiCamera.hpp"
#include "rovio/State.hpp"

#include "rovio_states.hpp"

namespace rovio {

/** \brief Class defining the overall filter state (state, prediction measurement, noise).
 *
 *  @tparam nMax      - Maximal number of considered features in the filter state.
 *  @tparam nLevels   - Total number of pyramid levels considered.
 *  @tparam patchSize - Edge length of the patches (in pixel). Must be a multiple of 2!
 *  @tparam nCam      - Used total number of cameras.
 *  @tparam nPose     - Additional 6D pose in state.
 */
class FilterState: public LWF::FilterState<State,PredictionMeas,PredictionNoise,0>{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef LWF::FilterState<State,PredictionMeas,PredictionNoise,0> Base;
  typedef typename Base::mtState mtState;  /**<Local Filter %State Type. \see LWF::FilterState*/
  using Base::state_;  /**<Filter State. \see LWF::FilterState*/
  using Base::cov_;  /**<Filter State Covariance Matrix. \see LWF::FilterState*/
  using Base::usePredictionMerge_;  /**<Whether multiple subsequent prediction steps should be merged into one.*/
  FeatureSetManager<conf::nLevels,conf::patchSize,conf::nCam,conf::nMax> fsm_;
  mutable TransformFeatureOutputCT transformFeatureOutputCT_;
  mutable FeatureOutput featureOutput_;
  mutable MXD featureOutputCov_;
  cv::Mat img_[conf::nCam];     /**<Mainly used for drawing.*/
  cv::Mat patchDrawing_;  /**<Mainly used for drawing.*/
  int drawPB_;  /**<Size of border around patch.*/
  int drawPS_;  /**<Size of patch with border for drawing.*/
  double imgTime_;        /**<Time of the last image, which was processed.*/
  int imageCounter_;      /**<Total number of images, used so far for updates. Same as total number of update steps.*/
  ImagePyramid<conf::nLevels> prevPyr_[conf::nCam]; /**<Previous image pyramid.*/
  bool plotPoseMeas_; /**<Should the pose measurement be plotted.*/
  mutable MultilevelPatch<conf::nLevels,conf::patchSize> mlpErrorLog_[conf::nMax];  /**<Multilevel patch containing log of error.*/

  /** \brief Constructor
   */
  FilterState():fsm_(nullptr), transformFeatureOutputCT_(nullptr), featureOutputCov_((int)(FeatureOutput::D_),(int)(FeatureOutput::D_)){
    usePredictionMerge_ = true;
    imgTime_ = 0.0;
    imageCounter_ = 0;
    plotPoseMeas_ = true;
    state_.initFeatureManagers(fsm_);
    fsm_.allocateMissing();
    drawPB_ = 1;
    drawPS_ = mtState::patchSize_*pow(2,mtState::nLevels_-1)+2*drawPB_;
  }

  /** \brief Destructor
   */
  virtual ~FilterState(){};

  /** \brief Sets the multicamera pointer
   *
   * @param mpMultiCamera - multicamera pointer;
   */
  void setCamera(MultiCamera<conf::nCam>* mpMultiCamera){
    fsm_.setCamera(mpMultiCamera);
    transformFeatureOutputCT_.mpMultiCamera_ = mpMultiCamera;
  }

  /** \brief Initializes the FilterState \ref Base::state_ with the IMU-Pose.
   *
   *  @param WrWM - Position Vector, pointing from the World-Frame to the IMU-Frame, expressed in World-Coordinates.
   *  @param qMW  - Quaternion, expressing World-Frame in IMU-Coordinates (World Coordinates->IMU Coordinates)
   */
  void initWithImuPose(V3D WrWM, QPD qMW){
    state_.WrWM() = WrWM;
    state_.qWM()  = qMW.inverse();
  }

  /** \brief Initializes the FilterState \ref Base::state_ with the Acceleration-Vector.
   *
   *  @param fMeasInit - Acceleration-Vector which should be used for initializing the attitude of the IMU.
   */
  void initWithAccelerometer(const V3D& fMeasInit){
    V3D unitZ(0,0,1);
    if(fMeasInit.norm()>1e-6){
      state_.qWM().setFromTwoVectors(fMeasInit,unitZ);
    } else {
      state_.qWM().setIdentity();
    }
  }

  /** \brief Resets the covariance of a feature
   *
   *  @param i       - Feature index.
   *  @param initCov - Initialization 3x3 Covariance-Matrix.
   */
  void resetFeatureCovariance(unsigned int i,const Eigen::Matrix<double,3,3>& initCov){
    cov_.template block<mtState::D_,1>(0,mtState::fea_idx(i)+2).setZero();
    cov_.template block<1,mtState::D_>(mtState::fea_idx(i)+2,0).setZero();
    cov_.template block<mtState::D_,2>(0,mtState::fea_idx(i)).setZero();
    cov_.template block<2,mtState::D_>(mtState::fea_idx(i),0).setZero();
    cov_.template block<1,1>(mtState::fea_idx(i)+2,mtState::fea_idx(i)+2) = initCov.block<1,1>(0,0);
    cov_.template block<1,2>(mtState::fea_idx(i)+2,mtState::fea_idx(i)) = initCov.block<1,2>(0,1);
    cov_.template block<2,1>(mtState::fea_idx(i),mtState::fea_idx(i)+2) = initCov.block<2,1>(1,0);
    cov_.template block<2,2>(mtState::fea_idx(i),mtState::fea_idx(i)) = initCov.block<2,2>(1,1);
  }

  /** \brief Get the median distance parameter values of the state features for each camera.
   *
   *  \note The distance parameter type depends on the set \ref DepthType.
   *  @param initDistanceParameter     - Depth parameter value which is set, if no median distance parameter could be
   *                                  computed from the state features for a specific camera.
   *  @param medianDistanceParameters  - Array, containing the median distance parameter values for each camera.
   *  @param maxUncertaintyToDistanceRatio  - Maximal uncertainty where feature gets considered
   * */
  void getMedianDepthParameters(double initDistanceParameter, std::array<double,conf::nCam>* medianDistanceParameters, const float maxUncertaintyToDistanceRatio) {
    // Fill array with initialization value.
    // The initialization value is set, if no median distance value can be computed for a given camera frame.
    medianDistanceParameters->fill(initDistanceParameter);
    // Collect the distance values of the features for each camera frame.
    std::vector<double> distanceParameterCollection[conf::nCam];
    for (unsigned int i = 0; i < conf::nMax; i++) {
      if (fsm_.isValid_[i]) {
        for(int camID = 0;camID<conf::nCam;camID++){
          transformFeatureOutputCT_.setFeatureID(i);
          transformFeatureOutputCT_.setOutputCameraID(camID);
          transformFeatureOutputCT_.transformState(state_, featureOutput_);
          if(featureOutput_.c().isInFront()){
            transformFeatureOutputCT_.transformCovMat(state_, cov_, featureOutputCov_);
            const double uncertainty = std::fabs(sqrt(featureOutputCov_(2,2))*featureOutput_.d().getDistanceDerivative());
            const double depth = featureOutput_.d().getDistance();
            if(uncertainty/depth < maxUncertaintyToDistanceRatio){
              distanceParameterCollection[camID].push_back(featureOutput_.d().p_);
            }
          }
        }
      }
    }
    // Compute and store the median distance parameter.
    int size;
    for (unsigned int i = 0; i < conf::nCam; i++) {
      size = distanceParameterCollection[i].size();
      if(size > 3) { // Require a minimum of three features
        std::nth_element(distanceParameterCollection[i].begin(), distanceParameterCollection[i].begin() + size / 2, distanceParameterCollection[i].end());
        (*medianDistanceParameters)[i] = distanceParameterCollection[i][size/2];
      }
    }
  }
};

}


#endif /* ROVIO_FILTERSTATES_HPP_ */
