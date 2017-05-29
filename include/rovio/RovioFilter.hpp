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

#ifndef ROVIO_ROVIO_FILTER_HPP_
#define ROVIO_ROVIO_FILTER_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/FilterBase.hpp"
#include "rovio/FilterStates.hpp"
#include "rovio/ImgUpdate.hpp"
#include "rovio/PoseUpdate.hpp"
#include "rovio/VelocityUpdate.hpp"
#include "rovio/ImuPrediction.hpp"
#include "rovio/MultiCamera.hpp"

namespace rovio {

class RovioFilter:public LWF::FilterBase<ImuPrediction,ImgUpdate,PoseUpdate,VelocityUpdate>{
 public:
  typedef LWF::FilterBase<ImuPrediction,ImgUpdate,PoseUpdate,VelocityUpdate> Base;
  using Base::init_;
  using Base::reset;
  using Base::predictionTimeline_;
  using Base::safe_;
  using Base::front_;
  using Base::readFromInfo;
  using Base::boolRegister_;
  using Base::intRegister_;
  using Base::doubleRegister_;
  using Base::mUpdates_;
  using Base::mPrediction_;
  using Base::stringRegister_;
  using Base::subHandlers_;
  using Base::updateToUpdateMeasOnly_;
  typedef typename Base::mtFilterState mtFilterState;
  typedef typename Base::mtPrediction mtPrediction;
  typedef typename Base::mtState mtState;
  rovio::MultiCamera<mtState::nCam_> multiCamera_;
  std::string cameraCalibrationFile_[mtState::nCam_];
  int depthTypeInt_;

  /** \brief Constructor. Initializes the filter.
   */
  RovioFilter(){
    updateToUpdateMeasOnly_ = true;
    std::get<0>(mUpdates_).setCamera(&multiCamera_);
    init_.setCamera(&multiCamera_);
    depthTypeInt_ = 1;
    subHandlers_.erase("Update0");
    subHandlers_["ImgUpdate"] = &std::get<0>(mUpdates_);
    subHandlers_.erase("Update1");
    subHandlers_["PoseUpdate"] = &std::get<1>(mUpdates_);
    subHandlers_.erase("Update2");
    subHandlers_["VelocityUpdate"] = &std::get<2>(mUpdates_);
    boolRegister_.registerScalar("Common.doVECalibration",init_.state_.aux().doVECalibration_);
    intRegister_.registerScalar("Common.depthType",depthTypeInt_);
    for(int camID=0;camID<mtState::nCam_;camID++){
      cameraCalibrationFile_[camID] = "";
      stringRegister_.registerScalar("Camera" + std::to_string(camID) + ".CalibrationFile",cameraCalibrationFile_[camID]);
      doubleRegister_.registerVector("Camera" + std::to_string(camID) + ".MrMC",init_.state_.aux().MrMC_[camID]);
      doubleRegister_.registerQuaternion("Camera" + std::to_string(camID) + ".qCM",init_.state_.aux().qCM_[camID]);
      doubleRegister_.removeScalarByVar(init_.state_.MrMC(camID)(0));
      doubleRegister_.removeScalarByVar(init_.state_.MrMC(camID)(1));
      doubleRegister_.removeScalarByVar(init_.state_.MrMC(camID)(2));
      doubleRegister_.removeScalarByVar(init_.state_.qCM(camID).w());
      doubleRegister_.removeScalarByVar(init_.state_.qCM(camID).x());
      doubleRegister_.removeScalarByVar(init_.state_.qCM(camID).y());
      doubleRegister_.removeScalarByVar(init_.state_.qCM(camID).z());
      for(int j=0;j<3;j++){
        doubleRegister_.removeScalarByVar(init_.cov_(mtState::vep_idx(camID)+j,mtState::vep_idx(camID)+j));
        doubleRegister_.removeScalarByVar(init_.cov_(mtState::vea_idx(camID)+j,mtState::vea_idx(camID)+j));
        doubleRegister_.registerScalar("Init.Covariance.vep",init_.cov_(mtState::vep_idx(camID)+j,mtState::vep_idx(camID)+j));
        doubleRegister_.registerScalar("Init.Covariance.vea",init_.cov_(mtState::vea_idx(camID)+j,mtState::vea_idx(camID)+j));
      }
      doubleRegister_.registerVector("Camera" + std::to_string(camID) + ".MrMC",init_.state_.MrMC(camID));
      doubleRegister_.registerQuaternion("Camera" + std::to_string(camID) + ".qCM",init_.state_.qCM(camID));
    }
    for(int i=0;i<mtState::nPose_;i++){
      doubleRegister_.removeScalarByVar(init_.state_.poseLin(i)(0));
      doubleRegister_.removeScalarByVar(init_.state_.poseLin(i)(1));
      doubleRegister_.removeScalarByVar(init_.state_.poseLin(i)(2));
      doubleRegister_.removeScalarByVar(init_.state_.poseRot(i).w());
      doubleRegister_.removeScalarByVar(init_.state_.poseRot(i).x());
      doubleRegister_.removeScalarByVar(init_.state_.poseRot(i).y());
      doubleRegister_.removeScalarByVar(init_.state_.poseRot(i).z());
      for(int j=0;j<3;j++){
        doubleRegister_.removeScalarByVar(init_.cov_(mtState::pop_idx(i)+j,mtState::pop_idx(i)+j));
        doubleRegister_.removeScalarByVar(init_.cov_(mtState::poa_idx(i)+j,mtState::poa_idx(i)+j));
      }
    }
    if(std::get<1>(mUpdates_).inertialPoseIndex_>=0){
      std::get<1>(mUpdates_).doubleRegister_.registerVector("IrIW",init_.state_.poseLin(std::get<1>(mUpdates_).inertialPoseIndex_));
      std::get<1>(mUpdates_).doubleRegister_.registerQuaternion("qWI",init_.state_.poseRot(std::get<1>(mUpdates_).inertialPoseIndex_));
      for(int j=0;j<3;j++){
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("init_cov_IrIW",init_.cov_(mtState::pop_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j,mtState::pop_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("init_cov_qWI",init_.cov_(mtState::poa_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j,mtState::poa_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("pre_cov_IrIW",mPrediction_.prenoiP_(mtPrediction::mtNoise::pop_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j,mtPrediction::mtNoise::pop_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("pre_cov_qWI",mPrediction_.prenoiP_(mtPrediction::mtNoise::poa_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j,mtPrediction::mtNoise::poa_idx(std::get<1>(mUpdates_).inertialPoseIndex_)+j));
      }
    }
    if(std::get<1>(mUpdates_).bodyPoseIndex_>=0){
      std::get<1>(mUpdates_).doubleRegister_.registerVector("MrMV",init_.state_.poseLin(std::get<1>(mUpdates_).bodyPoseIndex_));
      std::get<1>(mUpdates_).doubleRegister_.registerQuaternion("qVM",init_.state_.poseRot(std::get<1>(mUpdates_).bodyPoseIndex_));
      for(int j=0;j<3;j++){
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("init_cov_MrMV",init_.cov_(mtState::pop_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j,mtState::pop_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("init_cov_qVM",init_.cov_(mtState::poa_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j,mtState::poa_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("pre_cov_MrMV",mPrediction_.prenoiP_(mtPrediction::mtNoise::pop_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j,mtPrediction::mtNoise::pop_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j));
        std::get<1>(mUpdates_).doubleRegister_.registerScalar("pre_cov_qVM",mPrediction_.prenoiP_(mtPrediction::mtNoise::poa_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j,mtPrediction::mtNoise::poa_idx(std::get<1>(mUpdates_).bodyPoseIndex_)+j));
      }
    }
    int ind;
    for(int i=0;i<FilterState::mtState::nMax_;i++){
      ind = mtState::fea_idx(i);
      doubleRegister_.removeScalarByVar(init_.cov_(ind,ind));
      doubleRegister_.removeScalarByVar(init_.cov_(ind+1,ind+1));
      ind = mtState::fea_idx(i)+2;
      doubleRegister_.removeScalarByVar(init_.cov_(ind,ind));
      doubleRegister_.removeScalarByVar(init_.state_.dep(i).p_);
      doubleRegister_.removeScalarByVar(init_.state_.CfP(i).nor_.q_.w());
      doubleRegister_.removeScalarByVar(init_.state_.CfP(i).nor_.q_.x());
      doubleRegister_.removeScalarByVar(init_.state_.CfP(i).nor_.q_.y());
      doubleRegister_.removeScalarByVar(init_.state_.CfP(i).nor_.q_.z());
      std::get<0>(mUpdates_).intRegister_.registerScalar("statLocalQualityRange",init_.fsm_.features_[i].mpStatistics_->localQualityRange_);
      std::get<0>(mUpdates_).intRegister_.registerScalar("statLocalVisibilityRange",init_.fsm_.features_[i].mpStatistics_->localVisibilityRange_);
      std::get<0>(mUpdates_).intRegister_.registerScalar("statMinGlobalQualityRange",init_.fsm_.features_[i].mpStatistics_->minGlobalQualityRange_);
      std::get<0>(mUpdates_).boolRegister_.registerScalar("doPatchWarping",init_.state_.CfP(i).trackWarping_);
    }
    std::get<0>(mUpdates_).doubleRegister_.removeScalarByVar(std::get<0>(mUpdates_).outlierDetection_.getMahalTh(0));
    std::get<0>(mUpdates_).doubleRegister_.registerScalar("MahalanobisTh",std::get<0>(mUpdates_).outlierDetection_.getMahalTh(0));
    std::get<0>(mUpdates_).outlierDetection_.setEnabledAll(true);
    std::get<1>(mUpdates_).outlierDetection_.setEnabledAll(true);
    boolRegister_.registerScalar("Common.verbose",std::get<0>(mUpdates_).verbose_);
    mPrediction_.doubleRegister_.removeScalarByStr("alpha");
    mPrediction_.doubleRegister_.removeScalarByStr("beta");
    mPrediction_.doubleRegister_.removeScalarByStr("kappa");
    boolRegister_.registerScalar("PoseUpdate.doVisualization",init_.plotPoseMeas_);
    reset(0.0);
  }

  /** \brief Reloads the camera calibration for all cameras and resets the depth map type.
   */
  void refreshProperties(){
    if(std::get<0>(mUpdates_).useDirectMethod_){
      init_.mode_ = LWF::ModeIEKF;
    } else {
      init_.mode_ = LWF::ModeEKF;
    }
    for(int camID = 0;camID<mtState::nCam_;camID++){
      if (!cameraCalibrationFile_[camID].empty()) {
        multiCamera_.cameras_[camID].load(cameraCalibrationFile_[camID]);
      }
    }
    for(int i=0;i<FilterState::mtState::nMax_;i++){
      init_.state_.dep(i).setType(depthTypeInt_);
    }
  };

  /** \brief Destructor
   */
  virtual ~RovioFilter(){};
//  void resetToImuPose(V3D WrWM, QPD qMW, double t = 0.0){
//    init_.state_.initWithImuPose(WrWM,qMW);
//    reset(t);
//  }

  /** \brief Resets the filter with an accelerometer measurement.
   *
   *  @param fMeasInit - Accelerometer measurement.
   *  @param t         - Current time.
   */
  void resetWithAccelerometer(const V3D& fMeasInit, double t = 0.0){
    init_.initWithAccelerometer(fMeasInit);
    reset(t);
  }

  /** \brief Resets the filter with an external pose.
   *
   *  @param WrWM - Position Vector, pointing from the World-Frame to the IMU-Frame, expressed in World-Coordinates.
   *  @param qMW  - Quaternion, expressing World-Frame in IMU-Coordinates (World Coordinates->IMU Coordinates)
   *  @param t    - Current time.
   */
  void resetWithPose(V3D WrWM, QPD qMW, double t = 0.0) {
    init_.initWithImuPose(WrWM, qMW);
    reset(t);
  }

  /** \brief Sets the transformation between IMU and Camera.
   *
   *  @param R_VM  -  Rotation matrix, expressing the orientation of the IMU  in Camera Cooridinates (IMU Coordinates -> Camera Coordinates).
   *  @param VrVM  -  Vector, pointing from the camera frame to the IMU frame, expressed in IMU Coordinates.
   *  @param camID -  ID of the considered camera.
   */
  void setExtrinsics(const Eigen::Matrix3d& R_CM, const Eigen::Vector3d& CrCM, const int camID = 0){
    init_.state_.aux().qCM_[camID] = QPD(R_CM);
    init_.state_.aux().MrMC_[camID] = -(init_.state_.aux().qCM_[camID].inverse()*CrCM);
  }
};

}


#endif /* ROVIO_ROVIO_FILTER_HPP_ */
