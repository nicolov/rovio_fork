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

#ifndef ROVIO_POSEUPDATE_HPP_
#define ROVIO_POSEUPDATE_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/Update.hpp"
#include "lightweight_filtering/State.hpp"

#include "rovio_states.hpp"

namespace rovio {

struct PoseUpdateMeas : public PoseUpdateMeasManifold {
  inline Eigen::Matrix<double,6,6>& measuredCov() {
    return measuredCov_;
  }

  inline const Eigen::Matrix<double,6,6>& measuredCov() const {
    return measuredCov_;
  }
private:
  Eigen::Matrix<double,6,6> measuredCov_;
};

class PoseOutlierDetection: public LWF::OutlierDetection<LWF::ODEntry<PoseInnovation::pos_idx_,6>>{
 public:
  virtual ~PoseOutlierDetection(){};
};

/** \brief Class for adding poseUpdates to the filter.
 *
 * Coordinate frames overview:
 * I: Inertial frame of measured pose
 * V: Body frame of measured pose
 * W: Inertial frame of ROVIO
 * M: IMU-coordinate frame
 *
 *  @tparam FILTERSTATE         - FilterState
 *  @tparam inertialPoseIndex   - Index where the estimated relative pose between inertial frames is stored. Set -1 if this should not be estimated.
 *  @tparam bodyPoseIndex       - Index where the estimated relative pose between body frames is stored. Set -1 if this should not be estimated.
 */

static constexpr int inertialPoseIndex = -1;
static constexpr int bodyPoseIndex = -1;

class PoseUpdate: public LWF::Update<PoseInnovation,rovio::FilterState,PoseUpdateMeas,PoseUpdateNoise,PoseOutlierDetection,false>{
 public:
  typedef LWF::Update<PoseInnovation,rovio::FilterState,PoseUpdateMeas,PoseUpdateNoise,PoseOutlierDetection,false> Base;
  using Base::eval;
  using Base::boolRegister_;
  using Base::intRegister_;
  using Base::doubleRegister_;
  using Base::meas_;
  using Base::updnoiP_; // This is the update covariance as used by the Kalman functions.
  typedef typename Base::mtState mtState;
  typedef typename Base::mtFilterState mtFilterState;
  typedef typename Base::mtInnovation mtInnovation;
  typedef typename Base::mtMeas mtMeas;
  typedef typename Base::mtNoise mtNoise;
  typedef typename Base::mtOutlierDetection mtOutlierDetection;
  static constexpr int inertialPoseIndex_ = inertialPoseIndex;
  static constexpr int bodyPoseIndex_ = bodyPoseIndex;
  QPD qVM_;
  V3D MrMV_;
  QPD qWI_;
  V3D IrIW_;
  Eigen::MatrixXd defaultUpdnoiP_; // Configured update covariance, that will (optionally) be scaled by the measurement
  double timeOffset_;
  bool enablePosition_;
  bool enableAttitude_;
  bool noFeedbackToRovio_;
  bool doInertialAlignmentAtStart_;
  bool didAlignment_;
  bool useOdometryCov_;
  PoseUpdate() : defaultUpdnoiP_((int)(mtNoise::D_),(int)(mtNoise::D_)) {
    static_assert(mtState::nPose_>inertialPoseIndex_,"Please add enough poses to the filter state (templated).");
    static_assert(mtState::nPose_>bodyPoseIndex_,"Please add enough poses to the filter state (templated).");
    qVM_.setIdentity();
    MrMV_.setZero();
    qWI_.setIdentity();
    IrIW_.setZero();
    timeOffset_ = 0.0;
    enablePosition_ = true;
    enableAttitude_ = true;
    noFeedbackToRovio_ = true;
    doInertialAlignmentAtStart_ = true;
    didAlignment_ = false;
    useOdometryCov_ = false;
    doubleRegister_.registerVector("MrMV",MrMV_);
    doubleRegister_.registerQuaternion("qVM",qVM_);
    doubleRegister_.registerVector("IrIW",IrIW_);
    doubleRegister_.registerQuaternion("qWI",qWI_);
    intRegister_.removeScalarByStr("maxNumIteration");
    doubleRegister_.removeScalarByStr("alpha");
    doubleRegister_.removeScalarByStr("beta");
    doubleRegister_.removeScalarByStr("kappa");
    doubleRegister_.removeScalarByStr("updateVecNormTermination");
    doubleRegister_.registerScalar("timeOffset",timeOffset_);
    boolRegister_.registerScalar("enablePosition",enablePosition_);
    boolRegister_.registerScalar("enableAttitude",enableAttitude_);
    boolRegister_.registerScalar("noFeedbackToRovio",noFeedbackToRovio_);
    boolRegister_.registerScalar("doInertialAlignmentAtStart",doInertialAlignmentAtStart_);
    boolRegister_.registerScalar("useOdometryCov",useOdometryCov_);

    // Unregister configured covariance
    for (int i=0;i<6;i++) {
      doubleRegister_.removeScalarByVar(updnoiP_(i,i));
    }
    // Register configured covariance again under a different name
    mtNoise n;
    n.setIdentity();
    n.registerCovarianceToPropertyHandler_(defaultUpdnoiP_,this,"UpdateNoise.");
  }
  virtual ~PoseUpdate(){}
  const V3D& get_IrIW(const mtState& state) const{
    if(inertialPoseIndex_ >= 0){
      return state.poseLin(inertialPoseIndex_);
    } else {
      return IrIW_;
    }
  }
  const QPD& get_qWI(const mtState& state) const{
    if(inertialPoseIndex_ >= 0){
      return state.poseRot(inertialPoseIndex_);
    } else {
      return qWI_;
    }
  }
  const V3D& get_MrMV(const mtState& state) const{
    if(bodyPoseIndex_ >= 0){
      return state.poseLin(bodyPoseIndex_);
    } else {
      return MrMV_;
    }
  }
  const QPD& get_qVM(const mtState& state) const{
    if(bodyPoseIndex_ >= 0){
      return state.poseRot(bodyPoseIndex_);
    } else {
      return qVM_;
    }
  }
  void evalInnovation(mtInnovation& y, const mtState& state, const mtNoise& noise) const{
    // IrIV = IrIW + qWI^T*(WrWM + qWM*MrMV)
    // qVI = qVM*qWM^T*qWI
    if(enablePosition_){
      y.pos() = get_IrIW(state) + get_qWI(state).inverse()*V3D(state.WrWM()+state.qWM() * get_MrMV(state)) - meas_.pos() + noise.pos();
    } else {
      y.pos() = noise.pos();
    }
    if(enableAttitude_){
      QPD attNoise = exp_map(noise.att());
      y.att() = attNoise*get_qVM(state)*state.qWM().inverse()*get_qWI(state)*meas_.att().inverse();
    } else {
      QPD attNoise = exp_map(noise.att());
      y.att() = attNoise;
    }
  }
  void jacState(MXD& F, const mtState& state) const{
    F.setZero();
    if(enablePosition_){
      if(!noFeedbackToRovio_){
        F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) =
            MPD(get_qWI(state).inverse()).matrix();
        F.template block<3,3>(mtInnovation::pos_idx_,mtState::att_idx_) =
            -MPD(get_qWI(state).inverse()).matrix()*gSM(state.qWM() * get_MrMV(state));
      }
      if(inertialPoseIndex_ >= 0){
        F.template block<3,3>(mtInnovation::pos_idx_,mtState::pop_idx(inertialPoseIndex_)) =
            M3D::Identity();
        F.template block<3,3>(mtInnovation::pos_idx_,mtState::poa_idx(inertialPoseIndex_)) =
            gSM(get_qWI(state).inverse()*V3D(state.WrWM()+state.qWM()*get_MrMV(state)))*MPD(get_qWI(state).inverse()).matrix();
      }
      if(bodyPoseIndex_ >= 0){
        F.template block<3,3>(mtInnovation::pos_idx_,mtState::pop_idx(bodyPoseIndex_)) =
            MPD(get_qWI(state).inverse()*state.qWM()).matrix();
      }
    }
    if(enableAttitude_){
      if(!noFeedbackToRovio_){
        F.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) =
            -MPD(get_qVM(state)*state.qWM().inverse()).matrix();
      }
      if(inertialPoseIndex_ >= 0){
        F.template block<3,3>(mtInnovation::att_idx_,mtState::poa_idx(inertialPoseIndex_)) =
            MPD(get_qVM(state)*state.qWM().inverse()).matrix();
      }
      if(bodyPoseIndex_ >= 0){
        F.template block<3,3>(mtInnovation::att_idx_,mtState::poa_idx(bodyPoseIndex_)) =
            M3D::Identity();
      }
    }
  }
  void jacNoise(MXD& G, const mtState& state) const{
    G.setZero();
    G.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    G.template block<3,3>(mtInnovation::att_idx_,mtNoise::att_idx_) = M3D::Identity();
  }
  void preProcess(mtFilterState& filterstate, const mtMeas& meas, bool& isFinished){
    mtState& state = filterstate.state_;
    isFinished = false;
    if(!didAlignment_ && doInertialAlignmentAtStart_){
      // qWI = qWM*qVM^T*qVI;
      qWI_ = state.qWM()*get_qVM(state).inverse()*meas.att();
      if(inertialPoseIndex_ >= 0){
        state.poseRot(inertialPoseIndex_) = qWI_;
      }
      // IrIW = IrIV - qWI^T*(WrWM + qWM*MrMV);
      IrIW_ = meas.pos() - qWI_.inverse()*V3D(state.WrWM() + state.qWM() * get_MrMV(state));
      if(inertialPoseIndex_ >= 0){
        state.poseLin(inertialPoseIndex_) = IrIW_;
      }
      didAlignment_ = true;
    }

    // When enabled, scale the configured position covariance by the values in the measurement
    if(useOdometryCov_){
      updnoiP_ = defaultUpdnoiP_;
      updnoiP_ *= meas.measuredCov();

      // When either position or attitude are disabled, we need to make sure that the covariance matrix is block-diagonal,
      // otherwise the unused covariance block would affect the used one when inverting later on.
      if (enablePosition_ != enableAttitude_) {
        updnoiP_.template block<3,3>(mtInnovation::att_idx_, mtInnovation::pos_idx_).setZero();
        updnoiP_.template block<3,3>(mtInnovation::pos_idx_, mtInnovation::att_idx_).setZero();
      }
    } else {
      updnoiP_ = defaultUpdnoiP_;
    }
    /* std::cout << "Default\n" << defaultUpdnoiP_ << "\n\n"
              << "Meas\n" << meas.measuredCov() << "\n\n"
              << "Scaled (" << useOdometryCov_ << ")\n" << updnoiP_ << "\n\n"; */
  }
  void postProcess(mtFilterState& filterstate, const mtMeas& meas, const mtOutlierDetection& outlierDetection, bool& isFinished){
    mtState& state = filterstate.state_;
    isFinished = true;
    // WrWC = qWI*(IrIV - qWI^T*qWM*MrMV -IrIW) +qWM*MrMC
    state.aux().poseMeasLin_ = get_qWI(state)*V3D(meas.pos()-(get_qWI(state).inverse()*state.qWM())*get_MrMV(state)-get_IrIW(state))+state.att()*state.MrMC(0);
    // qCW = qCM*qVM^T*qVI*qWI^T;
    state.aux().poseMeasRot_ = state.qCM(0)*get_qVM(state).inverse()*meas.att()*get_qWI(state).inverse();
  }
};

}


#endif /* ROVIO_POSEUPDATE_HPP_ */
