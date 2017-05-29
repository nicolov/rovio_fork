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

#ifndef ROVIO_IMUPREDICTION_HPP_
#define ROVIO_IMUPREDICTION_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/Prediction.hpp"
#include "lightweight_filtering/State.hpp"
#include "rovio/FilterStates.hpp"

namespace rovio {

class ImuPrediction: public LWF::Prediction<rovio::FilterState>{
 public:
  typedef LWF::Prediction<rovio::FilterState> Base;
  using Base::meas_;
  using Base::prenoiP_;
  using Base::doubleRegister_;
  using Base::intRegister_;
  using Base::boolRegister_;
  using Base::disablePreAndPostProcessingWarning_;
  typedef typename Base::mtState mtState;
  typedef typename Base::mtFilterState mtFilterState;
  typedef typename Base::mtMeas mtMeas;
  typedef typename Base::mtNoise mtNoise;
  const V3D g_; /**<Gravity in inertial frame, always aligned with the z-axis.*/
  double inertialMotionRorTh_; /**<Threshold on the rotational rate for motion detection.*/
  double inertialMotionAccTh_; /**<Threshold on the acceleration for motion detection.*/
  mutable FeatureCoordinates oldC_;
  mutable FeatureDistance oldD_;
  mutable Eigen::Matrix2d bearingVectorJac_;
  ImuPrediction():g_(0,0,-9.81){
    int ind;
    inertialMotionRorTh_ = 0.1;
    inertialMotionAccTh_ = 0.1;
    doubleRegister_.registerScalar("MotionDetection.inertialMotionRorTh",inertialMotionRorTh_);
    doubleRegister_.registerScalar("MotionDetection.inertialMotionAccTh",inertialMotionAccTh_);
    for(int i=0;i<mtState::nMax_;i++){
      ind = mtNoise::fea_idx(i);
      doubleRegister_.removeScalarByVar(prenoiP_(ind,ind));
      doubleRegister_.removeScalarByVar(prenoiP_(ind+1,ind+1));
      doubleRegister_.registerScalar("PredictionNoise.nor",prenoiP_(ind,ind));
      doubleRegister_.registerScalar("PredictionNoise.nor",prenoiP_(ind+1,ind+1));
      ind = mtNoise::fea_idx(i)+2;
      doubleRegister_.removeScalarByVar(prenoiP_(ind,ind));
      doubleRegister_.registerScalar("PredictionNoise.dep",prenoiP_(ind,ind));
    }
    for(int camID=0;camID<mtState::nCam_;camID++){
      for(int j=0;j<3;j++){
        doubleRegister_.removeScalarByVar(prenoiP_(mtNoise::vep_idx(camID)+j,mtNoise::vep_idx(camID)+j));
        doubleRegister_.removeScalarByVar(prenoiP_(mtNoise::vea_idx(camID)+j,mtNoise::vea_idx(camID)+j));
        doubleRegister_.registerScalar("PredictionNoise.vep",prenoiP_(mtNoise::vep_idx(camID)+j,mtNoise::vep_idx(camID)+j));
        doubleRegister_.registerScalar("PredictionNoise.vea",prenoiP_(mtNoise::vea_idx(camID)+j,mtNoise::vea_idx(camID)+j));
      }
    }
    for(int i=0;i<mtState::nPose_;i++){
      for(int j=0;j<3;j++){
        doubleRegister_.removeScalarByVar(prenoiP_(mtNoise::pop_idx(i)+j,mtNoise::pop_idx(i)+j));
        doubleRegister_.removeScalarByVar(prenoiP_(mtNoise::poa_idx(i)+j,mtNoise::poa_idx(i)+j));
      }
    }
    disablePreAndPostProcessingWarning_ = true;
  };

  /** \brief Destructor
   */
  virtual ~ImuPrediction(){};

  /* /brief Evaluation of prediction
   *
   * @todo implement without noise for speed
   */
  void evalPrediction(mtState& output, const mtState& state, const mtNoise& noise, double dt) const{
    output.aux().MwWMmeas_ = meas_.gyr();
    output.aux().MwWMest_  = meas_.gyr()-state.gyb();
    const V3D imuRor = output.aux().MwWMest_+noise.att()/sqrt(dt);
    const V3D dOmega = dt*imuRor;
    QPD dQ = exp_map(dOmega);
    for(unsigned int i=0;i<mtState::nMax_;i++){
      const int camID = state.CfP(i).camID_;
      if(&output != &state){
        output.CfP(i) = state.CfP(i);
        output.dep(i) = state.dep(i);
      }
      if(camID >= 0 && camID < mtState::nCam_){
        const V3D camRor = state.qCM(camID)*imuRor;
        const V3D camVel = state.qCM(camID)*V3D(imuRor.cross(state.MrMC(camID))-state.MvM());
        oldC_ = state.CfP(i);
        oldD_ = state.dep(i);
        output.dep(i).p_ = oldD_.p_-dt*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*camVel + noise.fea(i)(2)*sqrt(dt);
        V3D dm = -dt*(gSM(oldC_.get_nor().getVec())*camVel/oldD_.getDistance()
            + (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())*camRor)
            + oldC_.get_nor().getN()*noise.fea(i).template block<2,1>(0,0)*sqrt(dt);
        QPD qm = exp_map(dm);
        output.CfP(i).set_nor(oldC_.get_nor().rotated(qm));
        // WARP corners
        if(state.CfP(i).trackWarping_){
          bearingVectorJac_ = output.CfP(i).get_nor().getM().transpose()*(dt*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*(
                                  -1.0/oldD_.getDistance()*gSM(camVel)
                                  - (M3D::Identity()*(oldC_.get_nor().getVec().dot(camRor))+oldC_.get_nor().getVec()*camRor.transpose()))
                              +MPD(qm).matrix())*oldC_.get_nor().getM();
          output.CfP(i).set_warp_nor(bearingVectorJac_*oldC_.get_warp_nor());
        }
      }
    }
    output.WrWM() = state.WrWM()-dt*(state.qWM()*state.MvM()-noise.pos()/sqrt(dt));
    output.MvM() = (M3D::Identity()-gSM(dOmega))*state.MvM()-dt*(meas_.acc()-state.acb()+state.qWM().inverse()*g_-noise.vel()/sqrt(dt));
    output.acb() = state.acb()+noise.acb()*sqrt(dt);
    output.gyb() = state.gyb()+noise.gyb()*sqrt(dt);
    output.qWM() = state.qWM()*dQ;
    for(unsigned int i=0;i<mtState::nCam_;i++){
      output.MrMC(i) = state.MrMC(i)+noise.vep(i)*sqrt(dt);
      dQ = exp_map(noise.vea(i)*sqrt(dt));
      output.qCM(i) = dQ*state.qCM(i);
    }
    for(unsigned int i=0;i<mtState::nPose_;i++){
      output.poseLin(i) = state.poseLin(i)+noise.pop(i)*sqrt(dt);
      dQ = exp_map(noise.poa(i)*sqrt(dt));
      output.poseRot(i) = dQ*state.poseRot(i);
    }
    output.aux().wMeasCov_ = prenoiP_.template block<3,3>(mtNoise::att_idx_,mtNoise::att_idx_)/dt;
    output.fix();
    if(detectInertialMotion(state,meas_)){
      output.aux().timeSinceLastInertialMotion_ = 0;
    } else {
      output.aux().timeSinceLastInertialMotion_ = output.aux().timeSinceLastInertialMotion_ + dt;
    }
    output.aux().timeSinceLastImageMotion_ = output.aux().timeSinceLastImageMotion_ + dt;
  }
  void noMeasCase(mtFilterState& filterState, mtMeas& meas_, double dt){
    meas_.gyr() = filterState.state_.gyb();
    meas_.acc() = filterState.state_.acb()-filterState.state_.qWM().inverse()*g_;
  }
  void jacPreviousState(MXD& F, const mtState& state, double dt) const{
    const V3D imuRor = meas_.gyr()-state.gyb();
    const V3D dOmega = dt*imuRor;
    F.setZero();
    F.template block<3,3>(mtState::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    F.template block<3,3>(mtState::pos_idx_,mtState::vel_idx_) = -dt*MPD(state.qWM()).matrix();
    F.template block<3,3>(mtState::pos_idx_,mtState::att_idx_) = dt*gSM(state.qWM()*state.MvM());
    F.template block<3,3>(mtState::vel_idx_,mtState::vel_idx_) = (M3D::Identity()-gSM(dOmega));
    F.template block<3,3>(mtState::vel_idx_,mtState::acb_idx_) = dt*M3D::Identity();
    F.template block<3,3>(mtState::vel_idx_,mtState::gyb_idx_) = -dt*gSM(state.MvM());
    F.template block<3,3>(mtState::vel_idx_,mtState::att_idx_) = -dt*MPD(state.qWM()).matrix().transpose()*gSM(g_);
    F.template block<3,3>(mtState::acb_idx_,mtState::acb_idx_) = M3D::Identity();
    F.template block<3,3>(mtState::gyb_idx_,mtState::gyb_idx_) = M3D::Identity();
    F.template block<3,3>(mtState::att_idx_,mtState::gyb_idx_) = -dt*MPD(state.qWM()).matrix()*Lmat(dOmega);
    F.template block<3,3>(mtState::att_idx_,mtState::att_idx_) = M3D::Identity();
    LWF::NormalVectorElement nOut;
    QPD qm;
    for(unsigned int i=0;i<mtState::nMax_;i++){
      const int camID = state.CfP(i).camID_;
      if(camID >= 0 && camID < mtState::nCam_){
        const V3D camRor = state.qCM(camID)*imuRor;
        const V3D camVel = state.qCM(camID)*V3D(imuRor.cross(state.MrMC(camID))-state.MvM());
        oldC_ = state.CfP(i);
        oldD_ = state.dep(i);
        V3D dm = -dt*(gSM(oldC_.get_nor().getVec())*camVel/oldD_.getDistance()
            + (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())*camRor);
        qm = exp_map(dm);
        nOut = oldC_.get_nor().rotated(qm);
        F(mtState::fea_idx(i)+2,mtState::fea_idx(i)+2) = 1.0 - dt*oldD_.getParameterDerivativeCombined()
                *oldC_.get_nor().getVec().transpose()*camVel;
        F.template block<1,3>(mtState::fea_idx(i)+2,mtState::vel_idx_) =
            dt*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*MPD(state.qCM(camID)).matrix();
        F.template block<1,3>(mtState::fea_idx(i)+2,mtState::gyb_idx_) =
            -dt*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*gSM(state.qCM(camID)*state.MrMC(camID))*MPD(state.qCM(camID)).matrix();
        F.template block<1,2>(mtState::fea_idx(i)+2,mtState::fea_idx(i)) =
            -dt*oldD_.getParameterDerivative()*camVel.transpose()*oldC_.get_nor().getM();
        F.template block<2,2>(mtState::fea_idx(i),mtState::fea_idx(i)) =
            nOut.getM().transpose()*(
                    dt*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*(
                        -1.0/oldD_.getDistance()*gSM(camVel)
                        - (M3D::Identity()*(oldC_.get_nor().getVec().dot(camRor))+oldC_.get_nor().getVec()*camRor.transpose()))
                    +MPD(qm).matrix()
            )*oldC_.get_nor().getM();
        F.template block<2,1>(mtState::fea_idx(i),mtState::fea_idx(i)+2) =
            -nOut.getM().transpose()*gSM(nOut.getVec())*Lmat(dm)
                *dt*gSM(oldC_.get_nor().getVec())*camVel*(oldD_.getDistanceDerivative()/(oldD_.getDistance()*oldD_.getDistance()));
        F.template block<2,3>(mtState::fea_idx(i),mtState::vel_idx_) =
            -nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)
                *dt/oldD_.getDistance()*gSM(oldC_.get_nor().getVec())*MPD(state.qCM(camID)).matrix();
        F.template block<2,3>(mtState::fea_idx(i),mtState::gyb_idx_) =
            nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*(
                - (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())
                +1.0/oldD_.getDistance()*gSM(oldC_.get_nor().getVec())*gSM(state.qCM(camID)*state.MrMC(camID))
            )*dt*MPD(state.qCM(camID)).matrix();
        if(state.aux().doVECalibration_){
          F.template block<1,3>(mtState::fea_idx(i)+2,mtState::vea_idx(camID)) =
              dt*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*gSM(camVel);
          F.template block<1,3>(mtState::fea_idx(i)+2,mtState::vep_idx(camID)) =
              -dt*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*MPD(state.qCM(camID)).matrix()*gSM(imuRor);

          F.template block<2,3>(mtState::fea_idx(i),mtState::vea_idx(camID)) =
              -nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*(
                  (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())
              )*dt*gSM(camRor)
              -nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)
                  *dt/oldD_.getDistance()*gSM(oldC_.get_nor().getVec())*gSM(camVel);
          F.template block<2,3>(mtState::fea_idx(i),mtState::vep_idx(camID)) =
              nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)
                  *dt/oldD_.getDistance()*gSM(oldC_.get_nor().getVec())*MPD(state.qCM(camID)).matrix()*gSM(imuRor);
        }
      }
    }
    for(unsigned int i=0;i<mtState::nCam_;i++){
      F.template block<3,3>(mtState::vep_idx(i),mtState::vep_idx(i)) = M3D::Identity();
      F.template block<3,3>(mtState::vea_idx(i),mtState::vea_idx(i)) = M3D::Identity();
    }
    for(unsigned int i=0;i<mtState::nPose_;i++){
      F.template block<3,3>(mtState::pop_idx(i),mtState::pop_idx(i)) = M3D::Identity();
      F.template block<3,3>(mtState::poa_idx(i),mtState::poa_idx(i)) = M3D::Identity();
    }
  }
  void jacNoise(MXD& G, const mtState& state, double dt) const{
    const V3D imuRor = meas_.gyr()-state.gyb();
    const V3D dOmega = dt*imuRor;
    G.setZero();
    G.template block<3,3>(mtState::pos_idx_,mtNoise::pos_idx_) = M3D::Identity()*sqrt(dt);
    G.template block<3,3>(mtState::vel_idx_,mtNoise::vel_idx_) = M3D::Identity()*sqrt(dt);
    G.template block<3,3>(mtState::vel_idx_,mtNoise::att_idx_) = gSM(state.MvM())*sqrt(dt);
    G.template block<3,3>(mtState::acb_idx_,mtNoise::acb_idx_) = M3D::Identity()*sqrt(dt);
    G.template block<3,3>(mtState::gyb_idx_,mtNoise::gyb_idx_) = M3D::Identity()*sqrt(dt);
    G.template block<3,3>(mtState::att_idx_,mtNoise::att_idx_) = MPD(state.qWM()).matrix()*Lmat(dOmega)*sqrt(dt);
    for(unsigned int i=0;i<mtState::nCam_;i++){
      G.template block<3,3>(mtState::vep_idx(i),mtNoise::vep_idx(i)) = M3D::Identity()*sqrt(dt);
      G.template block<3,3>(mtState::vea_idx(i),mtNoise::vea_idx(i)) = M3D::Identity()*sqrt(dt);
    }
    for(unsigned int i=0;i<mtState::nPose_;i++){
      G.template block<3,3>(mtState::pop_idx(i),mtNoise::pop_idx(i)) = M3D::Identity()*sqrt(dt);
      G.template block<3,3>(mtState::poa_idx(i),mtNoise::poa_idx(i)) = M3D::Identity()*sqrt(dt);
    }
    LWF::NormalVectorElement nOut;
    for(unsigned int i=0;i<mtState::nMax_;i++){
      const int camID = state.CfP(i).camID_;
      if(camID >= 0 && camID < mtState::nCam_){
        oldC_ = state.CfP(i);
        oldD_ = state.dep(i);
        const V3D camRor = state.qCM(camID)*imuRor;
        const V3D camVel = state.qCM(camID)*V3D(imuRor.cross(state.MrMC(camID))-state.MvM());
        V3D dm = -dt*(gSM(oldC_.get_nor().getVec())*camVel/oldD_.getDistance()
            + (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())*camRor);
        QPD qm = exp_map(dm);
        nOut = oldC_.get_nor().rotated(qm);
        G(mtState::fea_idx(i)+2,mtNoise::fea_idx(i)+2) = sqrt(dt);
        G.template block<1,3>(mtState::fea_idx(i)+2,mtNoise::att_idx_) =
            sqrt(dt)*oldD_.getParameterDerivative()*oldC_.get_nor().getVec().transpose()*gSM(state.qCM(camID)*state.MrMC(camID))*MPD(state.qCM(camID)).matrix();
        G.template block<2,2>(mtState::fea_idx(i),mtNoise::fea_idx(i)) =
            -nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*oldC_.get_nor().getN()*sqrt(dt);
        G.template block<2,3>(mtState::fea_idx(i),mtNoise::att_idx_) =
            -nOut.getM().transpose()*gSM(qm*oldC_.get_nor().getVec())*Lmat(dm)*(
                - (M3D::Identity()-oldC_.get_nor().getVec()*oldC_.get_nor().getVec().transpose())
                +1.0/oldD_.getDistance()*gSM(oldC_.get_nor().getVec())*gSM(state.qCM(camID)*state.MrMC(camID))
             )*sqrt(dt)*MPD(state.qCM(camID)).matrix();
      }
    }
  }
  bool detectInertialMotion(const mtState& state, const mtMeas& meas) const{
    const V3D imuRor = meas.gyr()-state.gyb();
    const V3D imuAcc = meas.acc()-state.acb()+state.qWM().inverse()*g_;
    return (imuRor.norm() > inertialMotionRorTh_) | (imuAcc.norm() > inertialMotionAccTh_);
  }
};

}


#endif /* ROVIO_IMUPREDICTION_HPP_ */
