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

#ifndef ROVIOOUTPUT_HPP_
#define ROVIOOUTPUT_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/CoordinateTransform.hpp"

#include "rovio_states.hpp"

namespace rovio {

class CameraOutputCT:public LWF::CoordinateTransform<State,StandardOutput>{
 public:
  typedef LWF::CoordinateTransform<State,StandardOutput> Base;
  typedef typename Base::mtInput mtInput;
  typedef typename Base::mtOutput mtOutput;
  int camID_;
  CameraOutputCT(){
    camID_ = 0;
  };
  virtual ~CameraOutputCT(){};
  void evalTransform(mtOutput& output, const mtInput& input) const{
    // WrWC = WrWM + qWM*(MrMC)
    // CwWC = qCM*MwWM
    // CvC  = qCM*(MvM + MwWM x MrMC)
    // qCW  = qCM*qWM^T
    output.WrWB() = input.WrWM()+input.qWM()*input.MrMC(camID_);
    output.BwWB() = input.qCM(camID_)*V3D(input.aux().MwWMmeas_-input.gyb());
    output.BvB()  = input.qCM(camID_)*V3D(-input.MvM() + gSM(V3D(input.aux().MwWMmeas_-input.gyb()))*input.MrMC(camID_));
    output.qBW()  = input.qCM(camID_)*input.qWM().inverse();
  }
  void jacTransform(MXD& J, const mtInput& input) const{
    J.setZero();
    J.template block<3,3>(mtOutput::WrWB_idx_,mtInput::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtOutput::WrWB_idx_,mtInput::att_idx_) =
        -gSM(input.qWM()*input.MrMC(camID_));
    J.template block<3,3>(mtOutput::qBW_idx_,mtInput::att_idx_) =
        -MPD(input.qCM(camID_)*input.qWM().inverse()).matrix();
    J.template block<3,3>(mtOutput::BvB_idx_,mtInput::vel_idx_) = -MPD(input.qCM(camID_)).matrix();
    J.template block<3,3>(mtOutput::BvB_idx_,mtInput::gyb_idx_) = MPD(input.qCM(camID_)).matrix()
        * gSM(input.MrMC(camID_));
    J.template block<3,3>(mtOutput::BwWB_idx_,mtInput::gyb_idx_) = -MPD(input.qCM(camID_)).matrix();
    if(input.aux().doVECalibration_){
      J.template block<3,3>(mtOutput::WrWB_idx_,mtInput::vep_idx(camID_)) =
          MPD(input.qWM()).matrix();
      J.template block<3,3>(mtOutput::BvB_idx_,mtInput::vep_idx(camID_)) =
          MPD(input.qCM(camID_)).matrix()*gSM(V3D(input.aux().MwWMmeas_-input.gyb()));
      J.template block<3,3>(mtOutput::BwWB_idx_,mtInput::vea_idx(camID_)) =
          -gSM(input.qCM(camID_)*V3D(input.aux().MwWMmeas_-input.gyb()));
      J.template block<3,3>(mtOutput::BvB_idx_,mtInput::vea_idx(camID_)) =
          -gSM(input.qCM(camID_)*V3D(-input.MvM() + gSM(V3D(input.aux().MwWMmeas_-input.gyb()))*input.MrMC()));
      J.template block<3,3>(mtOutput::qBW_idx_,mtInput::vea_idx(camID_)) =
          M3D::Identity();
    }
  }
  void postProcess(MXD& cov,const mtInput& input){
    cov.template block<3,3>(mtOutput::BwWB_idx_,mtOutput::BwWB_idx_) += input.aux().wMeasCov_;
  }
};

class ImuOutputCT:public LWF::CoordinateTransform<State,StandardOutput>{
 public:
  typedef LWF::CoordinateTransform<State,StandardOutput> Base;
  typedef typename Base::mtInput mtInput;
  typedef typename Base::mtOutput mtOutput;
  ImuOutputCT(){};
  virtual ~ImuOutputCT(){};
  void evalTransform(mtOutput& output, const mtInput& input) const{
    // WrWM = WrWM
    // MwWM = MwWM
    // MvM  = MvM
    // qMW  = qWM^T
    output.WrWB() = input.WrWM();
    output.BwWB() = input.aux().MwWMmeas_-input.gyb();
    output.BvB()  = -input.MvM(); // minus is required!
    output.qBW()  = input.qWM().inverse();
  }
  void jacTransform(MXD& J, const mtInput& input) const{
    J.setZero();
    J.template block<3,3>(mtOutput::WrWB_idx_,mtInput::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtOutput::qBW_idx_,mtInput::att_idx_) = -MPD(input.qWM().inverse()).matrix();
    J.template block<3,3>(mtOutput::BvB_idx_,mtInput::vel_idx_) = -M3D::Identity();
    J.template block<3,3>(mtOutput::BwWB_idx_,mtInput::gyb_idx_) = -M3D::Identity();
  }
  void postProcess(MXD& cov,const mtInput& input){
    cov.template block<3,3>(mtOutput::BwWB_idx_,mtOutput::BwWB_idx_) += input.aux().wMeasCov_;
  }
};

}


#endif /* ROVIOOUTPUT_HPP_ */
