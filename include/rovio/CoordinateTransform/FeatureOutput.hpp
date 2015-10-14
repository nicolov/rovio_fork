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

#ifndef ROVIO_FEATUREOUTPUT_HPP_
#define ROVIO_FEATUREOUTPUT_HPP_

#include "rovio/FilterStates.hpp"
#include "rovio/RobocentricFeatureElement.hpp"
#include "rovio/MultiCamera.hpp"
#include "lightweight_filtering/CoordinateTransform.hpp"

namespace rovio {

class FeatureOutput: public LWF::State<RobocentricFeatureElement>{
 public:
  static constexpr unsigned int _fea = 0;
  FeatureOutput(){
    static_assert(_fea+1==E_,"Error with indices");
    this->template getName<_fea>() = "fea";
  }
  ~FeatureOutput(){};

  //@{
  /** \brief Get/Set the distance parameter
   *
   *  @return a reference to distance parameter of the feature.
   */
  inline FeatureDistance& d(){
    return this->template get<_fea>().distance_;
  }
  inline const FeatureDistance& d() const{
    return this->template get<_fea>().distance_;
  }
  //@}

  //@{
  /** \brief Get/Set the feature coordinates
   *
   *  @return a reference to the feature coordinates
   */
  inline FeatureCoordinates& c() {
    return this->template get<_fea>().coordinates_;
  }
  inline const FeatureCoordinates& c() const{
    return this->template get<_fea>().coordinates_;
  }
  //@}
};

template<typename STATE>
class FeatureOutputCT:public LWF::CoordinateTransform<STATE,FeatureOutput,true>{
 public:
  typedef LWF::CoordinateTransform<STATE,FeatureOutput,true> Base;
  typedef typename Base::mtInput mtInput;
  typedef typename Base::mtOutput mtOutput;
  typedef typename Base::mtMeas mtMeas;
  typedef typename Base::mtJacInput mtJacInput;
  int ID_;
  FeatureOutputCT(){
    ID_ = -1;
  };
  ~FeatureOutputCT(){};
  void setFeatureID(int ID){
    ID_ = ID;
  }
  void eval(mtOutput& output, const mtInput& input, const mtMeas& meas, double dt = 0.0) const{
    output.template get<mtOutput::_fea>() = input.template get<mtInput::_fea>(ID_);
  }
  void jacInput(mtJacInput& J, const mtInput& input, const mtMeas& meas, double dt = 0.0) const{
    J.setIdentity();
  }
};

template<typename STATE>
class TransformFeatureOutputCT:public LWF::CoordinateTransform<STATE,FeatureOutput,true>{
 public:
  typedef LWF::CoordinateTransform<STATE,FeatureOutput,true> Base;
  typedef typename Base::mtInput mtInput;
  typedef typename Base::mtOutput mtOutput;
  typedef typename Base::mtMeas mtMeas;
  typedef typename Base::mtJacInput mtJacInput;
  int outputCamID_;
  int ID_;
  bool ignoreDistanceOutput_;
  MultiCamera<STATE::nCam_>* mpMultiCamera_;
  TransformFeatureOutputCT(MultiCamera<STATE::nCam_>* mpMultiCamera){
    mpMultiCamera_ = mpMultiCamera;
    outputCamID_ = 0;
    ID_ = -1;
    ignoreDistanceOutput_ = false;
  };
  ~TransformFeatureOutputCT(){};
  void setFeatureID(int ID){
    ID_ = ID;
  }
  void setOutputCameraID(int camID){
    outputCamID_ = camID;
  }
  void eval(mtOutput& output, const mtInput& input, const mtMeas& meas, double dt = 0.0) const{
    // D = Destination camera frame.
    // qDC = qDM*qCM^T
    // CrCD = qCM*(MrMD-MrMC)
    // DrDP = qDC*(d_in*nor_in-CrCD)
    // d_out = ||DrDP||
    // nor_out = DrDP/d_out
    input.updateMultiCameraExtrinsics(mpMultiCamera_);
    mpMultiCamera_->transformFeature(outputCamID_,input.CfP(ID_),input.dep(ID_),output.c(),output.d());
    output.c().camID_ = outputCamID_;
    output.c().mpCamera_ = &mpMultiCamera_->cameras_[outputCamID_];
  }
  void jacInput(mtJacInput& J, const mtInput& input, const mtMeas& meas, double dt = 0.0) const{
    input.updateMultiCameraExtrinsics(mpMultiCamera_);
    J.setZero();
    const int& camID = input.CfP(ID_).camID_;
    const QPD qDC = input.qCM(outputCamID_)*input.qCM(camID).inverted(); // TODO: avoid double computation
    const V3D CrCD = input.qCM(camID).rotate(V3D(input.MrMC(outputCamID_)-input.MrMC(camID)));
    const V3D CrCP = input.dep(ID_).getDistance()*input.CfP(ID_).get_nor().getVec();
    const V3D DrDP = qDC.rotate(V3D(CrCP-CrCD));
    const double d_out = DrDP.norm();
    const LWF::NormalVectorElement nor_out(DrDP);

    const Eigen::Matrix<double,1,3> J_d_DrDP = nor_out.getVec().transpose();
    const Eigen::Matrix<double,2,3> J_nor_DrDP = nor_out.getM().transpose()*(Eigen::Matrix3d::Identity()-nor_out.getVec()*nor_out.getVec().transpose())/d_out; // TODO: evtl delete projection
    const Eigen::Matrix<double,3,3> J_DrDP_qDC = gSM(DrDP);
    const Eigen::Matrix<double,3,3> J_DrDP_CrCP = MPD(qDC).matrix();
    const Eigen::Matrix<double,3,3> J_DrDP_CrCD = -MPD(qDC).matrix();

    const Eigen::Matrix<double,3,3> J_qDC_qDB = Eigen::Matrix3d::Identity();
    const Eigen::Matrix<double,3,3> J_qDC_qCB = -MPD(qDC).matrix();
    const Eigen::Matrix<double,3,3> J_CrCD_qCB = gSM(CrCD);
    const Eigen::Matrix<double,3,3> J_CrCD_BrBC = -MPD(input.qCM(camID)).matrix();
    const Eigen::Matrix<double,3,3> J_CrCD_BrBD = MPD(input.qCM(camID)).matrix();

    const Eigen::Matrix<double,3,1> J_CrCP_d = input.CfP(ID_).get_nor().getVec()*input.dep(ID_).getDistanceDerivative();
    const Eigen::Matrix<double,3,2> J_CrCP_nor = input.dep(ID_).getDistance()*input.CfP(ID_).get_nor().getM();

    J.template block<2,2>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_fea>(ID_)) = J_nor_DrDP*J_DrDP_CrCP*J_CrCP_nor;
    J.template block<2,1>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_fea>(ID_)+2) = J_nor_DrDP*J_DrDP_CrCP*J_CrCP_d;
    if(!ignoreDistanceOutput_){
      J.template block<1,2>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_fea>(ID_)) = J_d_DrDP*J_DrDP_CrCP*J_CrCP_nor;
      J.template block<1,1>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_fea>(ID_)+2) = J_d_DrDP*J_DrDP_CrCP*J_CrCP_d;
    }

    if(input.aux().doVECalibration_ && camID != outputCamID_){
      J.template block<2,3>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_vea>(camID)) = J_nor_DrDP*(J_DrDP_qDC*J_qDC_qCB+J_DrDP_CrCD*J_CrCD_qCB);
      J.template block<2,3>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_vea>(outputCamID_)) = J_nor_DrDP*J_DrDP_qDC*J_qDC_qDB;
      J.template block<2,3>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_vep>(camID)) = J_nor_DrDP*J_DrDP_CrCD*J_CrCD_BrBC;
      J.template block<2,3>(mtOutput::template getId<mtOutput::_fea>(),mtInput::template getId<mtInput::_vep>(outputCamID_)) = J_nor_DrDP*J_DrDP_CrCD*J_CrCD_BrBD;
      if(!ignoreDistanceOutput_){
        J.template block<1,3>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_vea>(camID)) = J_d_DrDP*(J_DrDP_qDC*J_qDC_qCB+J_DrDP_CrCD*J_CrCD_qCB);
        J.template block<1,3>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_vea>(outputCamID_)) = J_d_DrDP*J_DrDP_qDC*J_qDC_qDB;
        J.template block<1,3>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_vep>(camID)) = J_d_DrDP*J_DrDP_CrCD*J_CrCD_BrBC;
        J.template block<1,3>(mtOutput::template getId<mtOutput::_fea>()+2,mtInput::template getId<mtInput::_vep>(outputCamID_)) = J_d_DrDP*J_DrDP_CrCD*J_CrCD_BrBD;
      }
    }
  }
};

}


#endif /* ROVIO_FEATUREOUTPUT_HPP_ */
