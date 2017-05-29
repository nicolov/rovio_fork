/*
 * TestClasses.hpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Bloeschm
 */

#ifndef LWF_TestClasses_HPP_
#define LWF_TestClasses_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/State.hpp"
#include "lightweight_filtering/FilterState.hpp"
#include "lightweight_filtering/Update.hpp"
#include "lightweight_filtering/Prediction.hpp"
#include "lightweight_filtering/GIFPrediction.hpp"

#include "test_states.hpp"

namespace LWFTest{

namespace Nonlinear{

class OutlierDetectionExample: public LWF::OutlierDetection<LWF::ODEntry<0,3,1>>{
 public:
  virtual ~OutlierDetectionExample(){};
};
class FilterState: public LWF::FilterState<State,PredictionMeas,PredictionNoise,UpdateNoise::D_>{
 public:
  virtual ~FilterState(){};
};

class UpdateExample: public LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,false>{
 public:
  using mtState = LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,false>::mtState;
  typedef UpdateMeas mtMeas;
  typedef UpdateNoise mtNoise;
  typedef Innovation mtInnovation;
  UpdateExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~UpdateExample(){};
  void evalInnovation(mtInnovation& inn, const mtState& state, const mtNoise& noise) const{
    inn.pos_.get() = state.att_.get().rotate(state.pos_.get())-meas_.pos_.get()+noise.pos_.get();
    inn.att_.get() = (state.att_.get()*meas_.att_.get().inverted()).boxPlus(noise.att_.get());
  }
  void jacState(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = MPD(state.att_.get()).matrix();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::att_idx_) = -gSM(state.att_.get().rotate(state.pos_.get()));
    J.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = M3D::Identity();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtInnovation::att_idx_,mtNoise::att_idx_) = M3D::Identity();
  }
};

class PredictionExample: public LWF::Prediction<FilterState>{
 public:
  using mtState = LWF::Prediction<FilterState>::mtState;
  typedef PredictionMeas mtMeas;
  typedef PredictionNoise mtNoise;
  PredictionExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~PredictionExample(){};
  void evalPrediction(mtState& output, const mtState& state, const mtNoise& noise, double dt) const{
    V3D g_(0,0,-9.81);
    V3D dOmega = dt*(meas_.gyr_.get()-state.gyb_.get()-noise.att_.get()/sqrt(dt));
    QPD dQ = dQ.exponentialMap(dOmega);
    output.pos_.get() = state.pos_.get()+dt*state.att_.get().rotate(V3D(state.vel_.get()+noise.pos_.get()/sqrt(dt)));
    output.vel_.get() = (M3D::Identity()-gSM(dOmega))*state.vel_.get()
                       +dt*(meas_.acc_.get()-state.acb_.get()+state.att_.get().inverseRotate(g_)-noise.vel_.get()/sqrt(dt));
    output.acb_.get() = state.acb_.get()+noise.acb_.get()*sqrt(dt);
    output.gyb_.get() = state.gyb_.get()+noise.gyb_.get()*sqrt(dt);
    output.att_.get() = state.att_.get()*dQ;
    output.att_.get().fix();
  }
  void jacPreviousState(Eigen::MatrixXd& J, const mtState& state, double dt) const{
    V3D g_(0,0,-9.81);
    V3D dOmega = dt*(meas_.gyr_.get()-state.gyb_.get());
    J.setZero();
    J.template block<3,3>(mtState::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtState::pos_idx_,mtState::vel_idx_) = dt*MPD(state.att_.get()).matrix();
    J.template block<3,3>(mtState::pos_idx_,mtState::att_idx_) = -dt*gSM(V3D(state.att_.get().rotate(state.vel_.get())));
    J.template block<3,3>(mtState::vel_idx_,mtState::vel_idx_) = (M3D::Identity()-gSM(dOmega));
    J.template block<3,3>(mtState::vel_idx_,mtState::acb_idx_) = -dt*M3D::Identity();
    J.template block<3,3>(mtState::vel_idx_,mtState::gyb_idx_) = -dt*gSM(state.vel_.get());
    J.template block<3,3>(mtState::vel_idx_,mtState::att_idx_) = dt*MPD(state.att_.get()).matrix().transpose()*gSM(g_);
    J.template block<3,3>(mtState::acb_idx_,mtState::acb_idx_) = M3D::Identity();
    J.template block<3,3>(mtState::gyb_idx_,mtState::gyb_idx_) = M3D::Identity();
    J.template block<3,3>(mtState::att_idx_,mtState::gyb_idx_) = -dt*MPD(state.att_.get()).matrix()*Lmat(dOmega);
    J.template block<3,3>(mtState::att_idx_,mtState::att_idx_) = M3D::Identity();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state, double dt) const{
    mtNoise noise;
    V3D g_(0,0,-9.81);
    V3D dOmega = dt*(meas_.gyr_.get()-state.gyb_.get());
    J.setZero();
    J.template block<3,3>(mtState::pos_idx_,mtNoise::pos_idx_) = MPD(state.att_.get()).matrix()*sqrt(dt);
    J.template block<3,3>(mtState::vel_idx_,mtNoise::vel_idx_) = -M3D::Identity()*sqrt(dt);
    J.template block<3,3>(mtState::vel_idx_,mtNoise::att_idx_) = -gSM(state.vel_.get())*sqrt(dt);
    J.template block<3,3>(mtState::acb_idx_,mtNoise::acb_idx_) = M3D::Identity()*sqrt(dt);
    J.template block<3,3>(mtState::gyb_idx_,mtNoise::gyb_idx_) = M3D::Identity()*sqrt(dt);
    J.template block<3,3>(mtState::att_idx_,mtNoise::att_idx_) = -MPD(state.att_.get()).matrix()*Lmat(dOmega)*sqrt(dt);
  }
};

class PredictAndUpdateExample: public LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,true>{
 public:
  using mtState = LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,true>::mtState;
  typedef UpdateMeas mtMeas;
  typedef UpdateNoise mtNoise;
  typedef Innovation mtInnovation;
  PredictAndUpdateExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~PredictAndUpdateExample(){};
  void evalInnovation(mtInnovation& inn, const mtState& state, const mtNoise& noise) const{
    inn.pos_.get() = state.att_.get().rotate(state.pos_.get())-meas_.pos_.get()+noise.pos_.get();
    inn.att_.get() = (state.att_.get()*meas_.att_.get().inverted()).boxPlus(noise.att_.get());
  }
  void jacState(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = MPD(state.att_.get()).matrix();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::att_idx_) = -gSM(state.att_.get().rotate(state.pos_.get()));
    J.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = M3D::Identity();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtInnovation::att_idx_,mtNoise::att_idx_) = M3D::Identity();
  }
};

class GIFPredictionExample: public LWF::GIFPrediction<FilterState,GIFInnovation,PredictionMeas,PredictionNoise>{
 public:
  using mtState = FilterState::mtState;
  typedef PredictionMeas mtMeas;
  typedef PredictionNoise mtNoise;
  typedef GIFInnovation mtInnovation;
  V3D g_;
  GIFPredictionExample(): g_(0,0,-9.81){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~GIFPredictionExample(){};
  void evalResidual(mtInnovation& inn, const mtState& state0, const mtState& state1, const mtNoise& noise, double dt) const{
    V3D dOmega = dt*(meas_.gyr_.get()-state0.gyb_.get()-noise.att_.get()/sqrt(dt));
    inn.pos_.get() = (state1.pos_.get()- state0.pos_.get())/dt - state0.att_.get().rotate(V3D(state0.vel_.get() + noise.pos_.get()/sqrt(dt)));
    inn.vel_.get() = (state1.vel_.get()- (M3D::Identity()-gSM(dOmega))*state0.vel_.get())/dt
        - (meas_.acc_.get()-state0.acb_.get()+state0.att_.get().inverseRotate(g_) - noise.vel_.get()/sqrt(dt));
    inn.acb_.get() = (state1.acb_.get()- state0.acb_.get())/dt + noise.acb_.get()/sqrt(dt);
    inn.gyb_.get() = (state1.gyb_.get()- state0.gyb_.get())/dt + noise.gyb_.get()/sqrt(dt);
    inn.att_.get() = (state0.att_.get().inverted()*state1.att_.get()).logarithmicMap()/dt - dOmega/dt;
  }
  void jacPreviousState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    V3D dOmega = -dt*(meas_.gyr_.get()-previousState.gyb_.get());
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::vel_idx_) = -MPD(previousState.att_.get()).matrix();
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::att_idx_) = gSM(V3D(previousState.att_.get().rotate(V3D(previousState.vel_.get()))));
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = -(M3D::Identity()+gSM(dOmega))/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::acb_idx_) = M3D::Identity();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::gyb_idx_) = gSM(previousState.vel_.get());
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::att_idx_) = -MPD(previousState.att_.get()).matrix().transpose()*gSM(g_);
    F.template block<3,3>(mtInnovation::acb_idx_,mtState::acb_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::gyb_idx_,mtState::gyb_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::att_idx_,mtState::gyb_idx_) = M3D::Identity();
    F.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = -Lmat((previousState.att_.get().inverted()*currentState.att_.get()).logarithmicMap()).inverse()*MPD(previousState.att_.get().inverted()).matrix()/dt;
  }
  void jacCurrentState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::acb_idx_,mtState::acb_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::gyb_idx_,mtState::gyb_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = Lmat((previousState.att_.get().inverted()*currentState.att_.get()).logarithmicMap()).inverse()*MPD(previousState.att_.get().inverted()).matrix()/dt;
  }
  void jacNoise(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    V3D dOmega = -dt*(meas_.gyr_.get()-previousState.gyb_.get());
    F.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = -MPD(previousState.att_.get()).matrix()/sqrt(dt);
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::vel_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::att_idx_) = gSM(previousState.vel_.get())/sqrt(dt);
    F.template block<3,3>(mtInnovation::acb_idx_,mtNoise::acb_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::gyb_idx_,mtNoise::gyb_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::att_idx_,mtNoise::att_idx_) = M3D::Identity()/sqrt(dt);
  }
  void getLinearizationPoint(mtState& currentState, const mtFilterState& filterState, const mtMeas& meas, double dt){
    V3D dOmega = dt*(meas_.gyr_.get()-filterState.state_.gyb_.get());
    QPD dQ = dQ.exponentialMap(dOmega);
    currentState.pos_.get() = filterState.state_.pos_.get()+dt*filterState.state_.att_.get().rotate(filterState.state_.vel_.get());
    currentState.vel_.get() = (M3D::Identity()-gSM(dOmega))*filterState.state_.vel_.get()         +dt*(meas_.acc_.get()-filterState.state_.acb_.get()+filterState.state_.att_.get().inverseRotate(g_));
    currentState.acb_.get() = filterState.state_.acb_.get();
    currentState.gyb_.get() = filterState.state_.gyb_.get();
    currentState.att_.get() = filterState.state_.att_.get()*dQ;
    currentState.att_.get().fix();
  };
};

class GIFPredictionExampleWithUpdate: public LWF::GIFPrediction<FilterState,GIFInnovationWithUpdate,GIFMeasWithUpdate,GIFNoiseWithUpdate>{
 public:
  using mtState = FilterState::mtState;
  typedef GIFMeasWithUpdate mtMeas;
  typedef GIFNoiseWithUpdate mtNoise;
  typedef GIFInnovationWithUpdate mtInnovation;
  V3D g_;
  GIFPredictionExampleWithUpdate(): g_(0,0,-9.81){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~GIFPredictionExampleWithUpdate(){};
  void evalResidual(mtInnovation& inn, const mtState& state0, const mtState& state1, const mtNoise& noise, double dt) const{
    V3D dOmega = dt*(meas_.gyr_.get()-state0.gyb_.get()-noise.att_.get()/sqrt(dt));
    inn.pos_.get() = (state1.pos_.get()- state0.pos_.get())/dt - state0.att_.get().rotate(V3D(state0.vel_.get() + noise.pos_.get()/sqrt(dt)));
    inn.vel_.get() = (state1.vel_.get()- (M3D::Identity()-gSM(dOmega))*state0.vel_.get())/dt
        - (meas_.acc_.get()-state0.acb_.get()+state0.att_.get().inverseRotate(g_) - noise.vel_.get()/sqrt(dt));
    inn.acb_.get() = (state1.acb_.get()- state0.acb_.get())/dt + noise.acb_.get()/sqrt(dt);
    inn.gyb_.get() = (state1.gyb_.get()- state0.gyb_.get())/dt + noise.gyb_.get()/sqrt(dt);
    inn.att_.get() = (state0.att_.get().inverted()*state1.att_.get()).logarithmicMap()/dt - dOmega/dt;
    inn.posu_.get() = state1.att_.get().rotate(state1.pos_.get())-meas_.pos_.get()+noise.posu_.get();
    inn.attu_.get() = (state1.att_.get()*meas_.att_.get().inverted()).boxPlus(noise.attu_.get());
  }
  void jacPreviousState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    V3D dOmega = -dt*(meas_.gyr_.get()-previousState.gyb_.get());
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::vel_idx_) = -MPD(previousState.att_.get()).matrix();
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::att_idx_) = gSM(V3D(previousState.att_.get().rotate(V3D(previousState.vel_.get()))));
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = -(M3D::Identity()+gSM(dOmega))/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::acb_idx_) = M3D::Identity();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::gyb_idx_) = gSM(previousState.vel_.get());
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::att_idx_) = -MPD(previousState.att_.get()).matrix().transpose()*gSM(g_);
    F.template block<3,3>(mtInnovation::acb_idx_,mtState::acb_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::gyb_idx_,mtState::gyb_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::att_idx_,mtState::gyb_idx_) = M3D::Identity();
    F.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = -Lmat((previousState.att_.get().inverted()*currentState.att_.get()).logarithmicMap()).inverse()*MPD(previousState.att_.get().inverted()).matrix()/dt;
  }
  void jacCurrentState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::acb_idx_,mtState::acb_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::gyb_idx_,mtState::gyb_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::att_idx_,mtState::att_idx_) = Lmat((previousState.att_.get().inverted()*currentState.att_.get()).logarithmicMap()).inverse()*MPD(previousState.att_.get().inverted()).matrix()/dt;
    F.template block<3,3>(mtInnovation::posu_idx_,mtState::pos_idx_) = MPD(currentState.att_.get()).matrix();
    F.template block<3,3>(mtInnovation::posu_idx_,mtState::att_idx_) = -gSM(currentState.att_.get().rotate(currentState.pos_.get()));
    F.template block<3,3>(mtInnovation::attu_idx_,mtState::att_idx_) = M3D::Identity();
  }
  void jacNoise(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    V3D dOmega = -dt*(meas_.gyr_.get()-previousState.gyb_.get());
    F.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = -MPD(previousState.att_.get()).matrix()/sqrt(dt);
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::vel_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::att_idx_) = gSM(previousState.vel_.get())/sqrt(dt);
    F.template block<3,3>(mtInnovation::acb_idx_,mtNoise::acb_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::gyb_idx_,mtNoise::gyb_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::att_idx_,mtNoise::att_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::posu_idx_,mtNoise::posu_idx_) = M3D::Identity();
    F.template block<3,3>(mtInnovation::attu_idx_,mtNoise::attu_idx_) = M3D::Identity();
  }
  void getLinearizationPoint(mtState& currentState, const mtFilterState& filterState, const mtMeas& meas, double dt){
    V3D dOmega = dt*(meas_.gyr_.get()-filterState.state_.gyb_.get());
    QPD dQ = dQ.exponentialMap(dOmega);
    currentState.pos_.get() = filterState.state_.pos_.get()+dt*filterState.state_.att_.get().rotate(filterState.state_.vel_.get());
    currentState.vel_.get() = (M3D::Identity()-gSM(dOmega))*filterState.state_.vel_.get()        +dt*(meas_.acc_.get()-filterState.state_.acb_.get()+filterState.state_.att_.get().inverseRotate(g_));
    currentState.acb_.get() = filterState.state_.acb_.get();
    currentState.gyb_.get() = filterState.state_.gyb_.get();
    currentState.att_.get() = filterState.state_.att_.get()*dQ;
    currentState.att_.get().fix();
  }
};

}

namespace Linear{

class OutlierDetectionExample: public LWF::OutlierDetection<LWF::ODEntry<0,3,1>>{
 public:
  virtual ~OutlierDetectionExample(){};
};
class FilterState: public LWF::FilterState<State,PredictionMeas,PredictionNoise,UpdateNoise::D_>{
 public:
  virtual ~FilterState(){};
};

class UpdateExample: public LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,false>{
 public:
  using mtState = LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,false>::mtState;
  typedef UpdateMeas mtMeas;
  typedef UpdateNoise mtNoise;
  typedef Innovation mtInnovation;
  UpdateExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~UpdateExample(){};
  void evalInnovation(mtInnovation& inn, const mtState& state, const mtNoise& noise) const{
    inn.pos_.get() = state.pos_.get()-meas_.pos_.get()+noise.pos_.get();
    inn.hei_.get() = V3D(0,0,1).dot(state.pos_.get())-meas_.hei_.get()+noise.hei_.get();
  }
  void jacState(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    J.template block<1,3>(mtInnovation::hei_idx_,mtState::pos_idx_) = V3D(0,0,1).transpose();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    J(mtInnovation::hei_idx_,mtNoise::hei_idx_) = 1.0;
  }
};

class PredictionExample: public LWF::Prediction<FilterState>{
 public:
  using mtState = LWF::Prediction<FilterState>::mtState;
  typedef PredictionMeas mtMeas;
  typedef PredictionNoise mtNoise;
  PredictionExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~PredictionExample(){};
  void evalPrediction(mtState& output, const mtState& state, const mtNoise& noise, double dt) const{
    output.pos_.get() = state.pos_.get()+dt*state.vel_.get()+noise.vel_.get()*sqrt(dt);
    output.vel_.get() = state.vel_.get()+dt*meas_.acc_.get()+noise.acc_.get()*sqrt(dt);
  }
  void jacPreviousState(Eigen::MatrixXd& J, const mtState& state, double dt) const{
    J.setZero();
    J.template block<3,3>(mtState::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    J.template block<3,3>(mtState::pos_idx_,mtState::vel_idx_) = dt*M3D::Identity();
    J.template block<3,3>(mtState::vel_idx_,mtState::vel_idx_) = M3D::Identity();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state, double dt) const{
    J.setZero();
    J.template block<3,3>(mtState::pos_idx_,mtNoise::vel_idx_) = M3D::Identity()*sqrt(dt);
    J.template block<3,3>(mtState::vel_idx_,mtNoise::acc_idx_) = M3D::Identity()*sqrt(dt);
  }
};

class PredictAndUpdateExample: public LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,true>{
 public:
  using mtState = LWF::Update<Innovation,FilterState,UpdateMeas,UpdateNoise,OutlierDetectionExample,true>::mtState;
  typedef UpdateMeas mtMeas;
  typedef UpdateNoise mtNoise;
  typedef Innovation mtInnovation;
  PredictAndUpdateExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~PredictAndUpdateExample(){};
  void evalInnovation(mtInnovation& inn, const mtState& state, const mtNoise& noise) const{
    inn.pos_.get() = state.pos_.get()-meas_.pos_.get()+noise.pos_.get();
    inn.hei_.get() = V3D(0,0,1).dot(state.pos_.get())-meas_.hei_.get()+noise.hei_.get();
  }
  void jacState(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    J.template block<1,3>(mtInnovation::hei_idx_,mtState::pos_idx_) = V3D(0,0,1).transpose();
  }
  void jacNoise(Eigen::MatrixXd& J, const mtState& state) const{
    mtInnovation inn;
    J.setZero();
    J.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    J(mtInnovation::hei_idx_,mtNoise::hei_idx_) = 1.0;
  }
};

class GIFPredictionExample: public LWF::GIFPrediction<FilterState,GIFInnovation,PredictionMeas,PredictionNoise>{
 public:
  using mtState = FilterState::mtState;
  typedef PredictionMeas mtMeas;
  typedef PredictionNoise mtNoise;
  typedef GIFInnovation mtInnovation;
  GIFPredictionExample(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~GIFPredictionExample(){};
  void evalResidual(mtInnovation& inn, const mtState& state0, const mtState& state1, const mtNoise& noise, double dt) const{
    inn.vel_.get() = (state1.pos_.get()- state0.pos_.get())/dt - state0.vel_.get() + noise.vel_.get()/sqrt(dt);
    inn.acc_.get() = (state1.vel_.get()- state0.vel_.get())/dt - meas_.acc_.get() + noise.acc_.get()/sqrt(dt);
  }
  void jacPreviousState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::pos_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = -M3D::Identity();
    F.template block<3,3>(mtInnovation::acc_idx_,mtState::vel_idx_) = -M3D::Identity()/dt;
  }
  void jacCurrentState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::pos_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::acc_idx_,mtState::vel_idx_) = M3D::Identity()/dt;
  }
  void jacNoise(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::vel_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::acc_idx_,mtNoise::acc_idx_) = M3D::Identity()/sqrt(dt);
  }
};

class GIFPredictionExampleWithUpdate: public LWF::GIFPrediction<FilterState,GIFInnovationWithUpdate,GIFMeasWithUpdate,GIFNoiseWithUpdate>{
 public:
  using mtState = FilterState::mtState;
  typedef GIFMeasWithUpdate mtMeas;
  typedef GIFNoiseWithUpdate mtNoise;
  typedef GIFInnovationWithUpdate mtInnovation;
  GIFPredictionExampleWithUpdate(){
    disablePreAndPostProcessingWarning_ = true;
  };
  virtual ~GIFPredictionExampleWithUpdate(){};
  void evalResidual(mtInnovation& inn, const mtState& state0, const mtState& state1, const mtNoise& noise, double dt) const{
    inn.vel_.get() = (state1.pos_.get()- state0.pos_.get())/dt - state0.vel_.get() + noise.vel_.get()/sqrt(dt);
    inn.acc_.get() = (state1.vel_.get()- state0.vel_.get())/dt - meas_.acc_.get() + noise.acc_.get()/sqrt(dt);
    inn.pos_.get() = state1.pos_.get()-meas_.pos_.get()+noise.pos_.get();
    inn.hei_.get() = V3D(0,0,1).dot(state1.pos_.get())-meas_.hei_.get()+noise.hei_.get();
  }
  void jacPreviousState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::pos_idx_) = -M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = -M3D::Identity();
    F.template block<3,3>(mtInnovation::acc_idx_,mtState::vel_idx_) = -M3D::Identity()/dt;
  }
  void jacCurrentState(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::pos_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::acc_idx_,mtState::vel_idx_) = M3D::Identity()/dt;
    F.template block<3,3>(mtInnovation::pos_idx_,mtState::pos_idx_) = M3D::Identity();
    F.template block<1,3>(mtInnovation::hei_idx_,mtState::pos_idx_) = V3D(0,0,1).transpose();
  }
  void jacNoise(Eigen::MatrixXd& F, const mtState& previousState, const mtState& currentState, double dt) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtNoise::vel_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::acc_idx_,mtNoise::acc_idx_) = M3D::Identity()/sqrt(dt);
    F.template block<3,3>(mtInnovation::pos_idx_,mtNoise::pos_idx_) = M3D::Identity();
    F(mtInnovation::hei_idx_,mtNoise::hei_idx_) = 1.0;
  }
};

}

class NonlinearTest{
 public:
  static const int id_ = 0;
  typedef Nonlinear::State mtState;
  typedef Nonlinear::FilterState mtFilterState;
  typedef Nonlinear::UpdateMeas mtUpdateMeas;
  typedef Nonlinear::UpdateNoise mtUpdateNoise;
  typedef Nonlinear::Innovation mtInnovation;
  typedef Nonlinear::PredictionNoise mtPredictionNoise;
  typedef Nonlinear::PredictionMeas mtPredictionMeas;
  typedef Nonlinear::UpdateExample mtUpdateExample;
  typedef Nonlinear::PredictionExample mtPredictionExample;
  typedef Nonlinear::PredictAndUpdateExample mtPredictAndUpdateExample;
  typedef Nonlinear::OutlierDetectionExample mtOutlierDetectionExample;
  typedef Nonlinear::GIFPredictionExample mtGIFPredictionExample;
  typedef Nonlinear::GIFPredictionExampleWithUpdate mtGIFPredictionExampleWithUpdate;
  typedef Nonlinear::GIFInnovationWithUpdate mtGIFInnovationWithUpdate;
  typedef Nonlinear::GIFNoiseWithUpdate mtGIFNoiseWithUpdate;
  typedef Nonlinear::GIFMeasWithUpdate mtGIFMeasWithUpdate;
  void init(mtState& state,mtUpdateMeas& updateMeas,mtPredictionMeas& predictionMeas){
    state.pos_.get() = V3D(2.1,-0.2,-1.9);
    state.vel_.get() = V3D(0.3,10.9,2.3);
    state.acb_.get() = V3D(0.3,10.9,2.3);
    state.gyb_.get() = V3D(0.3,10.9,2.3);
    state.att_.get() = QPD(4.0/sqrt(30.0),3.0/sqrt(30.0),1.0/sqrt(30.0),2.0/sqrt(30.0));
    updateMeas.pos_.get() = V3D(-1.5,12,5.23);
    updateMeas.att_.get() = QPD(3.0/sqrt(15.0),-1.0/sqrt(15.0),1.0/sqrt(15.0),2.0/sqrt(15.0));
    predictionMeas.acc_.get() = V3D(-5,2,17.3);
    predictionMeas.gyr_.get() = V3D(15.7,0.45,-2.3);
  }
  static void mergePredictionAndUpdateMeas(mtGIFMeasWithUpdate& merged, mtPredictionMeas& predictionMeas, mtUpdateMeas& updateMeas){
    merged.acc_.get() = predictionMeas.acc_.get();
    merged.gyr_.get() = predictionMeas.gyr_.get();
    merged.pos_.get() = updateMeas.pos_.get();
    merged.att_.get() = updateMeas.att_.get();
  }
};

class LinearTest{
 public:
  static const int id_ = 1;
  typedef Linear::State mtState;
  typedef Linear::FilterState mtFilterState;
  typedef Linear::UpdateMeas mtUpdateMeas;
  typedef Linear::UpdateNoise mtUpdateNoise;
  typedef Linear::Innovation mtInnovation;
  typedef Linear::PredictionNoise mtPredictionNoise;
  typedef Linear::PredictionMeas mtPredictionMeas;
  typedef Linear::UpdateExample mtUpdateExample;
  typedef Linear::PredictionExample mtPredictionExample;
  typedef Linear::PredictAndUpdateExample mtPredictAndUpdateExample;
  typedef Linear::OutlierDetectionExample mtOutlierDetectionExample;
  typedef Linear::GIFPredictionExample mtGIFPredictionExample;
  typedef Linear::GIFPredictionExampleWithUpdate mtGIFPredictionExampleWithUpdate;
  typedef Linear::GIFInnovationWithUpdate mtGIFInnovationWithUpdate;
  typedef Linear::GIFNoiseWithUpdate mtGIFNoiseWithUpdate;
  typedef Linear::GIFMeasWithUpdate mtGIFMeasWithUpdate;
  void init(mtState& state,mtUpdateMeas& updateMeas,mtPredictionMeas& predictionMeas){
    state.pos_.get() = V3D(2.1,-0.2,-1.9);
    state.vel_.get() = V3D(0.3,10.9,2.3);
    updateMeas.pos_.get() = V3D(-1.5,12,5.23);
    updateMeas.hei_.get() = 0.5;
    predictionMeas.acc_.get() = V3D(-5,2,17.3);
  }
  static void mergePredictionAndUpdateMeas(mtGIFMeasWithUpdate& merged, mtPredictionMeas& predictionMeas, mtUpdateMeas& updateMeas){
    merged.acc_.get() = predictionMeas.acc_.get();
    merged.pos_.get() = updateMeas.pos_.get();
    merged.hei_.get() = updateMeas.hei_.get();
  }
};

}

#endif /* LWF_TestClasses_HPP_ */
