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

#ifndef ROVIO_ZEROVELOCITYUPDATE_HPP_
#define ROVIO_ZEROVELOCITYUPDATE_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/Update.hpp"
#include "lightweight_filtering/State.hpp"
#include "rovio/FilterStates.hpp"

#include "rovio_states.hpp"

namespace rovio {

/** \brief Outlier Detection.
 */
class ZeroVelocityOutlierDetection: public LWF::OutlierDetection<LWF::ODEntry<ZeroVelocityInnovation::vel_idx_,3>>{
 public:
  virtual ~ZeroVelocityOutlierDetection(){};
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** \brief Class, holding the zero velocity update
 */
class ZeroVelocityUpdate: public LWF::Update<ZeroVelocityInnovation,
                                             rovio::FilterState,
                                             ZeroVelocityUpdateMeas,
                                             ZeroVelocityUpdateNoise,
                                             ZeroVelocityOutlierDetection,
                                             false>{
 public:
  typedef LWF::Update<ZeroVelocityInnovation,rovio::FilterState,ZeroVelocityUpdateMeas,
      ZeroVelocityUpdateNoise,ZeroVelocityOutlierDetection,false> Base;
  using Base::doubleRegister_;
  using Base::intRegister_;
  typedef typename Base::mtState mtState;
  typedef typename Base::mtFilterState mtFilterState;
  typedef typename Base::mtInnovation mtInnovation;
  typedef typename Base::mtMeas mtMeas;
  typedef typename Base::mtNoise mtNoise;
  typedef typename Base::mtOutlierDetection mtOutlierDetection;

  /** \brief Constructor.
   *
   *   Loads and sets the needed parameters.
   */
  ZeroVelocityUpdate(){
    intRegister_.removeScalarByStr("maxNumIteration");
    doubleRegister_.removeScalarByStr("alpha");
    doubleRegister_.removeScalarByStr("beta");
    doubleRegister_.removeScalarByStr("kappa");
    doubleRegister_.removeScalarByStr("updateVecNormTermination");
  };

  /** \brief Destructor
   */
  virtual ~ZeroVelocityUpdate(){};

  /** \brief Compute the inovvation term
   *
   *  @param mtInnovation - Class, holding innovation data.
   *  @param state        - Filter %State.
   *  @param meas         - Not Used.
   *  @param noise        - Additive discrete Gaussian noise.
   *  @param dt           - Not used.
   */
  void evalInnovation(mtInnovation& y, const mtState& state, const mtNoise& noise) const{
    y.vel() = state.MvM()+noise.vel();
  }

  /** \brief Computes the Jacobian for the update step of the filter.
   *
   *  @param F     - Jacobian for the update step of the filter.
   *  @param state - Filter state.
   *  @param meas  - Not used.
   *  @param dt    - Not used.
   */
  void jacState(MXD& F, const mtState& state) const{
    F.setZero();
    F.template block<3,3>(mtInnovation::vel_idx_,mtState::vel_idx_) = Eigen::Matrix3d::Identity();
  }

  /** \brief Computes the Jacobian for the update step of the filter w.r.t. to the noise variables
   *
   *  @param G     - Jacobian for the update step of the filter.
   *  @param state - Filter state.
   *  @param meas  - Not used.
   *  @param dt    - Not used.
   */
  void jacNoise(MXD& G, const mtState& state) const{
    G.setZero();
    G.template block<3,3>(mtInnovation::vel_idx_,mtNoise::vel_idx_) = Eigen::Matrix3d::Identity();
  }
};

}


#endif /* ROVIO_ZEROVELOCITYUPDATE_HPP_ */
