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

#pragma once

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/FilterState.hpp"
#include <map>
#include <unordered_set>
#include "CoordinateTransform/FeatureOutput.hpp"
#include "rovio/FeatureManager.hpp"
#include "rovio/MultiCamera.hpp"

#include "rovio_states.hpp"

namespace rovio {

struct StateAuxiliary {
    StateAuxiliary() {
        MwWMest_.setZero();
        MwWMmeas_.setZero();
        wMeasCov_.setIdentity();
        for (unsigned int i = 0; i < conf::nMax; i++) {
            A_red_[i].setIdentity();
            b_red_[i].setZero();
        }
        doVECalibration_ = true;
        activeFeature_ = 0;
        activeCameraCounter_ = 0;
        timeSinceLastInertialMotion_ = 0.0;
        timeSinceLastImageMotion_ = 0.0;
        for (unsigned int i = 0; i < conf::nCam; i++) {
            qCM_[i].setIdentity();
            MrMC_[i].setZero();
        }
        poseMeasRot_.setIdentity();
        poseMeasLin_.setZero();
    };

    V3D MwWMest_;  /**<Estimated rotational rate.*/
    V3D MwWMmeas_;  /**<Measured rotational rate.*/
    M3D wMeasCov_;  /**<Covariance of the measured rotational rate.*/
    Eigen::Matrix2d A_red_[conf::nMax];  /**<Reduced Jacobian of the pixel intensities w.r.t. to pixel coordinates, needed for the multilevel patch alignment. \see rovio::MultilevelPatchFeature::A_ \see rovio::getLinearAlignEquationsReduced()*/
    Eigen::Vector2d b_red_[conf::nMax];  /**<Reduced intensity errors, needed for the multilevel patch alignment. \see rovio::MultilevelPatchFeature::A_ \see rovio::getLinearAlignEquationsReduced()*/
    FeatureCoordinates feaCoorMeas_[conf::nMax];  /**<Intermediate variable for storing the measured feature location.*/
    QPD qCM_[conf::nCam];  /**<Quaternion Array: IMU coordinates to camera coordinates.*/
    V3D MrMC_[conf::nCam];  /**<Position Vector Array: Vectors pointing from IMU to the camera frame, expressed in the IMU frame.*/
    bool doVECalibration_;  /**<Do Camera-IMU extrinsic parameter calibration?*/
    int activeFeature_;  /**< Active Feature ID. ID of the currently updated feature. Needed in the image update procedure.*/
    int activeCameraCounter_;  /**<Counter for iterating through the cameras, used such that when updating a feature we always start with the camId where the feature is expressed in.*/
    double timeSinceLastInertialMotion_;  /**<Time since the IMU showed motion last.*/
    double timeSinceLastImageMotion_;  /**<Time since the Image showed motion last.*/
    QPD poseMeasRot_; /**<Groundtruth attitude measurement. qMI.*/
    Eigen::Vector3d poseMeasLin_; /**<Groundtruth position measurement. IrIM*/
    FeatureManager<rovio::conf::nLevels, rovio::conf::patchSize, rovio::conf::nCam> *mpCurrentFeature_; /**<Pointer to active feature*/
};

struct State : public StateManifold {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static constexpr int nMax_ = conf::nMax;            /**<Max number of features.*/
    static constexpr int nLevels_ = conf::nLevels;      /**<Max number of image levels.*/
    static constexpr int patchSize_ = conf::patchSize;  /**<Patch size.*/
    static constexpr int nCam_ = conf::nCam;            /**<Total number of cameras.*/
    static constexpr int nPose_ = conf::nPose;          /**<Total number of addtional pose states.*/

    virtual ~State() {};

    StateAuxiliary aux_;

    inline StateAuxiliary &aux() {
        return this->aux_;
    }

    inline const StateAuxiliary &aux() const {
        return this->aux_;
    }

    /** \brief Initializes the feature manager pointer as far as possible
     *
     *  @param featureManager* A pointer to a featureManagerArray
     */
    void initFeatureManagers(FeatureManager<conf::nLevels, conf::patchSize, conf::nCam> **mpFeatureManager) {
        for (int i = 0; i < conf::nMax; i++) {
            mpFeatureManager[i]->mpCoordinates_ = &CfP(i);
            mpFeatureManager[i]->mpDistance_ = &dep(i);
        }
    }

    /** \brief Initializes the feature manager pointer as far as possible
     *
     *  @param featureManager* A pointer to a featureManagerArray
     */
    void initFeatureManagers(FeatureSetManager<conf::nLevels, conf::patchSize, conf::nCam, conf::nMax> &fsm) {
        for (int i = 0; i < conf::nMax; i++) {
            fsm.features_[i].mpCoordinates_ = &CfP(i);
            fsm.features_[i].mpDistance_ = &dep(i);
        }
    }

    //@{
    /** \brief Get/Set the position vector pointing from the World-Frame to the IMU-Frame, expressed in World-Coordinates (World->IMU, expressed in World).
     *
     *  @return a reference to the position vector WrWM (World->IMU, expressed in World).
     */
    inline V3D &WrWM() {
        return this->pos();
    }

    inline const V3D &WrWM() const {
        return this->pos();
    }
    //@}

    //@{
    /** \brief Get/Set the absolute velocity vector of the IMU-Frame MvM, expressed in IMU-Coordinates.
     *
     *  @return a reference to the absolute velocity vector of the IMU-Frame MvM, expressed in IMU-Coordinates.
     */
    inline V3D &MvM() {
        return this->vel();
    }

    inline const V3D &MvM() const {
        return this->vel();
    }
    //@}

    //@{
    /** \brief Get/Set the quaternion qWM, expressing the  IMU-Frame in World-Coordinates (IMU Coordinates->World Coordinates).
     *
     *  @return a reference to the quaternion qWM (IMU Coordinates->World Coordinates).
     */
    inline QPD &qWM() {
        return this->att();
    }

    inline const QPD &qWM() const {
        return this->att();
    }
    //@}

    //@{
    /** \brief Get/Set the feature coordinates belonging to a specific feature i.
     *
     *  @param i - Feature Index
     *  @return a reference to the feature coordinates of feature i.
     */
    inline FeatureCoordinates &CfP(const int i = 0) {
        assert(i < nMax_);
        return this->fea(i).coordinates_;
    }

    inline const FeatureCoordinates &CfP(const int i = 0) const {
        assert(i < nMax_);
        return this->fea(i).coordinates_;
    }
    //@}

    //@{
    /** \brief Get/Set the quaternion qCM, expressing the IMU-Frame in Camera-Coordinates (IMU Coordinates->%Camera Coordinates).
     *
     *  @param camID - %Camera ID
     *  @return a reference to the quaternion qCM (IMU Coordinates->%Camera Coordinates).
     */
    inline QPD &qCM(const int camID = 0) {
        assert(camID < nCam_);
        if (this->aux().doVECalibration_) {
            return this->vea(camID);
        } else {
            return this->aux().qCM_[camID];
        }
    }

    inline const QPD &qCM(const int camID = 0) const {
        assert(camID < nCam_);
        if (this->aux().doVECalibration_) {
            return this->vea(camID);
        } else {
            return this->aux().qCM_[camID];
        }
    }
    //@}

    //@{
    /** \brief Get/Set the position vector pointing from the IMU-Frame to the Camera-Frame, expressed in IMU-Coordinates (IMU->%Camera, expressed in IMU).
     *
     *  @param camID - %Camera ID
     *  @return a reference to the position vector MrMC (IMU->%Camera, expressed in IMU).
     */
    inline V3D &MrMC(const int camID = 0) {
        assert(camID < nCam_);
        if (this->aux().doVECalibration_) {
            return this->vep(camID);
        } else {
            return this->aux().MrMC_[camID];
        }
    }

    inline const V3D &MrMC(const int camID = 0) const {
        assert(camID < nCam_);
        if (this->aux().doVECalibration_) {
            return this->vep(camID);
        } else {
            return this->aux().MrMC_[camID];
        }
    }
    //@}

    //@{
    /** \brief Get the position vector pointing from the World-Frame to the Camera-Frame, expressed in World-Coordinates (World->%Camera, expressed in World).
     *
     *  @param camID - %Camera ID
     *  @return the position vector WrWC (World->%Camera, expressed in World).
     */
    inline V3D WrWC(const int camID = 0) const {
        assert(camID < nCam_);
        return this->pos() + this->att() * MrMC(camID);
    }

    //@}

    //@{
    /** \brief Get the quaternion qCW, expressing the World-Frame in Camera-Coordinates (World Coordinates->%Camera Coordinates).
     *
     *  @param camID - %Camera ID
     *  @return he quaternion qCW (World Coordinates->%Camera Coordinates).
     */
    inline QPD qCW(const int camID = 0) const {
        assert(camID < nCam_);
        return qCM(camID) * this->att().inverse();
    }
    //@}

    //@{
    /** \brief Get/Set the distance parameter of a specific feature i.
     *
     *  \note The distance parameter can be either defined as regular distance, inverse distance, logarithmic distance or hyperbolic distance.
     *        The kind of the distance encoding depends on the defined DepthMap.
     *
     *  @param i - Feature Index
     *  @return a reference to distance parameter of the feature.
     */
    inline FeatureDistance &dep(const int i) {
        assert(i < nMax_);
        return this->fea(i).distance_;
    }

    inline const FeatureDistance &dep(const int i) const {
        assert(i < nMax_);
        return this->fea(i).distance_;
    }
    //@}

    //@{
    /** \brief Get/Set the rotational part of the pose with index i.
     *
     *  @param i - Pose index
     *  @return  a reference to the rotational part of the pose with index i.
     */
    inline QPD &poseRot(const int i) {
        assert(i < nPose_);
        return this->poa(i);
    }

    inline const QPD &poseRot(const int i) const {
        assert(i < nPose_);
        return this->poa(i);
    }
    //@}

    //@{
    /** \brief Get/Set the linear part of the pose with index i.
     *
     *  @param i - Pose index
     *  @return a reference to the linear part of the pose with index i.
     */
    inline V3D &poseLin(const int i = 0) {
        assert(i < nPose_);
        return this->pop(i);
    }

    inline const V3D &poseLin(const int i = 0) const {
        assert(i < nPose_);
        return this->pop(i);
    }
    //@}

    //@{
    /** \brief Get the position vector pointing from the World-Frame to the Camera-Frame, expressed in World-Coordinates (World->%Camera, expressed in World).
     *
     *  @note - This is compute based on the external pose measurement
     *  @return the position vector WrWC (World->%Camera, expressed in World).
     */
    inline V3D WrWC_ext() const {
        return this->aux().poseMeasLin_;
    }

    //@}

    //@{
    /** \brief Get the quaternion qCW, expressing the World-Frame in Camera-Coordinates (World Coordinates->%Camera Coordinates).
     *
     *  @note - This is compute based on the external pose measurement
     *  @return he quaternion qCW (World Coordinates->%Camera Coordinates).
     */
    inline QPD qCW_ext() const {
        return this->aux().poseMeasRot_;
    }
    //@}

    /** \brief Update the extrinsics of the MultiCamera
     *
     *  @param mpMultiCamera - Pointer to the multicamera object
     */
    void updateMultiCameraExtrinsics(MultiCamera<conf::nCam> *mpMultiCamera) const {
        for (int i = 0; i < conf::nCam; i++) {
            mpMultiCamera->BrBC_[i] = MrMC(i);
            mpMultiCamera->qCB_[i] = qCM(i);
        }
    }
};

}