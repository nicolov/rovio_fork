/*
 * Common.hpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Bloeschm
 */

#ifndef LWF_COMMON_HPP_
#define LWF_COMMON_HPP_

#include <map>
#include <type_traits>
#include <tuple>
#include <Eigen/Dense>
#include <iostream>
#include "lightweight_filtering/PropertyHandler.hpp"

typedef Eigen::Quaterniond QPD;
typedef Eigen::Matrix3d MPD;
typedef Eigen::Vector3d V3D;
typedef Eigen::Matrix3d M3D;
typedef Eigen::VectorXd VXD;
typedef Eigen::MatrixXd MXD;

inline M3D gSM(const V3D& vec){
  M3D mat;
  mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
  return mat;
}

static void enforceSymmetry(MXD& mat){
  mat = 0.5*(mat+mat.transpose()).eval();
}

inline M3D Lmat (const V3D& a) {
  const double norm = a.norm();
  const M3D skewMatrix = gSM(a);
  if (norm < 1.0e-4) {
    return M3D::Identity() + 0.5*skewMatrix;
  }
  return M3D::Identity() + (double(1.0) - cos(norm))/(norm*norm)*skewMatrix + (norm - sin(norm))/(norm*norm*norm)*(skewMatrix*skewMatrix);
}

inline V3D log_map(const QPD& q) {
    using std::acos;
    using std::sqrt;

    // define these compile time constants to avoid std::abs:
    static const double twoPi = 2.0 * M_PI, NearlyOne = 1.0 - 1e-10,
    NearlyNegativeOne = -1.0 + 1e-10;

    V3D omega;

    const double qw = q.w();
    // See Quaternion-Logmap.nb in doc for Taylor expansions
    if (qw > NearlyOne) {
      // Taylor expansion of (angle / s) at 1
      // (2 + 2 * (1-qw) / 3) * q.vec();
      omega = ( 8. / 3. - 2. / 3. * qw) * q.vec();
    } else if (qw < NearlyNegativeOne) {
      // Taylor expansion of (angle / s) at -1
      // (-2 - 2 * (1 + qw) / 3) * q.vec();
      omega = (-8. / 3. - 2. / 3. * qw) * q.vec();
    } else {
      // Normal, away from zero case
      double angle = 2 * acos(qw), s = sqrt(1 - qw * qw);
      // Important:  convert to [-pi,pi] to keep error continuous
      if (angle > M_PI)
      angle -= twoPi;
      else if (angle < -M_PI)
      angle += twoPi;
      omega = (angle / s) * q.vec();
    }

    return omega;
}

inline QPD exp_map(const V3D& omega) {
    using std::cos;
    using std::sin;

    double theta2 = omega.dot(omega);
    if (theta2 > std::numeric_limits<double>::epsilon()) {
      double theta = std::sqrt(theta2);
      double ha = 0.5 * theta;
      V3D vec = (sin(ha) / theta) * omega;
      return QPD(cos(ha), vec.x(), vec.y(), vec.z());
    } else {
      // first order approximation sin(theta/2)/theta = 0.5
      V3D vec = 0.5 * omega;
      return QPD(1.0, vec.x(), vec.y(), vec.z());
    }
}

inline QPD box_plus(const QPD& q, const V3D& v) {
  return exp_map(v) * q;
}

inline V3D box_minus(const QPD& q1, const QPD& q2) {
  return log_map(q1 * q2.inverse());
}

namespace LWF{
  enum FilteringMode{
    ModeEKF,
    ModeUKF,
    ModeIEKF
  };
}

#endif /* LWF_COMMON_HPP_ */
