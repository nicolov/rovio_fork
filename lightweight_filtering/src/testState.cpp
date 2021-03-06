#include "lightweight_filtering/State.hpp"
#include "lightweight_filtering/common.hpp"
#include "gtest/gtest.h"
#include <assert.h>

#include "test_states.hpp"

// The fixture for testing class ScalarState
class ScalarElementTest : public virtual ::testing::Test {
 protected:
  ScalarElementTest():covMat_(1,1) {
    unsigned int s = 1;
    testElement1_.setRandom(s);
    testElement2_.setRandom(s);
  }
  virtual ~ScalarElementTest() {
  }
  LWF::ScalarElement testElement1_;
  LWF::ScalarElement testElement2_;
  LWF::ScalarElement::mtDifVec difVec_;
  MXD covMat_;
};

// Test constructors
TEST_F(ScalarElementTest, constructor) {
  LWF::ScalarElement testElement1;
}

// Test setIdentity and Identity
TEST_F(ScalarElementTest, setIdentity) {
  testElement1_.setIdentity();
  ASSERT_EQ(testElement1_.s_,0.0);
  ASSERT_EQ(LWF::ScalarElement::Identity().s_,0.0);
}

// Test plus and minus
TEST_F(ScalarElementTest, plusAndMinus) {
  testElement2_.boxMinus(testElement1_,difVec_);
  ASSERT_EQ(difVec_(0),testElement2_.s_-testElement1_.s_);
  LWF::ScalarElement testElement3;
  testElement1_.boxPlus(difVec_,testElement3);
  ASSERT_NEAR(testElement2_.s_,testElement3.s_,1e-6);
}

// Test minus Jacobian
TEST_F(ScalarElementTest, minusJac) {
  testElement2_.boxMinusJac(testElement1_,covMat_);
  ASSERT_NEAR((covMat_-Eigen::Matrix<double,1,1>::Identity()).norm(),0.0,1e-6);
}

// Test getValue
TEST_F(ScalarElementTest, accessors) {
  ASSERT_TRUE(testElement1_.get() == testElement1_.s_);
}

// Test operator=
TEST_F(ScalarElementTest, operatorEQ) {
  testElement2_ = testElement1_;
  ASSERT_NEAR(testElement2_.s_,testElement1_.s_,1e-6);
}

// The fixture for testing class VectorState
class VectorElementTest : public virtual ::testing::Test {
 protected:
  VectorElementTest():covMat_((int)N_,(int)N_) {
    unsigned int s = 1;
    testElement1_.setRandom(s);
    testElement2_.setRandom(s);
  }
  virtual ~VectorElementTest() {
  }
  static const unsigned int N_ = 4;
  LWF::VectorElement<N_> testElement1_;
  LWF::VectorElement<N_> testElement2_;
  LWF::VectorElement<N_>::mtDifVec difVec_;
  MXD covMat_;
};

// Test constructors
TEST_F(VectorElementTest, constructor) {
  LWF::VectorElement<N_> testElement1;
}

// Test setIdentity and Identity
TEST_F(VectorElementTest, setIdentity) {
  testElement1_.setIdentity();
  ASSERT_EQ(testElement1_.v_.norm(),0.0);
  ASSERT_EQ(LWF::VectorElement<N_>::Identity().v_.norm(),0.0);
}

// Test plus and minus
TEST_F(VectorElementTest, plusAndMinus) {
  testElement2_.boxMinus(testElement1_,difVec_);
  ASSERT_NEAR((difVec_-(testElement2_.v_-testElement1_.v_)).norm(),0.0,1e-6);
  LWF::VectorElement<N_> testElement3;
  testElement1_.boxPlus(difVec_,testElement3);
  ASSERT_NEAR((testElement2_.v_-testElement3.v_).norm(),0.0,1e-6);
}

// Test minus Jacobian
TEST_F(VectorElementTest, minusJac) {
  testElement2_.boxMinusJac(testElement1_,covMat_);
  ASSERT_NEAR((covMat_-Eigen::Matrix<double,N_,N_>::Identity()).norm(),0.0,1e-6);
}

// Test getValue
TEST_F(VectorElementTest, accessors) {
  ASSERT_TRUE(testElement1_.get() == testElement1_.v_);
}

// Test operator=
TEST_F(VectorElementTest, operatorEQ) {
  testElement2_ = testElement1_;
  ASSERT_NEAR((testElement2_.v_-testElement1_.v_).norm(),0.0,1e-6);
}

// The fixture for testing class QuaternionElementTest
class QuaternionElementTest : public virtual ::testing::Test {
 protected:
  QuaternionElementTest():covMat_(3,3) {
    unsigned int s = 1;
    testElement1_.setRandom(s);
    testElement2_.setRandom(s);
  }
  virtual ~QuaternionElementTest() {
  }
  LWF::QuaternionElement testElement1_;
  LWF::QuaternionElement testElement2_;
  LWF::QuaternionElement testElement3_;
  LWF::QuaternionElement::mtDifVec difVec_;
  LWF::QuaternionElement::mtDifVec difVecIn_;
  LWF::QuaternionElement::mtDifVec difVecOut1_;
  LWF::QuaternionElement::mtDifVec difVecOut2_;
  MXD covMat_;
};

// Test constructors
TEST_F(QuaternionElementTest, constructor) {
  LWF::QuaternionElement testElement1;
}

// Test setIdentity and Identity
TEST_F(QuaternionElementTest, setIdentity) {
  testElement1_.setIdentity();
  QPD idQ = QPD::Identity();

  ASSERT_TRUE(testElement1_.q_.isApprox(idQ, 1e-6));
  ASSERT_TRUE(LWF::QuaternionElement::Identity().q_.isApprox(idQ,1e-6));
}

// Test plus and minus
TEST_F(QuaternionElementTest, plusAndMinus) {
  testElement2_.boxMinus(testElement1_,difVec_);
  // Check that the difference vector matches the one computed on quaternions directly
  V3D quatDiffVec = log_map(testElement2_.q_ * testElement1_.q_.inverse());
  ASSERT_NEAR((difVec_-quatDiffVec).norm(),0.0,1e-6);
  testElement1_.boxPlus(difVec_,testElement3_);
  ASSERT_TRUE(testElement2_.q_.angularDistance(testElement3_.q_) < 1e-6);
}

// Test minus Jacobian
TEST_F(QuaternionElementTest, minusJac) {
  testElement2_.boxMinus(testElement1_,difVecOut1_);
  double d = 0.00001;
  testElement2_.boxMinusJac(testElement1_,covMat_);
  for(unsigned int i=0;i<3;i++){
    difVecIn_.setZero();
    difVecIn_(i) = d;
    testElement2_.boxPlus(difVecIn_,testElement3_);
    testElement3_.boxMinus(testElement1_,difVecOut2_);
    difVecOut2_ -= difVecOut1_;
    difVecOut2_ = difVecOut2_/d;
    ASSERT_NEAR((covMat_.col(i)-difVecOut2_).norm(),0.0,1e-5);
  }
}

// Test getValue
TEST_F(QuaternionElementTest, accessors) {
  ASSERT_TRUE(testElement1_.get().isApprox(testElement1_.q_,1e-6));
}

// Test operator=
TEST_F(QuaternionElementTest, operatorEQ) {
  testElement2_ = testElement1_;
  ASSERT_TRUE(testElement2_.q_.isApprox(testElement1_.q_,1e-6));
}

// Test LMat and derivative of rotation
TEST_F(QuaternionElementTest, LMat) {
  double d = 0.00001;
  LWF::QuaternionElement att;
  LWF::VectorElement<3> vec;
  vec.v_ = V3D(0.4,-0.2,1.7);
  M3D J;
  LWF::QuaternionElement attDisturbed;
  LWF::VectorElement<3> vecDisturbed;
  M3D I;
  I.setIdentity();
  V3D dif;
  I = d*I;
  att.q_ = exp_map(vec.v_);
  for(unsigned int i=0;i<3;i++){
    vec.boxPlus(I.col(i),vecDisturbed);
    attDisturbed.q_ = exp_map(vecDisturbed.v_);
    attDisturbed.boxMinus(att,dif);
    J.col(i) = dif*1/d;
  }
  ASSERT_NEAR((J-Lmat(vec.v_)).norm(),0.0,1e-5);
}

// Test LMat and derivative of rotation
TEST_F(QuaternionElementTest, DerivativeOfRotation) {
  double d = 1e-6;
  LWF::QuaternionElement att;
  unsigned int s = 1;
  att.setRandom(s);
  LWF::QuaternionElement attDisturbed;
  M3D J;
  V3D dif;
  V3D a(2.1,0.7,-1.7);
  V3D b = att.q_ * a;
  V3D b_disturbed;
  for(unsigned int i=0;i<3;i++){
    dif.setZero(); dif(i) = d;
    att.boxPlus(dif,attDisturbed);
    b_disturbed = attDisturbed.q_ * a;
    J.col(i) = (b_disturbed-b)/d;
  }
  ASSERT_NEAR((J+gSM(b)).norm(),0.0,1e-5);
}

// The fixture for testing class NormalVectorElementTest
class NormalVectorElementTest : public virtual ::testing::Test {
 protected:
  NormalVectorElementTest():covMat_(2,2) {
    unsigned int s = 1;
    testElement1_.setRandom(s);
    testElement2_.setRandom(s);
  }
  virtual ~NormalVectorElementTest() {
  }
  LWF::NormalVectorElement testElement1_;
  LWF::NormalVectorElement testElement2_;
  LWF::NormalVectorElement testElement3_;
  LWF::NormalVectorElement::mtDifVec difVec_;
  LWF::NormalVectorElement::mtDifVec difVecIn_;
  LWF::NormalVectorElement::mtDifVec difVecOut1_;
  LWF::NormalVectorElement::mtDifVec difVecOut2_;
  MXD covMat_;
};

// Test constructors
TEST_F(NormalVectorElementTest, constructor) {
  LWF::NormalVectorElement testElement1;
}

// Test setIdentity and Identity
TEST_F(NormalVectorElementTest, setIdentity) {
  testElement1_.setIdentity();
  ASSERT_TRUE(testElement1_.getVec() == V3D(0,0,1));
  ASSERT_TRUE(LWF::NormalVectorElement::Identity().getVec() == V3D(0,0,1));
}

// Test plus and minus
TEST_F(NormalVectorElementTest, plusAndMinus) {
  testElement2_.boxMinus(testElement1_,difVec_);
  LWF::NormalVectorElement testElement3;
  testElement1_.boxPlus(difVec_,testElement3);
  ASSERT_NEAR((testElement2_.getVec()-testElement3.getVec()).norm(),0.0,1e-6);
}

// Test minus Jacobian
TEST_F(NormalVectorElementTest, minusJac) {
  testElement2_.boxMinus(testElement1_,difVecOut1_);
  double d = 0.00001;
  testElement2_.boxMinusJac(testElement1_,covMat_);
  for(unsigned int i=0;i<2;i++){
    difVecIn_.setZero();
    difVecIn_(i) = d;
    testElement2_.boxPlus(difVecIn_,testElement3_);
    testElement3_.boxMinus(testElement1_,difVecOut2_);
    difVecOut2_ -= difVecOut1_;
    difVecOut2_ = difVecOut2_/d;
    ASSERT_NEAR((covMat_.col(i)-difVecOut2_).norm(),0.0,1e-4);
  }
}

// Test getValue
TEST_F(NormalVectorElementTest, accessors) {
  ASSERT_TRUE(testElement1_.get().getVec() == testElement1_.getVec());
}

// Test operator=
TEST_F(NormalVectorElementTest, operatorEQ) {
  testElement2_ = testElement1_;
  ASSERT_TRUE(testElement2_.getVec() == testElement1_.getVec());
}

// Test getPerpendiculars
TEST_F(NormalVectorElementTest, getTwoNormals) {
  ASSERT_NEAR(testElement1_.getPerp1().dot(testElement1_.getVec()),0.0,1e-6);
  ASSERT_NEAR(testElement1_.getPerp2().dot(testElement1_.getVec()),0.0,1e-6);
  ASSERT_NEAR(testElement1_.getPerp2().dot(testElement1_.getPerp1()),0.0,1e-6);
}

// Test set from vector
TEST_F(NormalVectorElementTest, setFromVector) {
  testElement2_.setFromVector(testElement1_.getVec());
  ASSERT_NEAR((testElement2_.getVec()-testElement1_.getVec()).norm(),0.0,1e-6);
}

// Test rotated
TEST_F(NormalVectorElementTest, rotated) {
  LWF::QuaternionElement q;
  unsigned int s = 1;
  q.setRandom(s);
  testElement2_.setFromVector(testElement1_.getVec());
  ASSERT_NEAR((q.q_ * testElement1_.getVec()-testElement1_.rotated(q.q_).getVec()).norm(),0.0,1e-6);
}

// Test inverted
TEST_F(NormalVectorElementTest, inverted) {
  ASSERT_NEAR((testElement1_.getVec()+testElement1_.inverted().getVec()).norm(),0.0,1e-6);
}

// Test derivative of vector
TEST_F(NormalVectorElementTest, derivative) {
  const double d  = 1e-6;
  difVec_.setZero();
  difVec_(0) = d;
  testElement1_.boxPlus(difVec_,testElement2_);
  ASSERT_NEAR(((testElement2_.getVec()-testElement1_.getVec())/d-testElement1_.getM().col(0)).norm(),0.0,1e-6);
  difVec_.setZero();
  difVec_(1) = d;
  testElement1_.boxPlus(difVec_,testElement2_);
  ASSERT_NEAR(((testElement2_.getVec()-testElement1_.getVec())/d-testElement1_.getM().col(1)).norm(),0.0,1e-6);
}

// Test getRotationFromTwoNormals
TEST_F(NormalVectorElementTest, getRotationFromTwoNormals) {
  V3D theta = LWF::NormalVectorElement::getRotationFromTwoNormals(testElement1_,testElement2_);
  QPD q = exp_map(theta);
  ASSERT_NEAR((testElement1_.rotated(q).getVec()-testElement2_.getVec()).norm(),0.0,1e-6);
  ASSERT_NEAR((LWF::NormalVectorElement::getRotationFromTwoNormals(testElement1_,testElement2_)+LWF::NormalVectorElement::getRotationFromTwoNormals(testElement2_,testElement1_)).norm(),0.0,1e-6);
}

// Test getRotationFromTwoNormalsJac
TEST_F(NormalVectorElementTest, getRotationFromTwoNormalsJac) {
  V3D testVec1 = testElement1_.getVec();
  V3D testVec2 = testElement2_.getVec();
  V3D perp = testElement1_.getPerp1();
  V3D theta = LWF::NormalVectorElement::getRotationFromTwoNormals(testVec1,testVec2,perp);
  covMat_ = LWF::NormalVectorElement::getRotationFromTwoNormalsJac(testVec1,testVec2);

  double d = 0.00001;
  for(unsigned int i=0;i<3;i++){
    V3D testVec3 = testVec1;
    testVec3(i) += d;
    V3D theta2 = LWF::NormalVectorElement::getRotationFromTwoNormals(testVec3,testVec2,perp);
    V3D diff = (theta2-theta)/d;
    ASSERT_NEAR((covMat_.col(i)-diff).norm(),0.0,1e-4);
  }
}

// Test derivative of boxminus (includes test of derivative of getRotationFromTwoNormals)
TEST_F(NormalVectorElementTest, derivativeBoxMinus) {
  const double d  = 1e-6;
  Eigen::Matrix2d J;
  LWF::NormalVectorElement::mtDifVec perturbation;
  LWF::NormalVectorElement::mtDifVec out;
  LWF::NormalVectorElement::mtDifVec out_perturbed;
  LWF::NormalVectorElement testElement2_perturbed;
  testElement2_.boxMinus(testElement1_,out);
  perturbation.setZero();
  perturbation(0) = d;
  testElement2_.boxPlus(perturbation,testElement2_perturbed);
  testElement2_perturbed.boxMinus(testElement1_,out_perturbed);
  J.col(0) = (out_perturbed-out)/d;
  perturbation.setZero();
  perturbation(1) = d;
  testElement2_.boxPlus(perturbation,testElement2_perturbed);
  testElement2_perturbed.boxMinus(testElement1_,out_perturbed);
  J.col(1) = (out_perturbed-out)/d;
  ASSERT_NEAR((testElement1_.getN().transpose()*-LWF::NormalVectorElement::getRotationFromTwoNormalsJac(testElement2_,testElement1_)*testElement2_.getM()-J).norm(),0.0,1e-5);
}

// The fixture for testing class ArrayElementTest
class ArrayElementTest : public virtual ::testing::Test {
 protected:
  ArrayElementTest():covMat_((int)N_*3,(int)N_*3) {
    unsigned int s = 1;
    testElement1_.setRandom(s);
    testElement2_.setRandom(s);
  }
  virtual ~ArrayElementTest() {
  }
  static const unsigned int N_ = 5;
  LWF::ArrayElement<LWF::QuaternionElement,N_> testElement1_;
  LWF::ArrayElement<LWF::QuaternionElement,N_> testElement2_;
  LWF::ArrayElement<LWF::QuaternionElement,N_> testElement3_;
  LWF::ArrayElement<LWF::QuaternionElement,N_>::mtDifVec difVec_;
  LWF::ArrayElement<LWF::QuaternionElement,N_>::mtDifVec difVecIn_;
  LWF::ArrayElement<LWF::QuaternionElement,N_>::mtDifVec difVecOut1_;
  LWF::ArrayElement<LWF::QuaternionElement,N_>::mtDifVec difVecOut2_;
  MXD covMat_;
};

// Test constructors
TEST_F(ArrayElementTest, constructor) {
  LWF::ArrayElement<LWF::QuaternionElement,N_> testElement1;
}

// Test setIdentity and Identity
TEST_F(ArrayElementTest, setIdentity) {
  testElement1_.setIdentity();
    QPD idQ = QPD::Identity();
  for(unsigned int i=0;i<N_;i++){
    ASSERT_TRUE(testElement1_.array_[i].q_.isApprox(idQ,1e-6));
    ASSERT_TRUE((LWF::ArrayElement<LWF::QuaternionElement,N_>::Identity().array_[i].q_.isApprox(idQ,1e-6)));
  }
}

// Test plus and minus
TEST_F(ArrayElementTest, plusAndMinus) {
  testElement2_.boxMinus(testElement1_,difVec_);
  for(unsigned int i=0;i<N_;i++){
    ASSERT_NEAR((difVec_.block<3,1>(i*3,0)-box_minus(testElement2_.array_[i].q_,testElement1_.array_[i].q_)).norm(),0.0,1e-6);
  }
  LWF::ArrayElement<LWF::QuaternionElement,N_> testElement3;
  testElement1_.boxPlus(difVec_,testElement3);
  for(unsigned int i=0;i<N_;i++){
    ASSERT_TRUE(testElement2_.array_[i].q_.angularDistance(testElement3.array_[i].q_) < 1e-6);
  }
}

// Test minus Jacobian
TEST_F(ArrayElementTest, minusJac) {
  testElement2_.boxMinus(testElement1_,difVecOut1_);
  double d = 0.00001;
  testElement2_.boxMinusJac(testElement1_,covMat_);
  for(unsigned int i=0;i<3*N_;i++){
    difVecIn_.setZero();
    difVecIn_(i) = d;
    testElement2_.boxPlus(difVecIn_,testElement3_);
    testElement3_.boxMinus(testElement1_,difVecOut2_);
    difVecOut2_ -= difVecOut1_;
    difVecOut2_ = difVecOut2_/d;
    ASSERT_NEAR((covMat_.col(i)-difVecOut2_).norm(),0.0,1e-5);
  }
}

// Test getValue
TEST_F(ArrayElementTest, accessors) {
  for(unsigned int i=0;i<N_;i++){
    ASSERT_TRUE(testElement1_.get(i).isApprox(testElement1_.array_[i].q_,1e-6));
  }
}

// Test operator=
TEST_F(ArrayElementTest, operatorEQ) {
  testElement2_ = testElement1_;
  for(unsigned int i=0;i<N_;i++){
    ASSERT_TRUE(testElement2_.array_[i].q_.isApprox(testElement1_.array_[i].q_,1e-6));
  }
}

class StateTesting : public virtual ::testing::Test {
 protected:
  StateTesting():covMat_((int)mtState::D_,(int)mtState::D_) {
    testScalar1_ = 4.5;
    testScalar2_ = -17.34;
    testVector1_[0] << 2.1, -0.2, -1.9;
    testVector2_[0] << -10.6, 0.2, -105.2;
    for(int i=1;i<4;i++){
      testVector1_[i] = testVector1_[i-1] + V3D(0.3,10.9,2.3);
      testVector2_[i] = testVector2_[i-1] + V3D(-1.5,12,1785.23);
    }
    testQuat1_[0] = QPD(4.0/sqrt(30.0),3.0/sqrt(30.0),1.0/sqrt(30.0),2.0/sqrt(30.0));
    testQuat2_[0] = QPD(0.0,0.36,0.48,0.8);
    for(int i=1;i<4;i++){
      testQuat1_[i] = box_plus(testQuat1_[i-1],testVector1_[i-1]);
      testQuat2_[i] = box_plus(testQuat2_[i-1],testVector2_[i-1]);
    }
    testState1_.sca_.get() = testScalar1_;
    testState1_.vec0_.get() = testVector1_[0];
    testState1_.vec1_.get() = testVector1_[1];
    testState1_.vec2_.get() = testVector1_[2];
    testState1_.vec3_.get() = testVector1_[3];
    testState1_.qua0_.get(0) = testQuat1_[0];
    testState1_.qua0_.get(1) = testQuat1_[1];
    testState1_.qua1_.get(0) = testQuat1_[2];
    testState1_.qua1_.get(1) = testQuat1_[3];
    testState2_.sca_.get() = testScalar2_;
    testState2_.vec0_.get() = testVector2_[0];
    testState2_.vec1_.get() = testVector2_[1];
    testState2_.vec2_.get() = testVector2_[2];
    testState2_.vec3_.get() = testVector2_[3];
    testState2_.qua0_.get(0) = testQuat2_[0];
    testState2_.qua0_.get(1) = testQuat2_[1];
    testState2_.qua1_.get(0) = testQuat2_[2];
    testState2_.qua1_.get(1) = testQuat2_[3];
  }
  virtual ~StateTesting() {
  }
  using mtState = LWFTest::TestState;
  mtState testState1_;
  mtState testState2_;
  mtState testState3_;
  mtState::mtDifVec difVec_;
  mtState::mtDifVec difVecIn_;
  mtState::mtDifVec difVecOut1_;
  mtState::mtDifVec difVecOut2_;
  MXD covMat_;
  double testScalar1_;
  double testScalar2_;
  V3D testVector1_[4];
  V3D testVector2_[4];
  QPD testQuat1_[4];
  QPD testQuat2_[4];
};

// Test constructors
TEST_F(StateTesting, constructors) {
  mtState testState1;
  mtState testState2(testState2_);
  testState2.boxMinus(testState2_,difVec_);
  ASSERT_NEAR(difVec_.norm(),0.0,1e-6);
}

// Test setIdentity and Identity
TEST_F(StateTesting, setIdentity) {
  testState1_.setIdentity();
  ASSERT_EQ(testState1_.sca_.get(),0);
  ASSERT_EQ(testState1_.vec0_.get().norm(),0);
  ASSERT_EQ(testState1_.vec1_.get().norm(),0);
  ASSERT_EQ(testState1_.vec2_.get().norm(),0);
  ASSERT_EQ(testState1_.vec3_.get().norm(),0);
  QPD idQuat = QPD::Identity();
  ASSERT_EQ(box_minus(testState1_.qua0(0), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(testState1_.qua0(1), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(testState1_.qua1(0), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(testState1_.qua1(1), idQuat).norm(), 0.0);
  ASSERT_EQ(mtState::Identity().sca_.get(),0);
  ASSERT_EQ(mtState::Identity().vec0_.get().norm(),0);
  ASSERT_EQ(mtState::Identity().vec1_.get().norm(),0);
  ASSERT_EQ(mtState::Identity().vec2_.get().norm(),0);
  ASSERT_EQ(mtState::Identity().vec3_.get().norm(),0);
  ASSERT_EQ(box_minus(mtState::Identity().qua0(0), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(mtState::Identity().qua0(1), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(mtState::Identity().qua1(0), idQuat).norm(), 0.0);
  ASSERT_EQ(box_minus(mtState::Identity().qua1(1), idQuat).norm(), 0.0);
}

// Test plus and minus
TEST_F(StateTesting, plusAndMinus) {
  testState2_.boxMinus(testState1_,difVec_);
  unsigned int index=0;
  ASSERT_EQ(difVec_(index),testScalar2_-testScalar1_);
  index ++;
  for(int i=0;i<4;i++){
    ASSERT_EQ((difVec_.block<3,1>(index,0)-(testVector2_[i]-testVector1_[i])).norm(),0);
    index = index + 3;
  }
  for(int i=0;i<4;i++){
    ASSERT_EQ((difVec_.block<3,1>(index,0)-box_minus(testQuat2_[i],testQuat1_[i])).norm(),0);
    index = index + 3;
  }
  testState1_.boxPlus(difVec_,testState2_);
  ASSERT_NEAR(testState2_.sca_.get(),testScalar2_,1e-10);
  ASSERT_NEAR((testState2_.vec0_.get()-testVector2_[0]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec1_.get()-testVector2_[1]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec2_.get()-testVector2_[2]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec3_.get()-testVector2_[3]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua0(0),testQuat2_[0]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua0(1),testQuat2_[1]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua1(0),testQuat2_[2]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua1(1),testQuat2_[3]).norm(),0,1e-10);
}

// Test minus Jacobian
TEST_F(StateTesting, minusJac) {
  testState2_.boxMinus(testState1_,difVecOut1_);
  double d = 0.00001;
  testState2_.boxMinusJac(testState1_,covMat_);
  for(unsigned int i=0;i<mtState::D_;i++){
    difVecIn_.setZero();
    difVecIn_(i) = d;
    testState2_.boxPlus(difVecIn_,testState3_);
    testState3_.boxMinus(testState1_,difVecOut2_);
    difVecOut2_ -= difVecOut1_;
    difVecOut2_ = difVecOut2_/d;
    ASSERT_NEAR((covMat_.col(i)-difVecOut2_).norm(),0.0,1e-5);
  }
}

// Test getValue, getId, getElementId
TEST_F(StateTesting, accessors) {
  ASSERT_NEAR(testState1_.sca_.get(),testScalar1_,1e-10);
  ASSERT_NEAR((testState1_.vec0_.get()-testVector1_[0]).norm(),0,1e-10);
  ASSERT_NEAR((testState1_.vec1_.get()-testVector1_[1]).norm(),0,1e-10);
  ASSERT_NEAR((testState1_.vec2_.get()-testVector1_[2]).norm(),0,1e-10);
  ASSERT_NEAR((testState1_.vec3_.get()-testVector1_[3]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState1_.qua0(0), testQuat1_[0]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState1_.qua0(1), testQuat1_[1]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState1_.qua1(0), testQuat1_[2]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState1_.qua1(1), testQuat1_[3]).norm(),0,1e-10);

  ASSERT_TRUE(testState1_.sca_idx_ == 0);
  ASSERT_TRUE(testState1_.vec0_idx_ == 1);
  ASSERT_TRUE(testState1_.vec1_idx_ == 4);
  ASSERT_TRUE(testState1_.vec2_idx_ == 7);
  ASSERT_TRUE(testState1_.vec3_idx_ == 10);
  ASSERT_TRUE(testState1_.qua0_idx(0) == 13);
  ASSERT_TRUE(testState1_.qua0_idx(1) == 16);
  ASSERT_TRUE(testState1_.qua1_idx(0) == 19);
  ASSERT_TRUE(testState1_.qua1_idx(1) == 22);
//
//  ASSERT_TRUE(testState1_.getElementId(0) == _sca);
//  ASSERT_TRUE(testState1_.getElementId(1) == _vec0);
//  ASSERT_TRUE(testState1_.getElementId(4) == _vec1);
//  ASSERT_TRUE(testState1_.getElementId(7) == _vec2);
//  ASSERT_TRUE(testState1_.getElementId(10) == _vec3);
//  ASSERT_TRUE(testState1_.getElementId(13) == _qua0);
//  ASSERT_TRUE(testState1_.getElementId(16) == _qua0);
//  ASSERT_TRUE(testState1_.getElementId(19) == _qua1);
//  ASSERT_TRUE(testState1_.getElementId(22) == _qua1);
//  std::cout << "CHECK (ERROR: Exceeded state size)" << std::endl;
//  ASSERT_TRUE(testState1_.getElementId(25) == _aux+1); // Exceeds state size
}

// Test operator=
TEST_F(StateTesting, operatorEQ) {
  testState2_ = testState1_;
  ASSERT_NEAR(testState2_.sca_.get(),testScalar1_,1e-10);
  ASSERT_NEAR((testState2_.vec0_.get()-testVector1_[0]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec1_.get()-testVector1_[1]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec2_.get()-testVector1_[2]).norm(),0,1e-10);
  ASSERT_NEAR((testState2_.vec3_.get()-testVector1_[3]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua0(0),testQuat1_[0]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua0(1),testQuat1_[1]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua1(0),testQuat1_[2]).norm(),0,1e-10);
  ASSERT_NEAR(box_minus(testState2_.qua1(1),testQuat1_[3]).norm(),0,1e-10);
}

// Test Constness
TEST_F(StateTesting, Constness) {
  const mtState testState1(testState1_);
  ASSERT_NEAR(box_minus(testState1.qua0(0),testQuat1_[0]).norm(),0,1e-10);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
