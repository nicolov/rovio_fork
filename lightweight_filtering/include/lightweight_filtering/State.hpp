/*
 * State.hpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Bloeschm
 */

#ifndef LWF_STATE_HPP_
#define LWF_STATE_HPP_

#include "lightweight_filtering/common.hpp"
#include "lightweight_filtering/PropertyHandler.hpp"
#include <random>

namespace LWF{

template<typename DERIVED, typename GET, unsigned int D, unsigned int E = D>
class ElementBase{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ElementBase(){};
  static constexpr unsigned int D_ = D;
  static constexpr unsigned int E_ = E;
  typedef Eigen::Matrix<double,D_,1> mtDifVec;
  typedef GET mtGet;
  std::string name_;
  static DERIVED Identity(){
    DERIVED identity;
    identity.setIdentity();
    return identity;
  }
  DERIVED& operator=(DERIVED other){
    other.swap(*this);
    return *this;
  }
  void registerCovarianceToPropertyHandler(Eigen::MatrixXd& cov, PropertyHandler* mpPropertyHandler, const std::string& str, int j){
    assert(j+D_<=cov.cols());
    for(unsigned int i=0;i<DERIVED::D_;i++){
      mpPropertyHandler->doubleRegister_.registerScalar(str + name_ + "_" + std::to_string(i), cov(j+i,j+i));
    }
  }
};

class ScalarElement: public ElementBase<ScalarElement,double,1>{
 public:
  double s_;
  ScalarElement(){}
  ScalarElement(const ScalarElement& other){
    s_ = other.s_;
  }
  void boxPlus(const mtDifVec& vecIn, ScalarElement& stateOut) const{
    stateOut.s_ = s_ + vecIn(0);
  }
  void boxMinus(const ScalarElement& stateIn, mtDifVec& vecOut) const{
    vecOut(0) = s_ - stateIn.s_;
  }
  void boxMinusJac(const ScalarElement& stateIn, MXD& matOut) const{
    matOut.setIdentity();
  }
  void print() const{
    std::cout << s_ << std::endl;
  }
  void setIdentity(){
    s_ = 0.0;
  }
  void setRandom(unsigned int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    s_ = distribution(generator);
    s++;
  }
  void fix(){
  }
  mtGet& get(unsigned int i = 0){
    assert(i==0);
    return s_;
  }
  const mtGet& get(unsigned int i = 0) const{
    assert(i==0);
    return s_;
  }
  void registerElementToPropertyHandler(PropertyHandler* mpPropertyHandler, const std::string& str){
    mpPropertyHandler->doubleRegister_.registerScalar(str + name_,s_);
  }
};

template<unsigned int N>
class VectorElement: public ElementBase<VectorElement<N>,Eigen::Matrix<double,N,1>,N>{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef ElementBase<VectorElement<N>,Eigen::Matrix<double,N,1>,N> Base;
  using typename Base::mtDifVec;
  using typename Base::mtGet;
  using Base::name_;
  static const unsigned int N_ = N;
  Eigen::Matrix<double,N_,1> v_;
  VectorElement(){}
  VectorElement(const VectorElement<N>& other){
    v_ = other.v_;
  }
  void boxPlus(const mtDifVec& vecIn, VectorElement<N>& stateOut) const{
    stateOut.v_ = v_ + vecIn;
  }
  void boxMinus(const VectorElement<N>& stateIn, mtDifVec& vecOut) const{
    vecOut = v_ - stateIn.v_;
  }
  void boxMinusJac(const VectorElement<N>& stateIn, MXD& matOut) const{
    matOut.setIdentity();
  }
  void print() const{
    std::cout << v_.transpose() << std::endl;
  }
  void setIdentity(){
    v_.setZero();
  }
  void setRandom(unsigned int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    for(unsigned int i=0;i<N_;i++){
      v_(i) = distribution(generator);
    }
    s++;
  }
  void fix(){
  }
  mtGet& get(unsigned int i = 0){
    assert(i==0);
    return v_;
  }
  const mtGet& get(unsigned int i = 0) const{
    assert(i==0);
    return v_;
  }
  void registerElementToPropertyHandler(PropertyHandler* mpPropertyHandler, const std::string& str){
    for(unsigned int i = 0;i<N_;i++){
      mpPropertyHandler->doubleRegister_.registerScalar(str + name_ + "_" + std::to_string(i), v_(i));
    }
  }
};

class QuaternionElement: public ElementBase<QuaternionElement,QPD,3>{
 public:
  QPD q_;
  QuaternionElement(){}
  QuaternionElement(const QuaternionElement& other){
    q_ = other.q_;
  }
  void boxPlus(const mtDifVec& vecIn, QuaternionElement& stateOut) const{
    stateOut.q_ = box_plus(q_, vecIn);
  }
  void boxMinus(const QuaternionElement& stateIn, mtDifVec& vecOut) const{
    vecOut = box_minus(q_, stateIn.q_);
  }
  void boxMinusJac(const QuaternionElement& stateIn, MXD& matOut) const{
    mtDifVec diff;
    boxMinus(stateIn,diff);
    matOut = Lmat(diff).inverse();
  }
  void print() const{
    std::cout << q_.matrix() << std::endl;
  }
  void setIdentity(){
    q_.setIdentity();
  }
  void setRandom(unsigned int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    q_.w() = distribution(generator);
    q_.x() = distribution(generator);
    q_.y() = distribution(generator);
    q_.z() = distribution(generator);
    fix();
    s++;
  }
  void fix(){
    q_.normalize();
  }
  mtGet& get(unsigned int i = 0){
    assert(i==0);
    return q_;
  }
  const mtGet& get(unsigned int i = 0) const{
    assert(i==0);
    return q_;
  }
  void registerElementToPropertyHandler(PropertyHandler* mpPropertyHandler, const std::string& str){
    mpPropertyHandler->doubleRegister_.registerQuaternion(str + name_, q_);
  }
};

class NormalVectorElement: public ElementBase<NormalVectorElement,NormalVectorElement,2>{
 public:
  QPD q_;
  const V3D e_x;
  const V3D e_y;
  const V3D e_z;
  NormalVectorElement(): e_x(1,0,0), e_y(0,1,0), e_z(0,0,1){}
  NormalVectorElement(const NormalVectorElement& other): e_x(1,0,0), e_y(0,1,0), e_z(0,0,1){
    q_ = other.q_;
  }
  NormalVectorElement(const V3D& vec): e_x(1,0,0), e_y(0,1,0), e_z(0,0,1){
    setFromVector(vec);
  }
  NormalVectorElement(const QPD& q): e_x(1,0,0), e_y(0,1,0), e_z(0,0,1){
    q_ = q;
  }
  V3D getVec() const{
    return q_ * e_z;
  }
  V3D getPerp1() const{
    return q_ * e_x;
  }
  V3D getPerp2() const{
    return q_ * e_y;
  }
  NormalVectorElement& operator=(const NormalVectorElement& other){
    q_ = other.q_;
    return *this;
  }
  static V3D getRotationFromTwoNormals(const V3D& a, const V3D& b, const V3D& a_perp){
    const V3D cross = a.cross(b);
    const double crossNorm = cross.norm();
    const double c = a.dot(b);
    const double angle = std::acos(c);
    if(crossNorm<1e-6){
      if(c>0){
        return cross;
      } else {
        return a_perp*M_PI;
      }
    } else {
      return cross*(angle/crossNorm);
    }
  }
  static V3D getRotationFromTwoNormals(const NormalVectorElement& a, const NormalVectorElement& b){
    return getRotationFromTwoNormals(a.getVec(),b.getVec(),a.getPerp1());
  }
  static M3D getRotationFromTwoNormalsJac(const V3D& a, const V3D& b){
    const V3D cross = a.cross(b);
    const double crossNorm = cross.norm();
    V3D crossNormalized = cross/crossNorm;
    M3D crossNormalizedSqew = gSM(crossNormalized);
    const double c = a.dot(b);
    const double angle = std::acos(c);
    if(crossNorm<1e-6){
      if(c>0){
        return -gSM(b);
      } else {
        return M3D::Zero();
      }
    } else {
      return -1/crossNorm*(crossNormalized*b.transpose()-(crossNormalizedSqew*crossNormalizedSqew*gSM(b)*angle));
    }
  }
  static M3D getRotationFromTwoNormalsJac(const NormalVectorElement& a, const NormalVectorElement& b){
    return getRotationFromTwoNormalsJac(a.getVec(),b.getVec());
  }
  void setFromVector(V3D vec){
    const double d = vec.norm();
    if(d > 1e-6){
      vec = vec/d;
      q_ = exp_map(getRotationFromTwoNormals(e_z,vec,e_x));
    } else {
      q_.setIdentity();
    }
  }
  NormalVectorElement rotated(const QPD& q) const{
    return NormalVectorElement(q*q_);
  }
  NormalVectorElement inverted() const{
    QPD q = exp_map(M_PI * getPerp1());
    return NormalVectorElement(q*q_);
  }
  void boxPlus(const mtDifVec& vecIn, NormalVectorElement& stateOut) const{
    QPD q = exp_map(vecIn(0)*getPerp1()+vecIn(1)*getPerp2());
    stateOut.q_ = q*q_;
  }
  void boxMinus(const NormalVectorElement& stateIn, mtDifVec& vecOut) const{
    vecOut = stateIn.getN().transpose()*getRotationFromTwoNormals(stateIn,*this);
  }
  void boxMinusJac(const NormalVectorElement& stateIn, MXD& matOut) const{
    matOut = -stateIn.getN().transpose()*getRotationFromTwoNormalsJac(*this,stateIn)*this->getM();
  }
  void print() const{
    std::cout << getVec().transpose() << std::endl;
  }
  void setIdentity(){
    q_.setIdentity();
  }
  void setRandom(unsigned int& s){
    std::default_random_engine generator (s);
    std::normal_distribution<double> distribution (0.0,1.0);
    q_.w() = distribution(generator);
    q_.x() = distribution(generator);
    q_.y() = distribution(generator);
    q_.z() = distribution(generator);
    q_.normalize();
    s++;
  }
  void fix(){
    q_.normalize();
  }
  mtGet& get(unsigned int i = 0){
    assert(i==0);
    return *this;
  }
  const mtGet& get(unsigned int i = 0) const{
    assert(i==0);
    return *this;
  }
  void registerElementToPropertyHandler(PropertyHandler* mpPropertyHandler, const std::string& str){
    mpPropertyHandler->doubleRegister_.registerQuaternion(str + name_, q_);
  }
  Eigen::Matrix<double,3,2> getM() const {
    Eigen::Matrix<double,3,2> M;
    M.col(0) = -getPerp2();
    M.col(1) = getPerp1();
    return M;
  }
  Eigen::Matrix<double,3,2> getN() const {
    Eigen::Matrix<double,3,2> M;
    M.col(0) = getPerp1();
    M.col(1) = getPerp2();
    return M;
  }
};

template<typename Element, unsigned int M>
class ArrayElement: public ElementBase<ArrayElement<Element,M>,typename Element::mtGet,M*Element::D_,Element::D_>{
 public:
  typedef ElementBase<ArrayElement<Element,M>,typename Element::mtGet,M*Element::D_,Element::D_> Base;
  using typename Base::mtDifVec;
  using typename Base::mtGet;
  using Base::name_;
  using Base::D_;
  static const unsigned int M_ = M;
  Element array_[M_];
  mutable MXD boxMinusJacMat_;
  ArrayElement(): boxMinusJacMat_((int)Element::D_,(int)Element::D_){
    for(unsigned int i=0; i<M_;i++){
      array_[i].name_ = "";
    }
  }
  ArrayElement(const ArrayElement& other): boxMinusJacMat_((int)Element::D_,(int)Element::D_){
    for(unsigned int i=0; i<M_;i++){
      array_[i] = other.array_[i];
    }
  }
  void boxPlus(const mtDifVec& vecIn, ArrayElement& stateOut) const{
    boxPlus_(vecIn,stateOut);
  }
  template<bool b = (D_>0), typename std::enable_if<b>::type* = nullptr>
  void boxPlus_(const mtDifVec& vecIn, ArrayElement& stateOut) const{
    for(unsigned int i=0; i<M_;i++){
      array_[i].boxPlus(vecIn.template block<Element::D_,1>(Element::D_*i,0),stateOut.array_[i]);
    }
  }
  template<bool b = (D_>0), typename std::enable_if<!b>::type* = nullptr>
  void boxPlus_(const mtDifVec& vecIn, ArrayElement& stateOut) const{}
  void boxMinus(const ArrayElement& stateIn, mtDifVec& vecOut) const{
    boxMinus_(stateIn,vecOut);
  }
  template<bool b = (D_>0), typename std::enable_if<b>::type* = nullptr>
  void boxMinus_(const ArrayElement& stateIn, mtDifVec& vecOut) const{
    typename Element::mtDifVec difVec;
    for(unsigned int i=0; i<M_;i++){
      array_[i].boxMinus(stateIn.array_[i],difVec);
      vecOut.template block<Element::D_,1>(Element::D_*i,0) = difVec;
    }
  }
  template<bool b = (D_>0), typename std::enable_if<!b>::type* = nullptr>
  void boxMinus_(const ArrayElement& stateIn, mtDifVec& vecOut) const{}
  void boxMinusJac(const ArrayElement& stateIn, MXD& matOut) const{
    boxMinusJac_(stateIn,matOut);
  }
  template<bool b = (D_>0), typename std::enable_if<b>::type* = nullptr>
  void boxMinusJac_(const ArrayElement& stateIn, MXD& matOut) const{
    matOut.setZero();
    for(unsigned int i=0; i<M_;i++){
      array_[i].boxMinusJac(stateIn.array_[i],boxMinusJacMat_);
      matOut.template block<Element::D_,Element::D_>(Element::D_*i,Element::D_*i) = boxMinusJacMat_;
    }
  }
  template<bool b = (D_>0), typename std::enable_if<!b>::type* = nullptr>
  void boxMinusJac_(const ArrayElement& stateIn, MXD& matOut) const{}
  void print() const{
    for(unsigned int i=0; i<M_;i++){
      array_[i].print();
    }
  }
  void setIdentity(){
    for(unsigned int i=0; i<M_;i++){
      array_[i].setIdentity();
    }
  }
  void setRandom(unsigned int& s){
    for(unsigned int i=0; i<M_;i++){
      array_[i].setRandom(s);
    }
  }
  void fix(){
    for(unsigned int i=0; i<M_;i++){
      array_[i].fix();
    }
  }
  mtGet& get(unsigned int i){
    assert(i<M_);
    return array_[i].get();
  }
  const mtGet& get(unsigned int i) const{
    assert(i<M_);
    return array_[i].get();
  }
  void registerElementToPropertyHandler(PropertyHandler* mpPropertyHandler, const std::string& str){
    for(unsigned int i=0; i<M_;i++){
      array_[i].registerElementToPropertyHandler(mpPropertyHandler,str + name_ + "_" + std::to_string(i));
    }
  }
};

}

#endif /* LWF_STATE_HPP_ */
