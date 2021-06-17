#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  // Prediction
  x_=F_*x_;
  P_=F_*P_*F_.transpose()+Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  // Kalman Filter Equations:  
  VectorXd y =z-H_*x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
    */
 // importing the state variables
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2); 
  double vy = x_(3);   
  
  //Adjusting too small values of px and py
  if (abs(px)<0.0001){
    if (px>0){
      px=0.0001;
    } 
    else{
      px=-0.0001;
    }
  }
  if (abs(py)<0.0001){
    if (py>0){
      py=0.0001;
    } 
    else{
      py=-0.0001;
    }
  } 
  
  //Polar coordinates
  double rho = pow(px*px+py*py,0.5);
  double phi = atan2(py,px);
  double rhoDot = (px*vx + py*vy) / rho;
  

  VectorXd h(3);
  h<<rho, phi, rhoDot;
  VectorXd y = z - h;
  
  //Adjusting rho value
  while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    } else {
      y(1) += M_PI;
    }
  } 
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_* P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
