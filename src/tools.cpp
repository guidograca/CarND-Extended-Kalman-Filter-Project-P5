#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.   */
  
  //Initializing variable rmse
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //Checking sizes are equal and non zero
  if (!(estimations.size()==ground_truth.size())&&(ground_truth.size()!=0)){
      std::cout<<"Data size is not equal or zero.";
      return rmse;
  }
  
  // Accumulating squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = (residual.array() * residual.array());
    rmse =   rmse + residual;
  }

  // Calculating the mean
  rmse=rmse/estimations.size();
  
  
  // Calculating the squared root
  rmse=rmse.array().sqrt();
  
  
  // Return the result
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian Matrix.
   */
  
  // Initializing variable Hj
  MatrixXd Hj(3,4);
  
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Coefficients
  double power = px*px+py*py;
  double square = pow(power,0.5);
  double crossed = vx*py-vy*px;
  double divider = power*square;

  
  // check division by zero
  if(power<0.0001){
      std::cout<<"Division by 0";
      power=0.0001;
  }

  // Compute the Jacobian matrix
  Hj<<px/square,py/square,0,0,
      -py/power, px/power,0,0,
      py*crossed/divider,-px*crossed/divider,px/square,py/square;
  return Hj;  
}
