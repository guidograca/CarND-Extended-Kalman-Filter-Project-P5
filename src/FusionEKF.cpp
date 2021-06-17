#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

	// Creating the H_laser_ and P matrixes: 
	H_laser_ << 1,0,0,0,
  				0,1,0,0;

	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_<< 1,0,0,0,
  			  0,1,0,0,
  			0,0,10000,0,
  			0,0,0,10000;     
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
 
  
  /*==========================================================================


   			* * * Initialization * * *
   
   
   ==========================================================================*/
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // First measurement
    cout << "EKF: First measurement" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
	
	previous_timestamp_ = measurement_pack.timestamp_;

    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      cout << "First mesuarement is RADAR" << endl;
      double rho, theta;
      double px,py;
      rho=measurement_pack.raw_measurements_[0];
      theta = measurement_pack.raw_measurements_[1];
      px = rho*cos(theta);
      py = rho*sin(theta);
      if (px<0.0001){
        px=0.0001;
      }
      if (py<0.0001){
        py=0.0001;
      }
      
	  ekf_.x_ << px, 
              	 py, 
                 0, 
                 0;
}
      
      
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      cout << "First mesuarement is LASER" << endl;
	  ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;
    }      

    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
        }
  
  /*==========================================================================


   			* * * Prediction * * *
   
   
   ==========================================================================*/

  // Measuring time 
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Updating F matrix
  VectorXd z_laser(2), z_radar(3);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;   

  //    Calculating Q
  double noise_ax = 9;
  double noise_ay = 9;
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  
  // Running the prediction
  ekf_.Predict();

  
  /*==========================================================================


   			* * * Update * * *
   
   
   ==========================================================================*/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
	z_radar << measurement_pack.raw_measurements_[0], 
    		   measurement_pack.raw_measurements_[1],
    		   measurement_pack.raw_measurements_[2]; 
    if (abs(tan(measurement_pack.raw_measurements_[1]))<10){
	  
      //This is used in order to filter out angles too close to 90°. It is better to rely only on radar measurements
      cout << "Next mesuarement is RADAR" <<endl;
      Hj_=tools.CalculateJacobian(ekf_.x_);
      ekf_.H_=Hj_;
      ekf_.R_=R_radar_;

      ekf_.UpdateEKF(z_radar);
    }
  } else {
    cout << "Next mesuarement is LASER" << endl;
	z_laser <<measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1];
    ekf_.H_=H_laser_;
    ekf_.R_=R_laser_;
    ekf_.Update(z_laser);

  }

  // Print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
