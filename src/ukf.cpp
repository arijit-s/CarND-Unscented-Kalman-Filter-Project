#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  //TODO: Tune the paramters.
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  //P Value from udacity lessons
  //  P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
  //         -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
  //         0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
  //         -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
  //         -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4.5;
  //std_a_ = 0.26;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 0.7;
  std_yawdd_ = 0.39;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set State dimension
  n_x_ = 5;
  //set augmented dimension
  n_aug_ = 7;
  //Set predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //Set weights_
  weights_ = VectorXd(2*n_aug_+1);
  
  lambda_ = 3 - n_aug_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout << "Initializing UKF: " << endl;
    x_ << 1, 1, 1, 1, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      //TODO:Tune the state initialization.
      //float v = 0.0198;
      //float yaw = 0.000345566;
      //float yaw_rate = 0.01382155;
      float v = 0;
      float yaw = 0;
      float yaw_rate = 0;

      x_ << px,py,v,yaw,yaw_rate;
      // cout << "First RADAR " << endl;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //For Laser 
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      // cout << "First LIDAR " << endl;

    }
    // cout << "Init Done " << endl;

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  
  // cout << "Starting prediction" << endl;

  Prediction(dt);

  // cout << "End prediction" << endl;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // cout << "Starting RADAR Update" << endl;
      UpdateRadar(meas_package);
      // cout << "Ending RADAR Update" << endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // cout << "Starting LIDAR Update" << endl;
      UpdateLidar(meas_package);
      // cout << "Ending LIDAR Update" << endl;
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
// cout << "Generate sigma points" << endl;
  //Generate sigma points
  MatrixXd Xsig = MatrixXd(n_x_,2*n_x_ + 1);
  MatrixXd A = P_.llt().matrixL();
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  VectorXd x_aug = VectorXd(7);
// cout << "create augmented state covariance" << endl;
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

// cout << "Sigma point prediction: "<<delta_t << endl;

  //Sigma point prediction
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
// cout << "   ->check division by 0 :"<<yawd<< endl;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        // cout << "   ->False"<<yawd<< endl;
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
        // cout << "   ->True"<<yawd<< endl;
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    // cout << "write predicted sigma point into right column"<< endl;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

// cout << "Predicted mean and covariance" << endl;

  //Predicted mean and covariance

  //create vector for predicted state
  VectorXd x_pred = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

// cout << "predicted state mean" << endl;

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

// cout << "predicted state covariance" << endl;
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
// cout << "End inside predict" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateLidarOld(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd R_laser_ = MatrixXd(2,2);
  // cout<<"Initialize z: "<<meas_package.raw_measurements_<<endl;

  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  MatrixXd H_ = MatrixXd(2,5);
  //TODO: Change parameters
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;
  // cout<<"H size and x size"<<H_.size()<<x_.size()<<endl;
  VectorXd z_pred = H_ * x_;
  // cout<<"z_pred size"<<z_pred.size()<<endl;
  VectorXd z = meas_package.raw_measurements_;

  // cout<<"Z and Z_pred sizezs: "<<z.size()<<z_pred<<endl;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

  float NIS_Lidar = y.transpose()*S.inverse()*y;
  //cout<<"NIS : "<<eta<<endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z_ = 2;
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  //add measurement noise covariance matrix
  MatrixXd R_laser_ = MatrixXd(n_z_,n_z_);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;


  z_pred.fill(0.0);
  S.fill(0.0);

  // //transform sigma points into measurement space
  // for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

  //   // extract values for better readibility
  //   double p_x = Xsig_pred_(0,i);
  //   double p_y = Xsig_pred_(1,i);
  //   //double v  =  Xsig_pred_(2,i);
  //   //double yaw = Xsig_pred_(3,i);

  //   // double v1 = cos(yaw)*v;
  //   // double v2 = sin(yaw)*v;

  //   // measurement model
  //   Zsig(0,i) = p_x;
  //   Zsig(1,i) = p_y;
  // }

  //mean predicted measurement
  for (int i=0; i < 2*n_aug_+1; i++) {
      Zsig(0,i) = Xsig_pred_(0,i);
      Zsig(1,i) = Xsig_pred_(1,i);
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }


  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R_laser_;


  //UKF Update
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = meas_package.raw_measurements_;
 
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


  float NIS_Lidar = z_diff.transpose()*S.inverse()*z_diff;
  //cout<<"NIS : "<<eta<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /*
  //Predict RADAR measurements.
  */

  //create matrix for sigma points in measurement space
  int n_z_ = 3;
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);
  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;


  z_pred.fill(0.0);
  S.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  =  Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }


  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  S = S + R;

  //UKF Update
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = meas_package.raw_measurements_;
 
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //Calculate NIS for LIDAR:
  float NIS_Radar = z_diff.transpose()*S.inverse()*z_diff;
  //cout<<"NIS : "<<eta<<endl;
 
}