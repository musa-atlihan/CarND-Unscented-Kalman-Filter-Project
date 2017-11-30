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
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // Parameters above this line are scaffolding, do not modify

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // model is not initialized yet, so set false
  is_initialized_ = false;

  // set time in miliseconds
  time_us_ = 0.f;

  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = 7;

  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // the matrix for predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //initial augmented state vector
  x_aug_ = VectorXd(n_aug_);

  // augmented covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  NIS_radar_ = 0;
  NIS_laser_ = 0;
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

  if (! is_initialized_)
  {
    // init state vector
    x_ << 0, 0, 0, 0, 0;

    // init covariance matrix
    P_ << 1.f, 0, 0, 0, 0,
          0, 1.f, 0, 0, 0,
          0, 0, 1.f, 0, 0,
          0, 0, 0, 1.f, 0,
          0, 0, 0, 0, 1.f;

    // init timestamp
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR  && use_radar_)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      float vx = rho_dot*cos(phi);
      float vy = rho_dot*sin(phi);
      // init state vector
      x_ << rho*cos(phi), rho*sin(phi), sqrt(vx*vx + vy*vy), 0.f, 0.f;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER  && use_laser_)
    {
      // init state vector
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.f, 0.f, 0.f;
    }

    // handle if Px and Py are zero
    if (fabs(x_[0]) < 0.001 && fabs(x_[1]) < 0.001)
    {
        x_[0] = 0.001;
        x_[1] = 0.001;
    }

    // init augmented state vector
    x_aug_.head(5) = x_;
    x_aug_[5] = 0.f;
    x_aug_[6] = 0.f;

    // init augmented state covariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_,n_x_) = P_;
    P_aug_(5,5) = std_a_*std_a_;
    P_aug_(6,6) = std_yawdd_*std_yawdd_;

    // init weights
    weights_(0) = lambda_/(lambda_ + n_aug_);
    for (int i=1; i<2*n_aug_+1; i++)
    {
      weights_(i) = 0.5/(n_aug_+lambda_);
    }

    // prev time
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return; // done initializing
  }

  // PREDICTION STEP //
  // delta t measured in seconds
  float dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.f;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  // UPDATE STEP //
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
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

  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_[5] = 0.f;
  x_aug_[6] = 0.f;

  // create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = std_a_*std_a_;
  P_aug_(6, 6) = std_yawdd_*std_yawdd_;

  // create square root matrix of P_aug_
  MatrixXd A = P_aug_.llt().matrixL();

  // augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  // generate augmented sigma points
  Xsig_aug.col(0)  = x_aug_;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug_ + sqrt(lambda_+n_aug_)*A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_)*A.col(i);
  }

  // predict sigma points
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd*( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd*( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // set weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++)
  {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // predicted state mean
  // iterate over sigma points
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  // iterate over sigma points
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i)*x_diff*x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // set measurement dimension (px and py)
  int n_z = 2;

  // get measurements
  VectorXd z = meas_package.raw_measurements_;

  // init Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
  Zsig.fill(0.0);

  // init S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // init measurement prediction vector
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    // sigma point predictions in process space
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    // sigma point predictions in measurement space
    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  // mean predicted measurement
  z_pred = Zsig*weights_;

  // measurement covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  // create measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  // update x_ and P_
  MatrixXd K = Tc*S.inverse();
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  //NIS Lidar Update
  NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;
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

  //set measurement dimension (rho, phi, rho_dot)
  int n_z = 3;

  // get measurements
  VectorXd z = meas_package.raw_measurements_;

  // init Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
  Zsig.fill(0.0);

  // init S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  // init measurement prediction vector
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // if zeros
    if (fabs(px) < 0.001)
    {
     px = 0.001;
    }
    if (fabs(py) < 0.001)
    {
     py = 0.001;
    }

    // rho
    Zsig(0,i) = sqrt(px*px + py*py);
    // phi
    Zsig(1,i) = atan2(py,px);
    // r_dot
    Zsig(2,i) = (px*v1 + py*v2 ) / sqrt(px*px + py*py);
  }

  // mean predicted measurement
  z_pred = Zsig*weights_;

  // measurement covariance matrix S
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  // create measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S = S + R;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  // update x_ and P_
  MatrixXd K = Tc*S.inverse();
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  // NIS Update
  NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;
}
