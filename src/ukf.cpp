#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.4; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; // 30;
  
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
  time_us_ = -1.0;

  // State dimension
  n_x_ = x_.rows();

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Sigma Points Weights
  weights_ = VectorXd(2*this->n_aug_+1);

  //set weights
  for(int i = 0; i < weights_.rows(); i++) {
    weights_(i) = (i == 0) ? lambda_ / (lambda_ + n_aug_) : 1 / (2*(lambda_ + n_aug_));
  }

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

  if(!is_initialized_)
  {
    Tools tools;

    // Set the timestamp as initial timestamp
    this->time_us_ = meas_package.timestamp_;

    // Intialize the state with first measurement
    float px = 0, py = 0, v = 0;
    if (MeasurementPackage::RADAR == meas_package.sensor_type_)
    {
      float r = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float r_dot = meas_package.raw_measurements_[2];

      px = r * cos(phi);
      py = r * sin(phi);
      v = r_dot;

    } else if (MeasurementPackage::LASER == meas_package.sensor_type_) {
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
    }

    x_ << px , py , v, 0, 0;
    P_ <<   1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 100, 0, 0,
        0, 0, 0, 100, 0,
        0, 0, 0, 0, 1;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // Compute delta time
  double delta_t = (meas_package.timestamp_ - this->time_us_) / 1000000.0;

  // Check if this is radar sensor data
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && this->use_radar_)
  {
    // Do the Prediction Step
    this->Prediction(delta_t);

    // Do the Update Step
    this->UpdateRadar(meas_package);

    // This was the last measurement, so remember this as the last timestamp
    this->time_us_ = meas_package.timestamp_;
  }

  // Check if this is laser sensor data
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && this->use_laser_)
  {
    // Do the Prediction Step
    // this->Prediction(delta_t);

    // This was the last measurement, so remember this as the last timestamp
    // this->time_us_ = meas_package.timestamp_;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Compute the Augmented Sigma Points
  this->AugmentedSigmaPoints();

  // Prediction of Sigma Points
  this->SigmaPointPrediction(delta_t);

  // Prediction of Mean and Covariance
  this->PredictMeanAndCovariance(&this->x_, &this->P_);
}

void UKF::AugmentedSigmaPoints() {
  // Get the reference to the state vector
  const VectorXd& x = this->x_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(this->n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(this->n_aug_, this->n_aug_);

  //create sigma point matrix
  this->Xsig_ = MatrixXd(this->n_aug_, 2 * this->n_aug_ + 1);

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(x.rows()) = x;

  MatrixXd Q = MatrixXd(2, 2);
  Q << (this->std_a_ * this->std_a_), 0, 0, (this->std_yawdd_ * this->std_yawdd_);

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(this->P_.rows(), this->P_.cols()) = this->P_;
  P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  double lambda_sqrt = std::sqrt(this->lambda_ + this->n_aug_);

  // calculate sqrt((lambda + n_aug) * P)
  A *= lambda_sqrt;

  // Adding up and subtracting mean
  MatrixXd B(A);
  MatrixXd C(A);
  for(int c = 0; c < A.cols(); c ++) {
      B.col(c) += x_aug;
      C.col(c) = x_aug - C.col(c);
  }

  // Push it all into one matrix
  this->Xsig_ << x_aug, B, C;
}


void UKF::SigmaPointPrediction(const double & delta_t) {
  //create matrix with predicted sigma points as columns
  this->Xsig_pred_ = MatrixXd(this->n_x_, 2 * this->n_aug_ + 1);

  double delta_t2 = (delta_t * delta_t);

  for (int col = 0; col < this->Xsig_.cols(); col++) {
    //predict sigma points
    VectorXd x_aug = VectorXd(this->Xsig_.col(col));
    VectorXd x_pred = VectorXd(this->n_x_);

    double &ypsi = x_aug(2);
    double &psi = x_aug(3);
    double &psi_dot = x_aug(4);
    double &nu_a = x_aug(5);
    double &nu_psi_dd = x_aug(6);

    VectorXd f = VectorXd(this->n_x_);
    VectorXd nu_k = VectorXd(this->n_x_); // Process noise vector nu k

    double epsilon = 1e-5;
    //avoid division by zero
    if ((psi_dot < epsilon) && (psi_dot > -epsilon)) // nearly zero
    {
       f << ypsi * std::cos(psi) * delta_t,
            ypsi * std::sin(psi) * delta_t,
            0,
            psi_dot * delta_t,
            0;

       nu_k <<
            0.5 * delta_t2 * std::cos(psi) * nu_a,
            0.5 * delta_t2 * std::sin(psi) * nu_a,
            delta_t * nu_a,
            0.5 * delta_t2 * nu_psi_dd,
            delta_t * nu_psi_dd;
    }
    else
    {
       f << (ypsi / psi_dot)*(std::sin(psi + psi_dot * delta_t) - sin(psi)),
            (ypsi / psi_dot)*(-std::cos(psi + psi_dot * delta_t) + cos(psi)),
            0,
            psi_dot * delta_t,
            0;
       nu_k <<
            0.5 * delta_t2 * std::cos(psi) * nu_a,
            0.5 * delta_t2 * std::sin(psi) * nu_a,
            delta_t * nu_a,
            0.5 * delta_t2 * nu_psi_dd,
            delta_t * nu_psi_dd;
    }

    //write predicted sigma points into right column
    x_pred = x_aug.topRows(this->n_x_) + f + nu_k;
    this->Xsig_pred_.col(col) = x_pred;
  }
}


void UKF::PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred) {

  //create vector for predicted state
  VectorXd x = VectorXd(this->n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(this->n_x_, this->n_x_);

  //predict state mean
  x.fill(0.0);
  for(int i = 0; i < this->Xsig_pred_.cols(); i++) {
    x += this->weights_(i) * this->Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P.fill(0.0);

  VectorXd tmp;
  for(int i = 0; i < this->Xsig_pred_.cols(); i++) {
    // temporary buffer
    tmp = this->Xsig_pred_.col(i) - x;

    // Angle normalization
    while (tmp(3) > M_PI) {
        tmp(3) -= 2.0 * M_PI;
    }
    while (tmp(3) < -M_PI) {
        tmp(3) += 2.0 * M_PI;
    }

    P += this->weights_(i) * (tmp * tmp.transpose());
  }

  *x_pred = x;
  *P_pred = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {

  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * this->n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  VectorXd wi = VectorXd(3);
  wi << this->std_radr_, this->std_radphi_, this->std_radrd_;

  MatrixXd R = MatrixXd(3, 3);
  R << (this->std_radr_ * this->std_radr_), 0, 0,
       0, (this->std_radphi_ * this->std_radphi_), 0,
       0, 0, (this->std_radrd_ * this->std_radrd_);

  // For every predicted sigma point, apply the measurement model
  // to transform sigma points into measurement space
  VectorXd transformed(n_z);
  for(int i = 0; i < this->Xsig_pred_.cols(); i++) {

      const VectorXd &Xsig_pred_c = this->Xsig_pred_.col(i);

      const double &px   = Xsig_pred_c(0);
      const double &py   = Xsig_pred_c(1);
      const double &ypsi = Xsig_pred_c(2);
      const double &psi  = Xsig_pred_c(3);

      transformed(0) = std::sqrt(px*px + py*py);
      transformed(1) = std::atan2(py, px);
      transformed(2) = ((px*std::cos(psi)*ypsi) + (py*std::sin(psi)*ypsi)) / transformed(0);

      Zsig.col(i) = transformed;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < Zsig.cols(); i++) {
    z_pred += this->weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  VectorXd tmp(n_z);
  for(int i = 0; i < Zsig.cols(); i++) {
     tmp = Zsig.col(i) - z_pred;

     //angle normalization
     while (tmp(1)> M_PI) tmp(1)-=2.*M_PI;
     while (tmp(1)<-M_PI) tmp(1)+=2.*M_PI;

     S += this->weights_(i) * (tmp * tmp.transpose());
  }
  S += R;

  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}


void UKF::UpdateStateRadar(VectorXd& z, VectorXd& z_pred, MatrixXd& S, MatrixXd& Zsig, VectorXd* x_out, MatrixXd* P_out, MatrixXd* NIS_R_out)  {

  int n_z = 3;

  // copy of the predicted state mean
  VectorXd x = this->x_;

  // copy of the predicted covariance matrix
  MatrixXd P = this->P_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(this->n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < 2 * this->n_aug_ + 1; i++) {

      VectorXd tmp_z_pred_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (tmp_z_pred_diff(1)> M_PI) tmp_z_pred_diff(1)-=2.*M_PI;
      while (tmp_z_pred_diff(1)<-M_PI) tmp_z_pred_diff(1)+=2.*M_PI;

      VectorXd tmp_x_diff = this->Xsig_pred_.col(i) - x;

      //angle normalization
      while (tmp_x_diff(1)> M_PI) tmp_x_diff(1)-=2.*M_PI;
      while (tmp_x_diff(1)<-M_PI) tmp_x_diff(1)+=2.*M_PI;

      Tc += this->weights_(i) * ((tmp_x_diff) * (tmp_z_pred_diff).transpose());
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Tmp
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K * S * K.transpose();

  std::cout << x << std::endl;
  std::cout << P << std::endl;

  *x_out = x;
  *P_out = P;
  *NIS_R_out = z_diff.transpose() * S.inverse() * z_diff;
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
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z_pred;
  MatrixXd S, Zsig;
  MatrixXd NIS;

  this->PredictRadarMeasurement(&z_pred, &S, &Zsig);
  this->UpdateStateRadar(meas_package.raw_measurements_, z_pred, S, Zsig, &this->x_, &this->P_, &NIS);
}
