// src/ukf.cpp
#include "ukf.h"
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

UKF::UKF() {
  // sensor usage
  use_laser_ = true;
  use_radar_ = true;

  // initial state vector
  n_x_ = 5;
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // process noise
  std_a_ = 1.5;
  std_yawdd_ = 0.5;

  // measurement noise
  std_laspx_ = 0.15;
  std_laspy_ = 0.15;
  std_radr_   = 0.3;
  std_radphi_ = 0.03;
  std_radrd_  = 0.3;

  // augmented dimension
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // sigma point prediction matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  is_initialized_ = false;
  time_us_ = 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0.0, 0.0, 0.0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, 0.0, 0.0;
    }

    P_ = MatrixXd::Identity(n_x_, n_x_);
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // compute delta t
  double dt = (meas_package.timestamp_ - time_us_) / 1e6;
  time_us_ = meas_package.timestamp_;

  // predict
  Prediction(dt);

  // update
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  // create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // create augmented covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // square root
  MatrixXd L = P_aug.llt().matrixL();

  // sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i+1)        = x_aug + sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_n_aug * L.col(i);
  }

  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_p, py_p;

    if (fabs(yawd) > 1e-6) {
      px_p = px + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = py + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
    } else {
      px_p = px + v * cos(yaw) * delta_t;
      py_p = py + v * sin(yaw) * delta_t;
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p  += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;
  VectorXd z = meas_package.raw_measurements_;

  MatrixXd H = MatrixXd::Zero(n_z, n_x_);
  H(0,0) = 1.0;
  H(1,1) = 1.0;

  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R(0,0) = std_laspx_ * std_laspx_;
  R(1,1) = std_laspy_ * std_laspy_;

  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  // update
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K*H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    double rho = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double rho_dot = (rho > 1e-6) ? (px*vx + py*vy)/rho : 0.0;

    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_dot;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // measurement covariance
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  S += R;

  // cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  // update
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
