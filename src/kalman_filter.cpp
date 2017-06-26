#define _USE_MATH_DEFINES

#include "kalman_filter.h"

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
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	CalculateUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  double rho = sqrt((px * px) + (py * py));
  double phi = atan2(py, px);
  double rhodot = ((px * vx) + (py * vy)) / rho;

  VectorXd h_x = VectorXd(3);
  h_x << rho, phi, rhodot;
  VectorXd y = z - h_x;

  // adjust phi in y to between -PI and PI
  while(y(1) > M_PI) {
    y(1) = y(1) - 2 * M_PI;
  }
  while(y(1) < -1 * M_PI) {
    y(1) = y(1) + 2 * M_PI;
  }
  CalculateUpdate(y);
}

void KalmanFilter::CalculateUpdate(const VectorXd &y) {
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