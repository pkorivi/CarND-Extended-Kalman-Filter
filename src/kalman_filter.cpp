#include "kalman_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_in,
                        MatrixXd &R_radar_in, MatrixXd &Q_in){
  cout<<"kalman init"<<endl;
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_ = R_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
 x_ = F_*x_;
 MatrixXd Ft = F_.transpose();
 P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
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
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  Hj_ = tools.CalculateJacobian( x_ );
  VectorXd z_pred(3);
  float rho = sqrt( px*px + py*py );

  //Avoiding division by zero
  if( rho == 0.0)
    return;

  z_pred << rho, atan2( py, px ), ( px*vx + py*vy )/rho;

  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - z_pred;

  if( y[1] > M_PI )
    y[1] -= 2.0*M_PI;
  if( y[1] < -M_PI )
    y[1] += 2.0*M_PI;
  
  MatrixXd Hj_trans = Hj_.transpose();
  MatrixXd S = Hj_*P_*Hj_trans + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Hj_trans*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());
  P_ = ( I_ - K*Hj_ )*P_;
}
