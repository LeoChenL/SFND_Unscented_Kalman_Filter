#include "ukf.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::endl;
using std::cout;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    // use_laser_ = false;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    // use_radar_ = false;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 5.0;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.8;

    /**
     * DO NOT MODIFY measurement noise values below.
     * These are provided by the sensor manufacturer.
     */

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

    /**
     * End DO NOT MODIFY section for measurement noise values
     */

    // Initially set to false, set to true in first call of ProcessMeasurement
    is_initialized_ = false;
    // State dimension
    n_x_ = 5;
    // Augmented state dimension
    n_aug_ = 7;
    // measurement dimension
    n_z_radar_ = 3;
    n_z_lidar_ = 2;
    // Sigma point spreading parameter
    lambda_ = 3-n_aug_;
    // NIS
    NIS_radar_ = 0;
    // timestamp
    previous_timestamp_ = 0;


    // Prediction
    // initial state vector
    x_ = VectorXd::Zero(n_x_);

    // initial covariance matrix
    P_ = MatrixXd::Zero(n_x_, n_x_);

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

    // Weights of sigma points
    weights_ = VectorXd::Zero(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
        double weight = 0.5 / (lambda_ + n_aug_);
        weights_(i) = weight;
    }


    // Radar Measurement
    //create matrix for sigma points in measurement space
    Zsig_ = MatrixXd::Zero(n_z_radar_, 2 * n_aug_ + 1);

    //mean predicted measurement
    z_pred_ = VectorXd::Zero(n_z_radar_);

    //measurement covariance matrix S
    S_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);

    //measurement noise covariance matrix
    R_radar_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);
    R_radar_ << std_radr_ * std_radr_, 0                         , 0,
                0                    , std_radphi_ * std_radphi_ , 0,
                0                    , 0                         , std_radrd_*std_radrd_;

    // Lidar Measurement
    // lidar measurement matrix
    H_lidar_ = MatrixXd(2,5);
    H_lidar_ << 1, 0, 0, 0, 0,
  			        0, 1, 0, 0, 0;

    R_lidar_ = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
    R_lidar_ << std_laspx_ * std_laspx_, 0 ,
                0                      , std_laspy_ * std_laspy_ ;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
   cout << endl <<"ProcessMeasurement " << endl;
   // check availability of radar and lidar
   if ( (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) ||
      (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) ){
        return;
      }

   /*****************************************************************************
    * Initialization
    * init the state x.
    * init the process covariance matrix P.
    * init timestamp
    ****************************************************************************/
   if (!is_initialized_){
 		// init state x_
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      cout << "KF Initialization Lidar" << endl;
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      cout << "KF Initialization Radar" << endl;
      //         x
      //         ^phi/
      //         |  /
      //         | /
      // y<-------/
      double rho, phi; //rho_dot;
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      // rho_dot = meas_package.raw_measurements_[2];
      x_ << rho * cos(phi),
            rho * sin(phi),
            0, 0, 0;
    }

    // init covariance matrix P_
    P_ << 1, 0, 0, 		0,    0,
  			  0, 1, 0, 		0,    0,
  			  0, 0, 1000, 0,    0,
  			  0, 0, 0, 		1000, 0,
          0, 0, 0,    0,    1000;

    // init timestamp
 		previous_timestamp_ = meas_package.timestamp_;
 		is_initialized_ = true;
 		return;
 	}


 	//compute the time elapsed between the current and previous measurements
 	double dt = static_cast<double>(
    (meas_package.timestamp_ - previous_timestamp_) / 1000000.0 );	//dt - expressed in seconds
 	previous_timestamp_ = meas_package.timestamp_;

  // Improve the stability by subdividing the prediction step for large dt’s into incremental updates
  // Otherwise Cholesky decomposition may fail
  while (dt > 0.1) {
      constexpr double interval_t = 0.05;
      Prediction(interval_t);
      dt -= interval_t;
  }


  /*****************************************************************************
   * Prediction Step
   * Predict the state next step state x
   * Predict the process noise covariance matrix P
   ****************************************************************************/
  std::cout << "Prediction Step " << std::endl;
  Prediction(dt);


  /*****************************************************************************
   * Measurement Step
   * Update the state and covariance matrices.
   ****************************************************************************/
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
       std::cout << "Update Step Radar" << std::endl;
       PredictRadarMeasurement();
       UpdateRadarState(meas_package.raw_measurements_);

       // update NIS
       NIS_radar_ = (meas_package.raw_measurements_-z_pred_).transpose()*
                    S_.inverse()*(meas_package.raw_measurements_-z_pred_);
       std::cout << "NIS Radar = " << NIS_radar_ << std::endl;
   }
   else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
       std::cout << "Update Step Lidar" << std::endl;
       UpdateLidarState(meas_package.raw_measurements_);
   }
}// End ProcessMeasurement


void UKF::Prediction(double dt) {
    MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

    // generate augmented sigma points
    GenerateAugmentedSigmaPoints(Xsig_aug);
    // calculate predicted sigma points
    SigmaPointPrediction(Xsig_aug, dt);
    // calculate predicted state and covariance
    PredictMeanAndCovariance();
}


void UKF::GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug){
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // create augmented mean state
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  MatrixXd Q(2,2);
  Q << pow(std_a_,2), 0,
       0,             pow(std_yawdd_,2);
  P_aug << P_, MatrixXd::Zero(P_.rows(),Q.cols()),
           MatrixXd::Zero(Q.rows(),P_.cols()), Q;

  // generate augmented sigma points
  MatrixXd A_aug = P_aug.llt().matrixL();
  // MatrixXd A_aug = P_aug.ldlt().matrixL();
  // MatrixXd D = P_aug.ldlt().vectorD().asDiagonal();
  // std::cout << "A_aug:" << std::endl << A_aug << std::endl;
  // std::cout << "D:" << std::endl << D << std::endl;
  // std::cout << "At*A:" << std::endl << A_aug *D*A_aug.transpose()<< std::endl;
  MatrixXd Xsig_aug_block1 = x_aug.replicate(1,n_aug_) +
                             sqrt(lambda_+n_aug_)*A_aug;
  MatrixXd Xsig_aug_block2 = x_aug.replicate(1,n_aug_) -
                             sqrt(lambda_+n_aug_)*A_aug;

  Xsig_aug << x_aug, Xsig_aug_block1, Xsig_aug_block2;
} // End AugmentedSigmaPoints


void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug, double dt){
  VectorXd v_determined(n_x_);
  VectorXd v_stochasitc(n_x_);

  for (int i=0;i<Xsig_aug.cols();++i){
    double p_x     = Xsig_aug(0,i);
    double p_y     = Xsig_aug(1,i);
    double v       = Xsig_aug(2,i);
    double yaw     = Xsig_aug(3,i);
    double yawd    = Xsig_aug(4,i);
    double mu_a    = Xsig_aug(5,i);
    double mu_yawd = Xsig_aug(6,i);

    if(fabs(yawd) > 0.001){
      // check zero denominator
      v_determined << v/yawd*( sin(yaw+yawd*dt)-sin(yaw) ),
                      v/yawd*(-cos(yaw+yawd*dt)+cos(yaw) ),
                      0,
                      yawd*dt,
                      0;

    } else {
      v_determined << v*cos(yaw)*dt,
                      v*sin(yaw)*dt,
                      0,
                      yawd*dt,
                      0;
    }

    v_stochasitc << 0.5*pow(dt,2)*cos(yaw)*mu_a,
                    0.5*pow(dt,2)*sin(yaw)*mu_a,
                    dt*mu_a,
                    0.5*pow(dt,2)*mu_yawd,
                    dt*mu_yawd;

    Xsig_pred_.col(i) = Xsig_aug.block(0,i,n_x_,1) + v_determined + v_stochasitc;
  }
} // End SigmaPointPrediction


void UKF::PredictMeanAndCovariance(){
  // predict state mean
  x_.setZero();
  x_ = Xsig_pred_ * weights_;

  // predict state covariance matrix
  P_.setZero();
  for (int i=0;i<2*n_aug_+1;++i){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization to keep the angle in -2pi~2pi
    while(x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;}
    while(x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;}

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}// End PredictMeanAndCovariance


void UKF::PredictRadarMeasurement(){
 // transform sigma points into measurement space
 z_pred_.setZero();
 for (int i=0;i<2*n_aug_+1;++i){
   double px   = Xsig_pred_(0,i);
   double py   = Xsig_pred_(1,i);
   double v    = Xsig_pred_(2,i);
   double yaw  = Xsig_pred_(3,i);

   VectorXd zsig(n_z_radar_);
   zsig << sqrt(pow(px,2)+pow(py,2)),
           atan2(py,px),
           ( px*cos(yaw)+ py*sin(yaw) )*v/sqrt(pow(px,2)+pow(py,2));

   // calculate mean predicted measurement
   z_pred_ += weights_(i) * zsig;
   Zsig_.col(i) = zsig;
 }

 // calculate innovation covariance matrix S
 S_.setZero();
 for (int i = 0;i < 2*n_aug_+1;++i){
   VectorXd z_diff(n_z_radar_);
   z_diff = Zsig_.col(i)-z_pred_;

   // angle normalization
   while(z_diff(1)>M_PI){z_diff(1)-=2.*M_PI;}
   while(z_diff(1)<-M_PI){z_diff(1)+=2.*M_PI;}

   S_ += weights_(i)*z_diff*z_diff.transpose();
 }
 S_ += R_radar_;
} // End PredictRadarMeasurement


void UKF::UpdateRadarState(const VectorXd& z){
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_radar_);

  // calculate cross correlation matrix
  for (int i=0;i<2 * n_aug_ + 1;++i){
    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    VectorXd z_diff = Zsig_.col(i)-z_pred_;
    // angle normalization
    while(x_diff(3)> M_PI) {x_diff(3)-= 2*M_PI;};
    while(x_diff(3)<-M_PI) {x_diff(3)+= 2*M_PI;};
    while(z_diff(1)> M_PI) {z_diff(1)-= 2*M_PI;};
    while(z_diff(1)<-M_PI) {z_diff(1)+= 2*M_PI;};

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc*S_.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z-z_pred_;
  while(z_diff(1)> M_PI) {z_diff(1)-= 2*M_PI;};
  while(z_diff(1)<-M_PI) {z_diff(1)+= 2*M_PI;};

  x_ += K*(z_diff);
  P_ -= K*S_*K.transpose();
} // End UpdateRadarState


void UKF::UpdateLidarState(const VectorXd& z) {
   // kalman gain
   MatrixXd Ht = H_lidar_.transpose();
   MatrixXd S = H_lidar_ * P_ * Ht + R_lidar_;
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si;

   // update state mean and covariance matrix
   MatrixXd I = MatrixXd::Identity(n_x_, n_x_);

 	 x_ = x_ + K * (z - H_lidar_ * x_);
 	 P_ = (I - K * H_lidar_) * P_;
}// End UpdateLidar
