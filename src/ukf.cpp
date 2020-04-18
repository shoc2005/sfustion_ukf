#include "ukf.h"
#include "Eigen/Dense"


using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  time_us_ = 0;

  // initial state vector
  x_ = VectorXd(5);
  LidarM_ = VectorXd(2);
  RadarM_ = VectorXd(3);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   is_initialized_ = false ;
   x_.fill(0.0);
   P_.fill(0.0);

   n_x_ = x_.size();

   n_aug_ = n_x_ + 2;

   lambda_ = 3 - n_aug_;
   Xsig_aug_ = Eigen::MatrixXd(UKF::n_aug_, 2 * UKF::n_aug_ + 1);
   Xsig_aug_.fill(0.0);

   Xsig_pred_ = Eigen::MatrixXd(UKF::n_x_, 2 * UKF::n_aug_ + 1);
   Xsig_pred_.fill(0.0);

   // create vector for weights
   weights_ = Eigen::VectorXd(2 * UKF::n_aug_ + 1);
   weights_.fill(0.0);


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    if (!is_initialized_)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double range = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double px = range * cos(phi);
            double py = range * sin(phi);

            x_ <<  range * cos(phi), range * sin(phi), 0, 0, 0;

            P_(0,0) = std_radr_*std_radr_*cos(phi)*cos(phi) + range*range*sin(phi)*sin(phi)*std_radphi_*std_radphi_;
            P_(1,1) = std_radr_*std_radr_*sin(phi)*sin(phi) + range*range*cos(phi)*cos(phi)*std_radphi_*std_radphi_;
            P_(2,2) = 1.0;
            P_(3,3) = 1.0; //(2*M_PI)*(2*M_PI);
            P_(4,4) = 1.0;
            RadarM_ << 0.0, 0.0, 0.0;
        }

        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

            P_(0,0) = std_laspx_*std_laspx_;
            P_(1,1) = std_laspy_*std_laspy_;
            P_(2,2) = 1.0;
            P_(3,3) = 1.0; //(2.0*M_PI)*(2.0*M_PI);
            P_(4,4) = 1.0;
            LidarM_ << 0.0, 0.0;
        }

        double w1 = UKF::lambda_ / (UKF::lambda_ + UKF::n_aug_);
        double wn = 0.5 / (UKF::lambda_ + UKF::n_aug_);

          for(int i = 0; i < weights_.size(); ++i)
          {
              if (i == 0)
              {
                  weights_(i) = w1;
              } else {
                  weights_(i) = wn;
              }
          }


        time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
    }

    double delta_t = (meas_package.timestamp_ - time_us_) / 1.0e+6;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	UKF::Prediction(delta_t);


   if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {

   RadarM_(0) = meas_package.raw_measurements_(0);
   RadarM_(1) = meas_package.raw_measurements_(1);
   RadarM_(2) = meas_package.raw_measurements_(2);

   UKF::UpdateRadar();

   }

   if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {

    LidarM_(0) = meas_package.raw_measurements_(0);
    LidarM_(1) = meas_package.raw_measurements_(1);

    UKF::UpdateLidar();
   }


}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */


    Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(UKF::n_aug_, 2 * UKF::n_aug_ + 1);
    Xsig_aug.fill(0.0);

    AugmentedSigmaPoints(&Xsig_aug);

    UKF::Xsig_aug_.block(0, 0, Xsig_aug.rows(), Xsig_aug.cols()) = Xsig_aug.block(0, 0, Xsig_aug.rows(), Xsig_aug.cols());


    Eigen::MatrixXd Xsig_pred = Eigen::MatrixXd(UKF::n_x_, 2 * UKF::n_aug_ + 1);
    Xsig_pred.fill(0.0);

    SigmaPointPrediction(&Xsig_pred, delta_t);

    UKF::Xsig_pred_.block(0,0, Xsig_pred.rows(),Xsig_pred.cols()) = Xsig_pred.block(0,0, Xsig_pred.rows(),Xsig_pred.cols());

    Eigen::VectorXd x = Eigen::VectorXd(UKF::n_x_);
    Eigen::MatrixXd P = Eigen::MatrixXd(UKF::n_x_, UKF::n_x_);
    P.fill(0.0);
    x.fill(0.0);


    PredictMeanAndCovariance(&x, &P);

    UKF::x_ = x.head(UKF::n_x_);
    UKF::P_.block(0, 0, P_.rows(), P_.cols()) = P.block(0, 0, P.rows(), P.cols());

 }

void UKF::UpdateLidar(void) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

   Eigen::VectorXd z_pred = Eigen::VectorXd(2);
   z_pred.fill(0.0);
   Eigen::MatrixXd S = Eigen::MatrixXd(2, 2);
   S.fill(0.0);
   Eigen::MatrixXd Zsig = Eigen::MatrixXd(2, 2 * UKF::n_aug_ + 1);
   Zsig.fill(0.0);

   PredictLidarMeasurement(&z_pred, &S, &Zsig);
   UpdateLidarState(&z_pred, &S, &Zsig);
}

void UKF::UpdateRadar(void) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   Eigen::VectorXd z_pred = Eigen::VectorXd(3);
   Eigen::MatrixXd S = Eigen::MatrixXd(3, 3);
   Eigen::MatrixXd Zsig = Eigen::MatrixXd(2, 2 * UKF::n_aug_ + 1);;

   PredictRadarMeasurement(&z_pred, &S, &Zsig);

   UpdateRadarState(&z_pred, &S, &Zsig);

}

void UKF::GenerateSigmaPoints(Eigen::MatrixXd* Xsig_out) {

      // create sigma point matrix
      Eigen::MatrixXd Xsig = Eigen::MatrixXd(UKF::n_x_, 2 * UKF::n_x_ + 1);

      // calculate square root of P
      Eigen::MatrixXd A = UKF::P_.llt().matrixL();

      int col_pos = 0;
      Xsig.col(col_pos) = UKF::x_;

      col_pos++;

      Eigen::MatrixXd Psq = A  * sqrt(UKF::lambda_ + UKF::n_x_);


      for (int i = 0; i < Psq.cols(); ++i)
        {
            VectorXd tmp = Psq.col(i) + UKF::x_;
            Xsig.col(col_pos) = tmp;
            col_pos++;
        }

      for (int i = 0; i < Psq.cols(); ++i)
        {
            VectorXd tmp = UKF::x_ - Psq.col(i);
            Xsig.col(col_pos) = tmp;
            col_pos++;
        }

        *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(Eigen::MatrixXd* Xsig_out) {

      // create augmented mean vector
      Eigen::VectorXd x_aug = Eigen::VectorXd(UKF::n_aug_);
      x_aug.fill(0.0);

      // create augmented state covariance
      Eigen::MatrixXd P_aug = Eigen::MatrixXd(UKF::n_aug_, UKF::n_aug_);

      // create sigma point matrix
      Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(UKF::n_aug_, 2 * UKF::n_aug_ + 1);
      Xsig_aug.fill(0.0);
      //cout << "AA" << endl;

      // create augmented mean state
        x_aug.head(UKF::n_x_) = UKF::x_;
        x_aug(UKF::n_aug_ - 2) = 0.0;
        x_aug(UKF::n_aug_ - 1) = 0.0;

      // create augmented covariance matrix

        P_aug.fill(0.0);
        P_aug.topLeftCorner(UKF::n_x_, UKF::n_x_)  = UKF::P_.block(0,0,n_x_, n_x_);
        //std::cout << "P_aug = " << std::endl << P_aug << std::endl;

        Eigen::MatrixXd Q_noise(2, 2);
        Q_noise << UKF::std_a_ * UKF::std_a_, 0,
                    0   , UKF::std_yawdd_ * UKF::std_yawdd_;
        P_aug.bottomRightCorner(2, 2) = Q_noise;


      // create square root matrix
        Eigen::MatrixXd L = P_aug.llt().matrixL();

      int col_pos = 0;
      Xsig_aug.col(col_pos) = x_aug;

      col_pos++;

      Eigen::MatrixXd Psq = L  * sqrt(UKF::lambda_ + UKF::n_aug_);

      for (int i = 0; i < Psq.cols(); ++i)
        {
            VectorXd tmp = Psq.col(i) + x_aug;
            Xsig_aug.col(col_pos) = tmp;
            col_pos++;
        }

      for (int i = 0; i < Psq.cols(); ++i)
        {
            VectorXd tmp = x_aug - Psq.col(i);
            Xsig_aug.col(col_pos) = tmp;
            col_pos++;
        }

      *Xsig_out = Xsig_aug;
}



void UKF::SigmaPointPrediction(Eigen::MatrixXd* Xsig_out, double delta_t) {

    Eigen::MatrixXd Xsig_pred = Eigen::MatrixXd(UKF::n_x_, 2 * UKF::n_aug_ + 1);


  // predict sigma points
  for (int i = 0; i < UKF::Xsig_aug_.cols(); ++i)
  {
      Eigen::VectorXd sig_point = UKF::Xsig_aug_.col(i);
      Eigen::VectorXd xk(UKF::n_x_);
      Eigen::VectorXd sig_pred(UKF::n_x_);
      Eigen::VectorXd sig_new(UKF::n_x_);
      Eigen::VectorXd sig_noise_pred(UKF::n_x_);
      xk = sig_point.head(UKF::n_x_);


      double px = sig_point(0);
      double py = sig_point(1);
      double v = sig_point(2);
      double psi = sig_point(3);
      double psi_d = sig_point(4);
      double nu_acc = sig_point(5);
      double nu_psi_dd = sig_point(6);

      sig_noise_pred << 0.5 * (delta_t * delta_t) * cos(psi) * nu_acc,
                        0.5 * (delta_t * delta_t) * sin(psi) * nu_acc,
                        delta_t * nu_acc,
                        0.5 * (delta_t * delta_t) * nu_psi_dd,
                        delta_t * nu_psi_dd;
        // avoid division by zero
      if (fabs(psi_d) > 0.001) // psi_d is not equal to 0
      {
          sig_pred << (v/psi_d) * (sin(psi + psi_d * delta_t) - sin(psi)),
                     (v/psi_d) * (-cos(psi + psi_d * delta_t) + cos(psi)),
                     0.0,
                     psi_d * delta_t,
                     0.0;
      } else {
                    sig_pred << v * cos(psi) * delta_t,
                                v * sin(psi) * delta_t,
                                0.0,
                                psi_d * delta_t,
                                0.0;
      }

      sig_new = xk + sig_pred + sig_noise_pred;
      // write predicted sigma points into right column
      Xsig_pred.col(i) = sig_new;

  }

  // write result
  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(Eigen::VectorXd* x_out, Eigen::MatrixXd* P_out) {

  // create vector for predicted state
  Eigen::VectorXd x = Eigen::VectorXd(UKF::n_x_);


  // create covariance matrix for prediction
  Eigen::MatrixXd P = Eigen::MatrixXd(UKF::n_x_, UKF::n_x_);



  // predict state mean

  x.fill(0.0);
  P.fill(0.0);
  for (int i = 0; i < UKF::Xsig_pred_.cols(); ++i)
  {
      x += weights_(i) * UKF::Xsig_pred_.col(i);
  }

  // predict state covariance matrix

  Eigen::VectorXd tmp;

  for (int i = 0; i < UKF::Xsig_pred_.cols(); ++i)
  {
      tmp = UKF::Xsig_pred_.col(i) - x;
      P += weights_(i) * (tmp * tmp.transpose());
  }

  *x_out = x;
  *P_out = P;
}

void UKF::PredictRadarMeasurement(Eigen::VectorXd* z_out, Eigen::MatrixXd* S_out, Eigen::MatrixXd* Zsig_out) {

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * UKF::n_aug_ + 1);

  // mean predicted measurement
  Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < UKF::Xsig_pred_.cols(); ++i)
  {
      VectorXd sig_p = UKF::Xsig_pred_.col(i);
      double px = sig_p(0);
      double py = sig_p(1);
      double v = sig_p(2);
      double psi = sig_p(3);

      double ro = sqrt(px * px + py * py);
      double phi = atan2(py , px);
      double ro_d = (px * cos(psi) * v + py * sin(psi) * v) / ro;

      Zsig(0,i) = ro;
      Zsig(1,i) = phi;
      Zsig(2,i) = ro_d;


  }
  // calculate mean predicted measurement

  for (int i = 0; i < Zsig.cols(); ++i)
  {
      z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S

  Eigen::MatrixXd R(n_z, n_z);
  R << UKF::std_radr_ * UKF::std_radr_, 0, 0,
        0, UKF::std_radphi_ * UKF::std_radphi_, 0,
        0, 0, UKF::std_radrd_ * UKF::std_radrd_;

  for (int i = 0; i < Zsig.cols(); ++i)
  {
      MatrixXd A = Zsig.col(i) - z_pred;
      while (A(1)> M_PI) A(1)-=2.*M_PI;
      while (A(1)<-M_PI) A(1)+=2.*M_PI;

      S += weights_(i) * (A * A.transpose());
  }

  S += R;

    // write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;

}

void UKF::PredictLidarMeasurement(Eigen::VectorXd* z_out, Eigen::MatrixXd* S_out, Eigen::MatrixXd* Zsig_out) {

  // set measurement dimension, lida can measure x and y
  int n_z = 2;

  // create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * UKF::n_aug_ + 1);
  Zsig.fill(0.0);

  // mean predicted measurement
  Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  //cout << "HHH" << endl;
  for (int i = 0; i < UKF::Xsig_pred_.cols(); ++i)
  {
      VectorXd sig_p = UKF::Xsig_pred_.col(i);
      double px = sig_p(0);
      double py = sig_p(1);
      double v = sig_p(2);
      double psi = sig_p(3);

      VectorXd meas_vec(2);
      meas_vec << px, py;

      Zsig.col(i) = meas_vec;

  }
  // calculate mean predicted measurement


  for (int i = 0; i < Zsig.cols(); ++i)
  {
       z_pred+= weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S

  Eigen::MatrixXd R(n_z, n_z);
  R << UKF::std_laspx_ * UKF::std_laspx_, 0,
        0, UKF::std_laspy_ * UKF::std_laspy_;

  for (int i = 0; i < Zsig.cols(); ++i)
  {
      MatrixXd A = Zsig.col(i) - z_pred;

      S += weights_(i) * (A * A.transpose());
  }

  S += R;
  // write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::UpdateRadarState(Eigen::VectorXd* z_pred, Eigen::MatrixXd* S, Eigen::MatrixXd* Zsig) {

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for cross correlation Tc
  Eigen::MatrixXd Tc = Eigen::MatrixXd(UKF::n_x_, n_z);

  Tc.fill(0.0);
  // calculate cross correlation matrix


  for (int i = 0; i < (*Zsig).cols(); ++i)
  {
    VectorXd Ax = UKF::Xsig_pred_.col(i) - UKF::x_;

    while (Ax(3)> M_PI) Ax(3)-=2.*M_PI;
    while (Ax(3)<-M_PI) Ax(3)+=2.*M_PI;

    VectorXd Az = (*Zsig).col(i) - (*z_pred);

    while (Az(1)> M_PI) Az(1)-=2.*M_PI;
    while (Az(1)<-M_PI) Az(1)+=2.*M_PI;

    Tc += weights_(i) * (Ax * Az.transpose());
  }

  // calculate Kalman gain K;

  Eigen::MatrixXd K = Tc * (*S).inverse();

  // update state mean and covariance matrix

  Eigen::VectorXd diff_z = RadarM_ - (*z_pred);

  while (diff_z(1)> M_PI) diff_z(1)-=2.*M_PI;
  while (diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;

  UKF::x_ += K * diff_z;

  UKF::P_ = UKF::P_ - (K * ((*S) * K.transpose()));

}


void UKF::UpdateLidarState(Eigen::VectorXd* z_pred, Eigen::MatrixXd* S, Eigen::MatrixXd* Zsig) {

  // set measurement dimension, Lidar can measure position x and y
  int n_z = 2;

  // create matrix for cross correlation Tc
  Eigen::MatrixXd Tc = Eigen::MatrixXd(UKF::n_x_, n_z);

  Tc.fill(0.0);
  // calculate cross correlation matrix
  for (int i = 0; i < Zsig->cols(); ++i)
  {
    Eigen::VectorXd Ax = UKF::Xsig_pred_.col(i) - UKF::x_;
    Eigen::VectorXd Az = (*Zsig).col(i) - (*z_pred);

    Tc += weights_(i) * (Ax * Az.transpose());
  }

  // calculate Kalman gain K;

  Eigen::MatrixXd K = Tc * S->inverse();


  UKF::x_ = UKF::x_ + K * (LidarM_ - (*z_pred));

  UKF::P_ = UKF::P_ - (K * ((*S) * K.transpose()));

}
