#pragma once
#include  "matrix.hpp"
struct TInputData {
        Matrix F; // process evolution
        Matrix H; // observation model
        Matrix Q; // process noise covariance
        Matrix R; // observation noise covariance
        Matrix P; // estimation error covariance
        my_Vector x;
};

class KalmanFilter {
private:
    Matrix F; // process evolution
    Matrix H; // observation model
    Matrix Q; // process noise covariance
    Matrix R; // observation noise covariance
    Matrix P; // estimation error covariance
    my_Vector x;
public:
    KalmanFilter() = default;
    KalmanFilter(TInputData & settings) {
        this->F = settings.F;
        this->H = settings.H;
        this->Q = settings.Q;
        this->R = settings.R;
        this->P = settings.P;
        this->x = settings.x;
    }

    my_Vector update(my_Vector z) {
        Matrix F = this->F; // process evolution
        Matrix H = this->H; // observation model
        Matrix Q = this->Q; // process noise covariance
        Matrix R = this->R; // observation noise covariance
        Matrix P = this->P; // estimation error covariance
        my_Vector x = this->x ;
        // Predict
        my_Vector xPred = F * x;
        //xPred.print();
        Matrix F_T = F.transpose();
        Matrix PPred = F * P * F_T + Q;
        // Update
        // Innovation or measurement residual
        my_Vector y = z + (-1) * (H * xPred);
        //std::cout<<"y = "<<std::endl;
        //y.print();
        // Kalman gain
        Matrix tmp = H * PPred * H.transpose() + R;
        //tmp.print();
        Matrix tmp_inv_diag = tmp.inverseDiagonal();
        //tmp_inv_diag.print();
        Matrix K = PPred * (H.transpose()) * tmp_inv_diag;
        // Posteriori state estimate
        //K.print();
        my_Vector xEst = xPred + (K * y);
        //std::cout<<"XEST = "<<std::endl;
        //xEst.print();
        // Posteriori estimate covariance
        Matrix PEst = PPred - (K * H * PPred);

        this->x = xEst;
        this->P = PEst;

        //std::cout<<"Filter update"<<std::endl;
        //xEst.print();
        //PEst.print();
        return xEst;
    }

    my_Vector get_state() {
        return this->x;
    }

    void set_QMatrix(Matrix Q) {
        this->Q = Q;
    }

    Matrix get_QMatrix() {
        return this->Q;
    }

    void set_RMatrix(Matrix R) {
        this->R = R;
    }

    void set_HMatrix(Matrix H) {
        this->H = H;
    }

    void set_FMatrix(Matrix F) {
        this->F = F;
    }

    Matrix get_FMatrix() {
        return this->F;
    }
};