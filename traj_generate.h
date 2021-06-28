//
// Created by yc on 2021/6/24.
//

#ifndef MATRIX_DYNAMICS_TRAJ_GENERATE_H
#define MATRIX_DYNAMICS_TRAJ_GENERATE_H

/**
 * Fourier Series of Joint angle
 */

#define Ts 0.005
#define pi 3.1415926


/*
 *Summary:Generate Fourier series of (q,dq,ddq) in radius(rad,rad/s,rad/s^2)
 *
 *param: a0 :zero-order coefficient
 *param: a  :array of a coefficient
 *param: b  :array of b coefficient
 *param: w  :base angle freq of fourier
 */
double Fourier_series_velocity_10th(int curve_count, double a0, const double *a, const double *b, double w);

double Fourier_series_position_10th(int curve_count, double a0, const double *a, const double *b, double w);

double Fourier_series_acceleration_10th(int curve_count, double a0, const double *a, const double *b, double w);

double Fourier_series_velocity_16th(int curve_count, double a0, const double *a, const double *b, double w);

double Fourier_series_position_16th(int curve_count, double a0, const double *a, const double *b, double w);

double Fourier_series_acceleration_16th(int curve_count, double a0, const double *a, const double *b, double w);

#endif //MATRIX_DYNAMICS_TRAJ_GENERATE_H
