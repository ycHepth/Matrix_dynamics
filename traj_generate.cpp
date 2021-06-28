//
// Created by yc on 2021/6/24.
//
#include <cmath>
#include "traj_generate.h"

double Fourier_series_velocity_10th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<10;i++){
        sum += a[i]*cos((i+1)*curve_count*Ts*w) + b[i]*sin((i+1)*curve_count*Ts*w);
    }
    return sum;
}

double Fourier_series_position_10th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<10;i++){
        sum += a[i]/(w*(i+1))*sin((i+1)*curve_count*Ts*w) - b[i]/(w*(i+1))*cos((i+1)*curve_count*Ts*w);
    }
    return sum + a0;
}

double Fourier_series_acceleration_10th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<10;i++){
        sum += -a[i]*(w*(i+1))*sin((i+1)*curve_count*Ts*w) + b[i]*(w*(i+1))*cos((i+1)*curve_count*Ts*w);
    }
    return sum;
}


double Fourier_series_velocity_16th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<16;i++){
        sum += a[i]*cos((i+1)*curve_count*Ts*w) + b[i]*sin((i+1)*curve_count*Ts*w);
    }
    return sum;
}

double Fourier_series_position_16th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<16;i++){
        sum += a[i]/(w*(i+1))*sin((i+1)*curve_count*Ts*w) - b[i]/(w*(i+1))*cos((i+1)*curve_count*Ts*w);
    }
    return sum + a0;
}

double Fourier_series_acceleration_16th(int curve_count, double a0, const double* a, const double* b , double w){
    double sum = 0;
    for(int i=0; i<16;i++){
        sum += -a[i]*(w*(i+1))*sin((i+1)*curve_count*Ts*w) + b[i]*(w*(i+1))*cos((i+1)*curve_count*Ts*w);
    }
    return sum;
}