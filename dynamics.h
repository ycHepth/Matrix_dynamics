//
// Created by yc on 2021/6/24.
//

#ifndef MATRIX_DYNAMICS_DYNAMICS_H
#define MATRIX_DYNAMICS_DYNAMICS_H

#include <cmath>
#include <Eigen/Dense>

class dynamics {
private:
    double L1zz, L2zz, L3zz;
    double lx1, lx2, lx3;
    double Im1, Im2, Im3;
    double m1, m2, m3;
    double L1, L2;
    double m1_r, m2_r, m3_r;
    double N1, N2, N3;
    double grav_acc;
    double fv1, fc1, fv2, fc2, fv3, fc3;
    double fv1_m, fc1_m, fv2_m, fc2_m, fv3_m, fc3_m;
public:
    dynamics();

    Eigen::Matrix<double, 3, 3> Inertia_term(Eigen::Vector3d q);

    Eigen::Vector3d Coriolis_term(Eigen::Vector3d dq, Eigen::Vector3d q);

    Eigen::Vector3d Gravity_term(Eigen::Vector3d q);

    Eigen::Vector3d Friction_link(Eigen::Vector3d dq);

    Eigen::Vector3d Friction_motor(Eigen::Vector3d dtheta);

    Eigen::Vector3d
    coupling_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d dq, Eigen::Vector3d q, Eigen::Vector3d ddtheta,
                      Eigen::Vector3d dtheta);

    Eigen::Vector3d link_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d dq, Eigen::Vector3d q, Eigen::Vector3d ddtheta);

    Eigen::Vector3d motor_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d ddtheta, Eigen::Vector3d dtheta);
};

dynamics::dynamics() : L1zz(0.782), L2zz(0.287), L3zz(0.0033), lx1(2.0462), lx2(0.8661), lx3(0.0041) {
    m1 = 4.52;
    m2 = 1.95;
    m3 = 0.95;
    m1_r = 0.37;
    m2_r = 0.37;
    m3_r = 0.23;
    N1 = 120;
    N2 = 80;
    N3 = 80;
    grav_acc = 9.81;
    fv1 = 1.3045;
    fc1 = -0.5821;
    fv2 = -0.7169;
    fc2 = -0.5728;
    fv3 = -0.4565;
    fc3 = 2.3638;
    fv1_m = 19.0233;
    fc1_m = 10.7662;
    fv2_m = 7.2910;
    fc2_m = 5.1601;
    fv3_m = 3.7434;
    fc3_m = 2.2495;
    L1 = 0.4;
    L2 = 0.395;
    Im1 = 0.8786;
    Im2 = 0.4845;
    Im3 = 0.015;
}

Eigen::Matrix<double, 3, 3> dynamics::Inertia_term(Eigen::Vector3d q) {
    double q1 = q(0);
    double q2 = q(1);
    double q3 = q(2);

    double M11 = L1zz + L2zz + L3zz + lx3 * (2 * (L1 * cos(q2 + q3) + L2 * cos(q3))) + lx2 * (2 * L1 * cos(q2));
    double M12 = L2zz + L3zz + lx3 * (2 * L2 * cos(q3) + L1 * cos(q2 + q3)) + lx2 * L1 * cos(q2);
    double M13 = L3zz + lx3 * (L1 * cos(q2 + q3) + L2 * cos(q3));
    double M22 = L2zz + L3zz + lx3 * 2 * L2 * cos(q3);
    double M23 = L3zz + lx3 * L2 * cos(q3);
    double M33 = L3zz;

    Eigen::Matrix<double, 3, 3> Inertia;
    Inertia << M11, M12, M13,
            M12, M22, M23,
            M13, M23, M33;

    return Inertia;

}

Eigen::Vector3d dynamics::Coriolis_term(Eigen::Vector3d dq, Eigen::Vector3d q) {
    double dq1 = dq(0);
    double dq2 = dq(1);
    double dq3 = dq(2);
    double q1 = q(0);
    double q2 = q(1);
    double q3 = q(2);

    double C1, C2, C3;
    C1 = -lx3 *
         (L1 * sin(q2 + q3) * (pow(dq2, 2) + pow(dq3, 2)) + L2 * sin(q3) * pow(dq3, 2) + 2 * L2 * sin(q3) * dq1 * dq3 +
          2 * L2 * sin(q3) * dq2 * dq3 + 2 * L1 * sin(q2 + q3) * dq1 * dq2 +
          2 * L1 * sin(q2 + q3) * dq1 * dq3 + 2 * L1 * sin(q2 + q3) * dq2 * dq3) -
         lx2 * (L1 * sin(q1) * pow(dq2, 2) + 2 * L1 * sin(q2) * dq1 * dq2);

    C2 = lx3 * (L1 * sin(q1 + q2 + q3) * pow(dq1, 2) - L2 * sin(q3) * pow(dq3, 2) - 2 * L2 * sin(q3) * dq1 * dq3 -
                2 * L2 * sin(q3) * dq2 * dq3) + lx2 * L1 * sin(q2) * pow(dq1, 2);

    C3 = lx3 * (L1 * sin(q2 + q3) * pow(dq1, 2) + L2 * sin(q3) * pow(dq1, 2) + L2 * sin(q3) * pow(dq2, 2) +
                2 * L2 * sin(q3) * dq1 * dq2);

    Eigen::Vector3d Coriolis;
    Coriolis << C1, C2, C3;

    return Coriolis;
}

Eigen::Vector3d dynamics::Gravity_term(Eigen::Vector3d q) {

    double q1 = q(0);
    double q2 = q(1);
    double q3 = q(2);

    double G1, G2, G3;

    G1 = lx3 * sin(q1 + q2 + q3) * grav_acc + lx2 * sin(q1 + q2) * grav_acc + lx1 * sin(q1) * grav_acc;
    G2 = lx3 * sin(q1 + q2 + q3) * grav_acc + lx2 * sin(q1 + q2) * grav_acc;
    G3 = lx3 * sin(q1 + q2 + q3) * grav_acc;

    Eigen::Vector3d Gravity_term;
    Gravity_term << G1, G2, G3;
    return Gravity_term;
}


Eigen::Vector3d dynamics::Friction_link(Eigen::Vector3d dq) {
    double dq1 = dq(0);
    double dq2 = dq(1);
    double dq3 = dq(2);

    double friction_l_1, friction_l_2, friction_l_3;

    friction_l_1 = fv1 * dq1 + fc1 * tanh(2 * dq1 / 0.02);
    friction_l_2 = fv2 * dq2 + fc2 * tanh(2 * dq2 / 0.02);
    friction_l_3 = fv3 * dq3 + fc3 * tanh(2 * dq3 / 0.02);

    Eigen::Vector3d Friction_L;
    Friction_L << friction_l_1, friction_l_2, friction_l_3;
    return Friction_L;
}

Eigen::Vector3d dynamics::Friction_motor(Eigen::Vector3d dtheta) {
    double dtheta1 = dtheta(0);
    double dtheta2 = dtheta(1);
    double dtheta3 = dtheta(2);

    double friction_m_1, friction_m_2, friction_m_3;

    friction_m_1 = fv1_m * dtheta1 + fc1_m * tanh(2 * dtheta1 / 0.02);
    friction_m_2 = fv2_m * dtheta2 + fc1_m * tanh(2 * dtheta2 / 0.02);
    friction_m_3 = fv3_m * dtheta3 + fc1_m * tanh(2 * dtheta3 / 0.02);

    Eigen::Vector3d Friction_M;
    Friction_M << friction_m_1, friction_m_2, friction_m_3;
    return Friction_M;
}


Eigen::Vector3d
dynamics::coupling_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d dq, Eigen::Vector3d q, Eigen::Vector3d ddtheta,
                            Eigen::Vector3d dtheta) {

    double q1 = q(0);
    double q2 = q(1);
    double q3 = q(2);

    double dq1 = dq(0);
    double dq2 = dq(1);
    double dq3 = dq(2);

    double ddq1 = ddq(0);
    double ddq2 = ddq(1);
    double ddq3 = ddq(2);

    double dtheta1 = dtheta(0);
    double dtheta2 = dtheta(1);
    double dtheta3 = dtheta(2);

    double ddtheta1 = ddtheta(0);
    double ddtheta2 = ddtheta(1);
    double ddtheta3 = ddtheta(2);

    Eigen::Matrix<double, 3, 3> Inertia = dynamics::Inertia_term(q);
    Eigen::Vector3d Coriolis = dynamics::Coriolis_term(dq, q);
    Eigen::Vector3d Gravity = dynamics::Gravity_term(q);
    Eigen::Vector3d friction_l = dynamics::Friction_link(dq);
    Eigen::Vector3d friction_m = dynamics::Friction_motor(dtheta);

    Eigen::Matrix<double, 3, 3> Inertia_m;
    double Mr11, Mr12, Mr22;
    Mr11 = (m2_r + m3_r) * pow(L1, 2) + m3_r * pow(L2, 2) + 2 * m3_r * L1 * L2 * cos(q2);
    Mr12 = m3_r * pow(L2, 2) + m3_r * L1 * L2 * cos(q2);
    Mr22 = m3_r * pow(L2, 2);

    Inertia_m << Mr11, Mr12, 0,
            Mr12, Mr22, 0,
            0, 0, 0;

    Eigen::Matrix3d B;
    Eigen::Vector3d I_mot;
    I_mot << Im1, Im2, Im3;
    B = I_mot.asDiagonal();


    Eigen::Matrix3d S;
    S << 0, Im2 / N2, Im3 / N3,
            0, 0, Im3 / N3,
            0, 0, 0;

    Eigen::Matrix3d SBS;
    SBS << Im2 / pow(N2, 2) + Im3 / pow(N3, 2), Im3 / pow(N3, 2), 0,
            Im3 / pow(N3, 2), Im3 / pow(N3, 2), 0,
            0, 0, 0;

    Eigen::Matrix3d H;
    H = Inertia + SBS;

    Eigen::Vector3d tau_m;
    tau_m = (H + S.transpose()) * ddq + (S + B) * ddtheta + Coriolis + Gravity + friction_l + friction_m;

    return tau_m;
}


Eigen::Vector3d
dynamics::link_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d dq, Eigen::Vector3d q, Eigen::Vector3d ddtheta) {
    Eigen::Matrix<double, 3, 3> Inertia = dynamics::Inertia_term(q);
    Eigen::Vector3d Coriolis = dynamics::Coriolis_term(dq, q);
    Eigen::Vector3d Gravity = dynamics::Gravity_term(q);
    Eigen::Vector3d friction_l = dynamics::Friction_link(dq);

    Eigen::Matrix3d S;
    S << 0, Im2 / N2, Im3 / N3,
            0, 0, Im3 / N3,
            0, 0, 0;

    Eigen::Matrix3d SBS;
    SBS << Im2 / pow(N2, 2) + Im3 / pow(N3, 2), Im3 / pow(N3, 2), 0,
            Im3 / pow(N3, 2), Im3 / pow(N3, 2), 0,
            0, 0, 0;

    Eigen::Matrix3d H;
    H = Inertia + SBS;

    Eigen::Vector3d dyn_link;

    dyn_link = H * ddq + Coriolis + Gravity + S * ddtheta + friction_l;

    return dyn_link;
}

Eigen::Vector3d dynamics::motor_dynamics(Eigen::Vector3d ddq, Eigen::Vector3d ddtheta, Eigen::Vector3d dtheta) {
    Eigen::Vector3d I_mot;
    I_mot << Im1, Im2, Im3;
    I_mot.asDiagonal();

    /* TODO:
     *      Here need stiffness of spring,
     */
    return I_mot;
}

#endif //MATRIX_DYNAMICS_DYNAMICS_H
