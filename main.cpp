#include <iostream>
#include "dynamics.h"
#include "matplotlibcpp.h"
#include <vector>
#include "traj_generate.h"


int main() {

    double w = 2 * pi * 0.1;  // angle freq  = 2*pi*f where f is base freq = 0.1Hz

    /**
     * Link Fourier
     */
    double a_hip_16th[16] = {0.110280730317861, -0.0837407650946106, 0.383449444920027, 0.269874826780982,
                             0.0449525398897773, -0.116363758926572, -0.139824218255315, 0.0491741727812546,
                             -0.173146867919171, 0.209648483989100, 0.0795573685101979, 0.0178965334825925,
                             -0.0605647094675307, 0.230815144252535, -0.0589089142969703, 0.0963730503317145};
    double b_hip_16th[16] = {-0.0401159292237624, -0.107400354190721, -0.189507209710420, 0.0770229595356877,
                             -0.204303369639698, -0.209261935644037, 0.160567039655744, -0.0765735219793643,
                             -0.196200624738493, 0.00608487418892085, 0.165902683107064, -0.389438219427507,
                             -0.142950670590139, 0.0270174141926537, -0.109677331856646, -0.0411545432598413};
    double a0_hip_16th = 0.3000;

    double a_knee_16th[16] = {0.0537910264571673, -0.461494615589504, 0.237258931031145, 0.324399144615701,
                              0.476559263140522, -0.731664684587670, -0.140452419488474, 0.107517956978429,
                              0.221960660882815, -0.0343337591132675, -0.0247055475132788, 0.00553950426478707,
                              -0.0146310140544731, -0.0767083350281416, 0.0265420903349684, -0.0811607370308330};
    double b_knee_16th[16] = {0.0996916496787676, 0.0153102539246314, -0.0117252970830793, -0.551009591246249,
                              -0.405611186257982, -0.483511346676119, 0.0501176988875974, -0.401968060355294,
                              0.0216613784158418, 0.0384931977617943, -0.0448590393100210, 0.0713910443929834,
                              0.0780659876943586, -0.0492494006970769, 0.0555463778843669, 0.0333885685962466};
    double a0_knee_16th = -0.9135;

    double a_ankle_16th[16] = {0.230023687945328, -0.296448774521793, -0.128164209740401, 0.299394818720724,
                               -0.0507820226164138, 0.275469544401282, -0.720980602886862, 0.584351783428297,
                               -0.0558452386359274, -0.0987907654225857, 0.0263278093999989, -0.0383638763219277,
                               0.0328498804731919, -0.00515509748961038, 0.0229807632110711, -0.00574724907244697};

    double b_ankle_16th[16] = {-0.0498867763441066, 0.317156463663102, -0.0621841423508625, -0.239395596120155,
                               -1.09075464356417, -1.03845424890843, -0.443929255922571, 0.137379727124741,
                               -0.141476253832814, -0.0202428208627133, 0.0327102687375745, -0.0100264588978977,
                               -0.0296621890704325, 0.0147218785506873, -0.00965295135283498, 0.0365015478474211};
    double a0_ankle_16th = 1.7061;

    /**
     * Motor Fourier
     */
    double a_hip_16th_m[16] = {0.168207788489660, -0.167178990483885, 0.558384089077667, 0.410657095119088,
                               0.0556764386955203, -0.212893277048878, -0.216942270854112, 0.0612663472826918,
                               -0.260166615849425, 0.143884295063219, 0.000937275192964842, 0.000711631502457993,
                               0.000302710817865032, 0.000940620958158975, 0.000557756914811288, 0.000244378225285001};
    double b_hip_16th_m[16] = {-0.0649016060148844, -0.138975155904158, -0.265439959674068, 0.0652071750710400,
                               -0.315537345897107, -0.226327946409821, 0.0964038085938133, 0.0126069475027638,
                               -0.243006865662500, -0.0808287483607801, 0.000590641783628319, -0.000410655277169853,
                               4.35822147570221e-05, 0.000725347956738342, 0.00107053670729631, 0.000453855104377182};
    double a0_hip_16th_m = 0.3107;

    double a_knee_16th_m[16] = {0.0426311177072208, -0.422595104985484, 0.204756864417496, 0.290197185364587,
                                0.461627428410567, -0.713279585958261, -0.127556693844491, 0.0876091302316757,
                                0.200905795501292, -0.0107392057279383, 0.000809489045794428, 6.80502975000904e-05,
                                0.000613508953163489, 0.000689914533482125, 0.000766948971526059, 0.00129439195515987};
    double b_knee_16th_m[16] = {0.0963997338051407, 0.0144710072139928, 0.0112603176260297, -0.535068964600796,
                                -0.363031532952293, -0.475634285004112, 0.0646657552937975, -0.424586824426446,
                                0.0161421587703932, 0.0640723305090677, 0.000372031401135905, -0.00105398968483520,
                                -0.00109758674632488, 0.000362535760027869, 0.000707856094710142, 0.000252241957130555};
    double a0_knee_16th_m = -0.8460;

    double a_ankle_16th_m[16] = {0.237286144879554, -0.305171175811538, -0.131973729101943, 0.313864360278701,
                                 -0.0680603227498872, 0.291704922990977, -0.755471010040142, 0.615778908345386,
                                 -0.0862063609302558, -0.0758758003451952, 0.00135735659545185, -0.000504198590814871,
                                 0.000986119886108717, 0.000276425507884861, 0.000645400453449564, 0.00193910848716422};
    double b_ankle_16th_m[16] = {-0.0573153412156880, 0.334662602659168, -0.0776718233135361, -0.231285044137639,
                                 -1.13835132357674, -1.04234311457895, -0.463064161534706, 0.152375897731873,
                                 -0.144715310778429, -0.0299352644167070, 0.00110129313571922, -0.000215163130284042,
                                 -0.00195804068199299, 0.000336078772866822, -0.000425418108846921,
                                 0.00220582042094104};
    double a0_ankle_16th_m = 2.1049;


    double q1, q2, q3;
    double dq1, dq2, dq3;
    double ddq1, ddq2, ddq3;

    double theta1, theta2, theta3;
    double dtheta1, dtheta2, dtheta3;
    double ddtheta1, ddtheta2, ddtheta3;

    dynamics SEA_dynamics;

    int simulation_size = 2000;

    std::vector<double> q1_save(simulation_size);
    std::vector<double> q2_save(simulation_size);
    std::vector<double> q3_save(simulation_size);
    std::vector<double> theta1_save(simulation_size);
    std::vector<double> theta2_save(simulation_size);
    std::vector<double> theta3_save(simulation_size);

    std::vector<double> tau1(simulation_size);
    std::vector<double> tau2(simulation_size);
    std::vector<double> tau3(simulation_size);

    std::vector<double> t(simulation_size);

    Eigen::IOFormat PrettyPrint(4, 0, ",", "\n", "[", "]", "[", "]");

    for (int i = 0; i < simulation_size; i++) {
        theta1 = Fourier_series_position_16th(i, a0_hip_16th_m, a_hip_16th_m, b_hip_16th_m, w);
        theta2 = Fourier_series_position_16th(i, a0_knee_16th_m, a_knee_16th_m, b_knee_16th_m, w);
        theta3 = Fourier_series_position_16th(i, a0_ankle_16th_m, a_ankle_16th_m, b_ankle_16th_m, w);

        dtheta1 = Fourier_series_velocity_16th(i, a0_hip_16th_m, a_hip_16th_m, b_hip_16th_m, w);
        dtheta2 = Fourier_series_velocity_16th(i, a0_knee_16th_m, a_knee_16th_m, b_knee_16th_m, w);
        dtheta3 = Fourier_series_velocity_16th(i, a0_ankle_16th_m, a_ankle_16th_m, b_ankle_16th_m, w);

        ddtheta1 = Fourier_series_acceleration_16th(i, a0_hip_16th_m, a_hip_16th_m, b_hip_16th_m, w);
        ddtheta2 = Fourier_series_acceleration_16th(i, a0_knee_16th_m, a_knee_16th_m, b_knee_16th_m, w);
        ddtheta3 = Fourier_series_acceleration_16th(i, a0_ankle_16th_m, a_ankle_16th_m, b_ankle_16th_m, w);

        q1 = Fourier_series_position_16th(i, a0_hip_16th, a_hip_16th, b_hip_16th, w);
        q2 = Fourier_series_position_16th(i, a0_knee_16th, a_knee_16th, b_knee_16th, w);
        q3 = Fourier_series_position_16th(i, a0_ankle_16th, a_ankle_16th, b_ankle_16th, w);

        dq1 = Fourier_series_velocity_16th(i, a0_hip_16th, a_hip_16th, b_hip_16th, w);
        dq2 = Fourier_series_velocity_16th(i, a0_knee_16th, a_knee_16th, b_knee_16th, w);
        dq3 = Fourier_series_velocity_16th(i, a0_ankle_16th, a_ankle_16th, b_ankle_16th, w);

        ddq1 = Fourier_series_acceleration_16th(i, a0_hip_16th, a_hip_16th, b_hip_16th, w);
        ddq2 = Fourier_series_acceleration_16th(i, a0_knee_16th, a_knee_16th, b_knee_16th, w);
        ddq3 = Fourier_series_acceleration_16th(i, a0_ankle_16th, a_ankle_16th, b_ankle_16th, w);

        Eigen::Vector3d q_current, dq_current, ddq_current;
        q_current << q1, q2, q3;
        dq_current << dq1, dq2, dq3;
        ddq_current << ddq1, ddq2, ddq3;

        Eigen::Vector3d theta_current, dtheta_current, ddtheta_current;
        theta_current << theta1, theta2, theta3;
        dtheta_current << dtheta1, dtheta2, dtheta3;
        ddtheta_current << ddtheta1, ddtheta2, ddtheta3;

        Eigen::Vector3d tau_m = SEA_dynamics.coupling_dynamics(ddq_current,dq_current,q_current,ddtheta_current,dtheta_current);

        t[i] = i / 200.0; // assign to time
        q1_save[i] = q1;
        q2_save[i] = q2;
        q3_save[i] = q3;
        theta1_save[i] = theta1;
        theta2_save[i] = theta2;
        theta3_save[i] = theta3;

        tau1[i] = tau_m(0);
        tau2[i] = tau_m(1);
        tau3[i] = tau_m(2);
    }

    Eigen::Vector3d I;
    I << 1, 2, 3;
    Eigen::Matrix3d B;
    B = I.asDiagonal();
    std::cout << B.format(PrettyPrint) << std::endl;


    matplotlibcpp::figure(1);
    matplotlibcpp::subplot(3,1,1);
    matplotlibcpp::plot(t, q1_save, "r-"); // must specified the length of input
    matplotlibcpp::plot(t, theta1_save, "b-");
    matplotlibcpp::named_plot("q1", t, q1_save);
    matplotlibcpp::named_plot("theta1", t, theta1_save);
    matplotlibcpp::legend();
    matplotlibcpp::title("trajectory");
    matplotlibcpp::ylabel("position [rad]");

    matplotlibcpp::subplot(3,1,2);
    matplotlibcpp::plot(t, q2_save, "r-"); // must specified the length of input
    matplotlibcpp::plot(t, theta2_save, "b-");
    matplotlibcpp::named_plot("q2", t, q2_save);
    matplotlibcpp::named_plot("theta2", t, theta2_save);
    matplotlibcpp::legend();
    matplotlibcpp::ylabel("position [rad]");

    matplotlibcpp::subplot(3,1,3);
    matplotlibcpp::plot(t, q3_save, "r-"); // must specified the length of input
    matplotlibcpp::plot(t, theta3_save, "b-");
    matplotlibcpp::named_plot("q3", t, q3_save);
    matplotlibcpp::named_plot("theta3", t, theta3_save);
    matplotlibcpp::legend();
    matplotlibcpp::xlabel("time [s]");
    matplotlibcpp::ylabel("position [rad]");


    matplotlibcpp::figure(2);
    matplotlibcpp::subplot(3,1,1);
    matplotlibcpp::plot(t,tau1,"r-");
    matplotlibcpp::xlim(0,10);
    matplotlibcpp::ylabel("torque [Nm]");

    matplotlibcpp::subplot(3,1,2);
    matplotlibcpp::plot(t,tau2,"r-");
    matplotlibcpp::xlim(0,10);
    matplotlibcpp::ylabel("torque [Nm]");

    matplotlibcpp::subplot(3,1,3);
    matplotlibcpp::plot(t,tau3,"r-");
    matplotlibcpp::xlim(0,10);
    matplotlibcpp::xlabel("time [s]");
    matplotlibcpp::ylabel("torque [Nm]");

    matplotlibcpp::show();


    return 0;
}
