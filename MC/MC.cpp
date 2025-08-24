#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main()
{
    long long n_particles = 100000;
    double radius = 10.0;
	int n_bins = 100;
    double bin_width = radius / static_cast<double>(n_bins);
    double SigT = 1.0;

    mt19937 gen(1);
    uniform_real_distribution<double> dis(0.0, 1.0);

	VectorXd n_col = VectorXd::Zero(n_bins); // average
	VectorXd n_col2 = VectorXd::Zero(n_bins); // variance

    for (long long i = 0; i < n_particles; i++) {
        double x0 = 0.0, y0 = 0.0, z0 = 0.0;
        double u = 0.0, v = 0.0, w = 1.0;

        while (true)
        {
            double dist = -log(dis(gen)) / SigT;

            // ... calculate position after flight.
            double x1 = x0 + dist * u;
            double y1 = y0 + dist * v;
            double z1 = z0 + dist * w;


			double airl_dist = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

            // ... identify bin id for current position.
            int ibin = static_cast<int>(airl_dist / bin_width);


            if (ibin >= n_bins) {
                break;  // ... leaked.
            }
            n_col(ibin) += 1.0;   // ... collided.
            n_col2(ibin) += 1.0;  // Note: Original Python code has the same increment as n_col

            // ... determine flight direction (isotropic scattering).
            double costh = 2.0 * dis(gen) - 1.0;
            double phi = 2.0 * M_PI * dis(gen);
            double sinth = sqrt(1.0 - costh * costh);
            u = sinth * cos(phi);
            v = sinth * sin(phi);
            w = costh;

            // ... update position.
            x0 = x1; y0 = y1; z0 = z1;
        }
    }

    // ... calculate bin volumes.
    VectorXd vol(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        double r_outer = bin_width * static_cast<double>(i + 1);
        double r_inner = bin_width * static_cast<double>(i);
        vol(i) = (4.0 * M_PI / 3.0) * (pow(r_outer, 3) - pow(r_inner, 3));
    }

    VectorXd ave = n_col / static_cast<double>(n_particles);
    VectorXd ave2 = n_col2 / static_cast<double>(n_particles);
    VectorXd var = ave2.array() - ave.array().square();

    Eigen::VectorXd flux_mc = ave.array() / SigT / vol.array();
    Eigen::VectorXd flux_mc_stdev = (var / static_cast<double>(n_particles)).array().sqrt() / SigT / vol.array();

    // ... save results to a file for plotting.
    // 変更: ファイル名を "mc_results.txt" に変更
    std::ofstream mc_results("mc_results.txt");
    mc_results << "# x  flux_mc  flux_mc_stdev\n";
    for (int i = 0; i < n_bins; ++i) {
        double x = bin_width * (static_cast<double>(i) + 0.5);
        mc_results << x << "  " << flux_mc(i) << "  " << flux_mc_stdev(i) << "\n";
    }
    mc_results.close();

    // ... save analytical solution to a file.
    // 変更: ファイル名を "analytical_results.txt" に変更
    std::ofstream analytical_results("analytical_results.txt");
    analytical_results << "# x_dif  y_dif\n";
    for (int i = 0; i < 200; ++i) {
        double x_dif = (radius / 200.0) * (static_cast<double>(i) + 0.5);
        if (x_dif > 0) {
            double y_dif = 3.0 * SigT / (4.0 * M_PI) * (1.0 / x_dif - 1.0 / radius);
            analytical_results << x_dif << "  " << y_dif << "\n";
        }
    }
    analytical_results.close();

    std::cout << "Simulation finished. Results saved to '.txt' files." << std::endl;

    // 追加: Pythonスクリプトを呼び出して可視化
    std::cout << "Now launching visualization script..." << std::endl;
    int status = system("python Visualize.py");
    if (status != 0) {
        std::cerr << "Pythonスクリプトの実行に失敗しました。" << std::endl;
    }


}
