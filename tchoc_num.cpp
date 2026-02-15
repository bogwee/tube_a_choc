#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

// Paramètres de Sod
const double g = 1.4;

// Zone Gauche
double p_L = 10.0/g;
double rho_L = 8.0;
double a_L = sqrt(g * p_L / rho_L);
double U_L = 0.0;

// Zone Droite
double p_R = 1.0/g;
double rho_R = 1.0;
double a_R = sqrt(g * p_R / rho_R);
double U_R = 0.0;

double calc_dt(vector<vector<double>> w, double dx, double cfl) {
    double max_speed = 0.0;
    int N = w[0].size();
    vector<double> rho = w[0];
    vector<double> rho_U = w[1];
    vector<double> E = w[2];
    vector<double> U(N), p(N), a(N);
    
    for(int i=0; i<N; i++) {
        U[i] = rho_U[i] / rho[i];
        p[i] = (g - 1) * (E[i] - rho_U[i] * U[i] / 2.0);
        a[i] = sqrt(g * p[i] / rho[i]);
    }

    for(int i=0; i<N; i++) {
        double speed = fabs(U[i]) + a[i];
        if(speed > max_speed) {
            max_speed = speed;
        }
    }
    return cfl * dx / max_speed;
}

vector<double> flux_centre(vector<double> w) {
    vector<double> f(w.size());
    double rho = w[0];
    double rho_U = w[1];
    double E = w[2];
    double U = rho_U / rho;
    double p = (g - 1) * (E - rho_U * U / 2.0);
    f[0] = rho_U;
    f[1] = rho * U * U + p;
    f[2] = (E + p) * U;
    return f;
}

vector<vector<double>> lax_wendroff(vector<vector<double>> W, double T, double dx, double cfl=0.95) {
    int N = W[0].size();
    double t = 0.0;
    vector<vector<double>> W_tilde(3, vector<double>(N-1));
    vector<vector<double>> W_new = W;
    double dt = calc_dt(W, dx, cfl);
    while(t < T) {
        double dt = calc_dt(W, dx, cfl);
        if (t + dt > T) dt = T - t;
        for(int i=0; i<N-1; i++) {
            auto f_right = flux_centre({W[0][i+1], W[1][i+1], W[2][i+1]});
            auto f_left  = flux_centre({W[0][i], W[1][i], W[2][i]});
            for(int j=0; j<3; j++) {
                W_tilde[j][i] = (W[j][i]+W[j][i+1])/2.0 - dt/(2*dx) * (f_right[j] - f_left[j]);
            }
        }
        for(int i=1; i<N-1; i++) {
            auto f_tilde_right = flux_centre({W_tilde[0][i], W_tilde[1][i], W_tilde[2][i]});
            auto f_tilde_left  = flux_centre({W_tilde[0][i-1], W_tilde[1][i-1], W_tilde[2][i-1]});
            for(int j=0; j<3; j++) {
                W_new[j][i] = W[j][i] - dt/dx * (f_tilde_right[j] - f_tilde_left[j]);
            }
        }
        W = W_new;
        t += dt;
    }
    return W;
}

vector<vector<double>> maccormack(vector<vector<double>> W, double T, double dx, double cfl=0.95) {
    int N = W[0].size();
    double t = 0.0;
    vector<vector<double>> W_predict(3, vector<double>(N-1));
    vector<vector<double>> W_correct = W;
    while(t < T) {
        double dt = calc_dt(W, dx, cfl);
        if (t + dt > T) dt = T - t;

        for(int i=0; i<N-1; i++) {
            auto f_right = flux_centre({W[0][i+1], W[1][i+1], W[2][i+1]});
            auto f_left  = flux_centre({W[0][i], W[1][i], W[2][i]});
            for(int j=0; j<3; j++) {
                W_predict[j][i] = W[j][i] - dt/dx * (f_right[j] - f_left[j]);
            }
        }

        for(int i=1; i<N-1; i++) {
            auto f_predict_right = flux_centre({W_predict[0][i], W_predict[1][i], W_predict[2][i]});
            auto f_predict_left  = flux_centre({W_predict[0][i-1], W_predict[1][i-1], W_predict[2][i-1]});
            for(int j=0; j<3; j++) {
                W_correct[j][i] = 0.5 * (W[j][i] + W_predict[j][i]) - 0.5 * dt/dx * (f_predict_right[j] - f_predict_left[j]);
            }
        }
        W = W_correct;
        t += dt;
    }
    return W;
}

int main() {

    // Paramètres ex 10.1
    int n = 101;
    double x0 = 0.5;
    double t = 0.2;

    vector<vector<double>> W(3, vector<double>(n));
    for(int i=0; i<n; i++) {
        W[0][i] = (i < n/2) ? rho_L : rho_R;
        W[1][i] = (i < n/2) ? W[0][i] * U_L : W[0][i] * U_R;
        W[2][i] = (i < n/2) ? (p_L/(g-1) + 0.5 * W[1][i] * U_L) : (p_R/(g-1) + 0.5 * W[1][i] * U_R);
    }

    vector<vector<double>> result = maccormack(W, t, 1.0/(n-1));

    for(int i=0; i<n; i++) {
        result[1][i] /= result[0][i]; // Convertir rho_U en U
        result[2][i] = (g - 1) * (result[2][i] - 0.5 * result[1][i] * result[1][i] / result[0][i]); // Convertir E en p
    }


    ofstream file("output_num.csv");
    file << "x,rho,U,p\n";
    for(int i = 0; i < n; ++i) {
        file << i * (1.0 / (n - 1)) << ',' << result[0][i] << ',' << result[1][i] << ',' << result[2][i] << '\n';
    }
    
    file.close();
    return 0;
}