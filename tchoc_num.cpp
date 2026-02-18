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

vector<double> dissipation_roe(vector<double> wl, vector<double> wr) {
    double rho_l = wl[0];
    double rho_U_l = wl[1];
    double E_l = wl[2];
    double U_l = rho_U_l / rho_l;
    double p_l = (g - 1) * (E_l - rho_U_l * U_l / 2.0);
    double H_l = (E_l + p_l) / rho_l;
    double a_l = sqrt(g * p_l / rho_l);

    double rho_r = wr[0];
    double rho_U_r = wr[1];
    double E_r = wr[2];
    double U_r = rho_U_r / rho_r;
    double p_r = (g - 1) * (E_r - rho_U_r * U_r / 2.0);
    double H_r = (E_r + p_r) / rho_r;
    double a_r = sqrt(g * p_r / rho_r);

    double d_rho = rho_r - rho_l;
    double d_rho_U = rho_U_r - rho_U_l;
    double d_E = E_r - E_l;

    double R = sqrt(rho_r / rho_l);
    double rho_bar = R * rho_l;
    double U_bar = (R * U_r + U_l) / (R + 1);
    double H_bar = (R * H_r + H_l) / (R + 1);
    double a_bar = sqrt((g - 1) * (H_bar - U_bar * U_bar / 2.0));

    vector<double> l = {fabs(U_bar - a_bar), fabs(U_bar), fabs(U_bar + a_bar)};

    vector<vector<double>> v = {{1, U_bar - a_bar, H_bar - U_bar * a_bar}, 
                                {1, U_bar, 0.5 * U_bar * U_bar}, 
                                {1, U_bar + a_bar, H_bar + U_bar * a_bar}};

    double alpha1 = 0.5 * (g - 1) * (U_bar * U_bar)/(a_bar * a_bar);
    double alpha2 = (g - 1)/(a_bar * a_bar);

    vector<double> alpha = {0.5 * (alpha1 + U_bar/a_bar) * d_rho - 0.5 * (alpha2 * U_bar + 1/a_bar) * d_rho_U + 0.5 * alpha2 * d_E, 
                            (1 - alpha1) * d_rho + alpha2 * U_bar * d_rho_U - alpha2 * d_E, 
                            0.5 * (alpha1 - U_bar/a_bar) * d_rho - 0.5 * (alpha2 * U_bar - 1/a_bar) * d_rho_U + 0.5 * alpha2 * d_E};

    vector<double> dissip(3);
    for(int i=0; i<3; i++) {
        for (int j = 0; j < 3; j++){
            dissip[i] += l[j] * alpha[j] * v[j][i];
        }
        
    }

    return dissip;
}

vector<double> flux_roe(vector<double> w_L, vector<double> w_R) {
    vector<double> phi(3);
    for(int i=0; i<3; i++) {
        phi[i] = 0.5 * (flux_centre(w_L)[i] + flux_centre(w_R)[i] - dissipation_roe(w_L, w_R)[i]);
    }
    return phi;
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

vector<vector<double>> lax_wendroff_artificiel(vector<vector<double>> W, double T, double dx, double D, double cfl=0.95) {
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
                W_tilde[j][i] = (W[j][i]+W[j][i+1])/2.0 - dt/(2*dx) * (f_right[j] - f_left[j] + D * dx * (2*W[j][i] - W[j][i+1] - W[j][i-1]));
            }
        }
        for(int i=1; i<N-1; i++) {
            auto f_tilde_right = flux_centre({W_tilde[0][i], W_tilde[1][i], W_tilde[2][i]});
            auto f_tilde_left  = flux_centre({W_tilde[0][i-1], W_tilde[1][i-1], W_tilde[2][i-1]});
            for(int j=0; j<3; j++) {
                W_new[j][i] = W[j][i] - dt/dx * (f_tilde_right[j] - f_tilde_left[j]+ D * dx * (2*W_tilde[j][i] - W_tilde[j][i+1] - W_tilde[j][i-1]));
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

vector<vector<double>> maccormack_artificiel(vector<vector<double>> W, double T, double dx, double D, double cfl=0.95) {
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
                W_predict[j][i] = W[j][i] - dt/dx * (f_right[j] - f_left[j] + D * dx * (2*W[j][i] - W[j][i+1] - W[j][i-1]));
            }
        }

        for(int i=1; i<N-1; i++) {
            auto f_predict_right = flux_centre({W_predict[0][i], W_predict[1][i], W_predict[2][i]});
            auto f_predict_left  = flux_centre({W_predict[0][i-1], W_predict[1][i-1], W_predict[2][i-1]});
            for(int j=0; j<3; j++) {
                W_correct[j][i] = 0.5 * (W[j][i] + W_predict[j][i]) - 0.5 * dt/dx * (f_predict_right[j] - f_predict_left[j] + D * dx * (2*W_predict[j][i] - W_predict[j][i+1] - W_predict[j][i-1]));
            }
        }
        W = W_correct;
        t += dt;
    }
    return W;
}

vector<vector<double>> roe(vector<vector<double>> W, double T, double dx, double cfl=0.95) {
    int N = W[0].size();
    double t = 0.0;
    vector<vector<double>> W_new = W;
    while(t < T) {
        double dt = calc_dt(W, dx, cfl);
        if (t + dt > T) dt = T - t;

        for(int i=1; i<N-1; i++) {
            auto phi_right = flux_roe({W[0][i], W[1][i], W[2][i]}, {W[0][i+1], W[1][i+1], W[2][i+1]});
            auto phi_left  = flux_roe({W[0][i-1], W[1][i-1], W[2][i-1]}, {W[0][i], W[1][i], W[2][i]});
            for(int j=0; j<3; j++) {
                W_new[j][i] = W[j][i] - dt/dx * (phi_right[j] - phi_left[j]);
            }
        }
        W = W_new;
        t += dt;
    }
    return W;
}

int main() {

    // Paramètres ex 10.1
    int n = 101;
    double x0 = 0.5;
    double t = 0.2;
    double D = 8.0;
    double dx = 1.0/(n-1);

    vector<vector<double>> W(3, vector<double>(n));
    for(int i=0; i<n; i++) {
        W[0][i] = (i < n/2) ? rho_L : rho_R;
        W[1][i] = (i < n/2) ? W[0][i] * U_L : W[0][i] * U_R;
        W[2][i] = (i < n/2) ? (p_L/(g-1) + 0.5 * W[1][i] * U_L) : (p_R/(g-1) + 0.5 * W[1][i] * U_R);
    }

    vector<vector<double>> result = roe(W, t, dx);

    for(int i=0; i<n; i++) {
        result[1][i] /= result[0][i]; // Convertir rho_U en U
        result[2][i] = (g - 1) * (result[2][i] - 0.5 * result[1][i] * result[1][i] / result[0][i]); // Convertir E en p
    }


    ofstream file("output_num.csv");
    file << "x,rho,U,p\n";
    for(int i = 0; i < n; ++i) {
        file << i * dx << ',' << result[0][i] << ',' << result[1][i] << ',' << result[2][i] << '\n';
    }
    
    file.close();
    return 0;
}