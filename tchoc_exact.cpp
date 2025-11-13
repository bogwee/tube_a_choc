#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

// Paramètres de Sod
const double g = 1.4;

// Zone Gauche
double p_L = 8.0;
double rho_L = 10.0/g;
double a_L = sqrt(g * p_L / rho_L);
double U_L = 0.0;

// Zone Droite
double p_R = 1.0;
double rho_R = 1.0/g;
double a_R = sqrt(g * p_R / rho_R);
double U_R = 0.0;

double fms(double Ms) {
    return a_L * ((g + 1)/(g - 1)) * (1 - pow((p_R/p_L)*((2*g/(g+1))*Ms*Ms - (g-1)/(g+1)), (g-1)/(2*g))) + 1/Ms;
}

double point_fixe(double x0, double tol=1e-6, int max_iter=2000) {
    double Msk = x0;
    double Msk_1;
    for(int iter = 0; iter < max_iter; iter++) {
        Msk_1 = fms(Msk);
        if(fabs(Msk_1 - Msk) < tol) {
            return Msk_1;
        }
        Msk = Msk_1;
    }
    cerr << "Warning: point_fixe did not converge" << endl;
    return Msk_1;
}

vector<vector<double>> tchoc_exact(vector<double> x, double x0, double t) {
    vector<vector<double>> result;
    double Ms = point_fixe(2.0);

    // Zone 1 (après onde de choc)
    double rho_1 = rho_R * ((g + 1) * Ms * Ms) / ((g - 1) * Ms * Ms + 2);
    double U_1 = 2/(g+1)*(Ms - 1/Ms);
    double p_1 = p_R * (2*g/(g+1)*Ms*Ms - (g-1)/(g+1));
    double a_1 = sqrt(g * p_1 / rho_1);
    
    // Zone 2 (entre onde de choc et contact)
    double U_2 = U_1;
    double p_2 = p_1;
    double rho_2 = rho_L * pow((p_2/p_L), (1.0/g));
    double a_2 = a_L - (g - 1)/2 * U_2;

    // Positions des discontinuités
    double x_1 = x0 - a_L * t; // Front de l'onde rarefaction
    double x_2 = x0 + (U_2 - a_2) * t; // Fin de l'onde rarefaction
    double x_3 = x0 + U_2 * t; // Contact
    double x_4 = x0 + Ms * t; // Onde de choc

    // Remplissage des résultats
    vector<double> rho(x.size()), U(x.size()), p(x.size()), a(x.size());
    for(size_t i = 0; i < x.size(); i++) {
        if(x[i] < x_1) {
            U[i] = U_L;
            a[i] = a_L;
            p[i] = p_L;
            rho[i] = rho_L;
            continue;
        } else if(x[i] < x_2) {
            U[i] = 2/(g+1)*(a_L + (x[i]-x0)/t);
            a[i] = 2/(g+1)*(a_L - (g - 1)/2 * (x[i]-x0)/t);
            p[i] = p_L * pow((a[i]/a_L), (2.0*g/(g-1)));
            rho[i] = g * p[i] / (a[i] * a[i]);
            continue;
        } else if(x[i] < x_3) {
            U[i] = U_2;
            a[i] = a_2;
            p[i] = p_2;
            rho[i] = rho_2;
            continue;
        } else if(x[i] < x_4) {
            U[i] = U_1;
            a[i] = a_1;
            p[i] = p_1;
            rho[i] = rho_1;
            continue;
        } else {
            U[i] = U_R;
            a[i] = a_R;
            p[i] = p_R;
            rho[i] = rho_R;
        }
    }
    result.push_back(rho);
    result.push_back(U);
    result.push_back(p);
    return result;
}


int main() {

    // Paramètres ex 10.1
    int n = 81;
    double x0 = 0.5;
    double t = 0.2;

    vector<double> x(n);
    for(int i=0; i<n; i++) {
        x[i] = i * (1.0 / (n - 1));
    }

    vector<vector<double>> result = tchoc_exact(x, x0, t);


    ofstream file("output.csv");
    file << "x,rho,U,p\n";
    for(int i = 0; i < n; ++i) {
        file << x[i] << ',' << result[0][i] << ',' << result[1][i] << ',' << result[2][i] << '\n';
    }
    
    file.close();
    return 0;
}