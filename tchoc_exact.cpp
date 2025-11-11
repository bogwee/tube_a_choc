#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double gamma;

struct Zone_R {
    double p = 1;
    double rho = 1/gamma;
    double T = 1/gamma;
    double U = 0;
};

struct Zone_L {
    double p;
    double rho;
    double T;
    double U = 0;
};


vector<vector<double>> tchoc_exact(vector<double> x, double x0, double t) {
    vector<vector<double>> result;
    return result;
}

int main() {
    int n = 81;
    vector<double> x(n);
    for(int i=0; i<n; i++) {
        x[i] = i * (1.0 / (n - 1));
    }

    double x0 = 0.5;
    double t = 0.2;
    gamma = 1.4;
    Zone_L zone_L;
    zone_L.p = 8.0;
    zone_L.rho = 10.0/gamma;
    vector<vector<double>> result = tchoc_exact(x, x0, t);
    return 0;
}