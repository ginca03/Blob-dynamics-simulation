#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

// Physics paramaters
const double Lk = 5467.5;
const double g = 2.44e-3;
const double n0 = 0.8e13; //cm^-3
const double B = 0.5; // Tesla
const double D = 0.0015;
const double mu = 0.04;

// Simulation parameters
const int Nx = 100; //x resolution
const int Ny = 100; //y resolution
const double Lx = 20.0;
const double Ly = 10.0;
const double dx = Lx / (Nx - 1);
const double dy = Ly / (Ny - 1);
const double dt = 1e-5;
const int Nt = 1e4;
const double minDensity = 1e-4;

// Animation parameters
const int sampleRate = 100;
const int animationDuration = 5;
const int fps = 30;

// Global variables
extern vector<vector<vector<double>>> densityHistory;
extern vector<vector<vector<double>>> potentialHistory;
extern vector<vector<vector<double>>> vorticityHistory;
extern vector<vector<double>> n;
extern vector<vector<double>> phi;
extern vector<vector<double>> Omega;

// Initialize these variables (in the source file that includes this header)
vector<vector<vector<double>>> densityHistory(Nt, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));
vector<vector<vector<double>>> potentialHistory(Nt, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));
vector<vector<vector<double>>> vorticityHistory(Nt, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));
vector<vector<double>> n(Nx, vector<double>(Ny, 0.0));
vector<vector<double>> phi(Nx, vector<double>(Ny, 0.0));
vector<vector<double>> Omega(Nx, vector<double>(Ny, 0.0));

// Function prototypes for physics simulation
double laplacian(const vector<vector<double>>& field, int i, int j);
void updatePotential(int t = 0);
double drag(const vector<vector<double>>& field, int i, int j);
vector<vector<double>> rotate_potential(vector<vector<double>>& phi, double theta);
void initialize_blob();
void initialize_wave(double frac_x, double frac_y, double amplitude);
void solve();
void blob();

// Use finite differince for computing laplacian
inline double laplacian(const vector<vector<double>>& field, int i, int j) {
    if(i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1){
        return 0.0;}  // Assume the laplacian is zero on the border

    double lap_x = (field[i+1][j] - 2 * field[i][j] + field[i-1][j]) / (dx * dx);
    double lap_y = (field[i][j+1] - 2 * field[i][j] + field[i][j-1]) / (dy * dy);
    return lap_x + lap_y;
}

// Invert locally ∇^2Ω = ɸ using Jacobi matrix inversion for laplacian
inline void updatePotential(int t) {
    int max_iter = 10000; double tolerance = 1e-4;
    // Temporary matrix to store updated values of phi
    vector<vector<double>> phi_new(Nx, vector<double>(Ny, 0.0));
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0; // Maximum difference for convergence check

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {

                auto p = [&](int x) { return (x + Nx) % Nx; };
                
                // Jacobi update rule
                phi_new[i][j] = 0.25 * (phi[p(i+1)][p(j)] + phi[p(i-1)][p(j)]
                                     + phi[p(i)][p(j+1)] + phi[p(i)][p(j-1)]
                                     - dx * dx * Omega[i][j]);

                // Calculate the maximum difference between iterations, using max norm for computational efficiency (all norms converge in finite dimensions)
                max_diff = max(max_diff, abs(phi_new[i][j] - phi[i][j]));
            }
        }

        phi = phi_new;

        // Check for convergence
        if (max_diff < tolerance) {
            //cout << "Converged in " << iter + 1 << " iterations." << endl;
            return;
        }
    }

    // Warn if the method did not converge within the maximum iterations
    cerr << "Warning: Jacobi method did not converge at timestep " << t <<" within " << max_iter << " iterations." << endl;
}

// Function to rotate the electric field by an angle theta
inline vector<vector<double>> rotate_potential(vector<vector<double>>& phi, double theta) {
    vector<vector<double>> phi_rotated(Nx, vector<double>(Ny, 0.0));
    
    auto p = [&](int i) { return (i + Nx) % Nx; };

    auto gradient = [&](int i, int j) {
        // Compute partial derivatives using central differences with periodic boundary conditions
        double dphidx = (phi[p(i + 1)][j] - phi[p(i - 1)][j]) / (2 * dx);
        double dphidy = (phi[i][p(j + 1)] - phi[i][p(j - 1)]) / (2 * dy);
        return make_pair(dphidx, dphidy);
    };

    // Construct the rotated potential
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            // Compute the original gradient at point (i, j)
            auto [dphidx, dphidy] = gradient(i, j);

            // Rotate the electric field components
            double Ex_rot = dphidx * cos(theta) - dphidy * sin(theta);
            double Ey_rot = dphidx * sin(theta) + dphidy * cos(theta);

            // Integrate the rotated gradient to compute the new potential
            phi_rotated[i][j] = -Ex_rot * dx - Ey_rot * dy;
        }
    }
    return phi_rotated;
}

// Initialise density, potential and vorticity
inline void initialize_blob() {
    //parameters fot init blob
    double x0 = Lx / 2.0;
    double y0 = Ly / 2.0;
    double sigma_x = Lx / 20.0 *0.5;
    double sigma_y = Ly / 10.0;
    
    //use periodic indexing
    if(Nx != Ny){cerr << "If the grid is not a square fix periodic indexing\n";}
    auto p = [&](int i){return (i + Nx) % Nx;}; //p(i) is i mod Nx

    //init density as gaussian distr.
    for (int i = 0; i < Nx; ++i) {for (int j = 0; j < Ny; ++j) {
            double x = p(i) * dx;
            double y = p(j) * dy;
            n[p(i)][p(j)] = 1 + exp(-((x - x0) * (x - x0)) / (2 * sigma_x * sigma_x)
                          -((y - y0) * (y - y0)) / (2 * sigma_y * sigma_y));
        }}
    
    //init potential using eq 2: dΩ/dt=0 at t=0
    for (int i = 0; i < Nx; ++i) {for (int j = 0; j < Ny; ++j) {
            phi[p(i)][p(j)] = -Lk*g / n[p(i)][p(j)] * (n[p(i)][p(j+1)] - n[p(i)][p(j-1)]) / (2*dy); //assuming vorticity at a minimum at t=0 dΩ/dt=0;
        }}
    
    //rotate potential by 30 degrees
    //phi = rotate_potential(phi, -M_PI/3);
    
    //init vorticity as ∇^2(Phi)=Ω
    for (int i = 0; i < Nx; ++i) {for (int j = 0; j < Ny; ++j) {
        Omega[p(i)][p(j)] = laplacian(phi, p(i), p(j));
        }}
    
    //plot initialisation profiles
    //plot_profiles(n, phi, Omega);
}

// Initialise density, potential and vorticity with a sinusoidal wave
inline void initialize_wave(double frac_x, double frac_y, double amplitude) {
    // Compute wave numbers based on grid fraction
    double kx = 2 * M_PI / (Lx / frac_x); // Wave number in x
    double ky = 2 * M_PI / (Ly / frac_y); // Wave number in y
    
    // Use periodic indexing for Dirichelet periodic conditions
    if (Nx != Ny) { cerr << "If the grid is not a square fix periodic indexing\n"; }
    auto p = [&](int i) { return (i + Nx) % Nx; }; // p(i) is i mod Nx

    // Initialize density with a sinusoidal wave
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = p(i) * dx;
            double y = p(j) * dy;
            n[p(i)][p(j)] = 1 + amplitude * sin(kx * x) * sin(ky * y);
        }
    }
    
    // Initialize potential phi
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            phi[p(i)][p(j)] = -Lk * g / n[p(i)][p(j)] * (n[p(i)][p(j + 1)] - n[p(i)][p(j - 1)]) / (2 * dy); // Assuming dΩ/dt = 0 initially
        }
    }

    // Initialize vorticity Omega
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Omega[p(i)][p(j)] = laplacian(phi, p(i), p(j));
        }
    }
}

//Compute the term [(b/B)⨯∇ɸ] * ∇ = (-E_y * d/dx + E_x * d/dy] / B = drag(field)
inline double drag(const vector<vector<double>>& field, int i, int j){
    
    auto p = [&](int i){return (i + Nx) % Nx;}; //p(i) is i mod Nx
    
    double dPhi_dx = (phi[p(i+1)][p(j)] - phi[p(i-1)][p(j)]) / (2*dx); //central differnce scheme
    double dPhi_dy = (phi[p(i)][p(j+1)] - phi[p(i)][p(j-1)]) / (2*dy);
    
    double dField_dx = (field[p(i+1)][p(j)] - field[p(i-1)][p(j)]) / (2*dx);
    double dField_dy = (field[p(i)][p(j+1)] - field[p(i)][p(j-1)]) / (2*dy);
    
    return (-dPhi_dy * dField_dx + dPhi_dx * dField_dy) / B;
}

// Solve equations with Runge-Kutta 4th order and diffusion
inline void solve() {
    for (int t = 0; t < Nt; ++t) {
        vector<vector<double>> n_next = n;
        vector<vector<double>> Omega_next = Omega;
        
        
        // Check CFL condition for diffusion stability
        double dt_diff_n = dx * dx / (4.0 * D);
        double dt_diff_Omega = dx * dx / (4.0 * mu);
        if (dt > dt_diff_n || dt > dt_diff_Omega) {
            cerr << "Time step dt exceeds CFL conditions for diffusion. Reduce dt.\n";
            return;
        }

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                // Avoid dividing by zero
                double nSafe = max(n[i][j], minDensity);

                auto p = [&](int x) { return (x + Nx) % Nx; };
                
                // Central difference for derivatives at (i, j)
                auto compute_der = [&](double n_val, double Omega_val, int i, int j) {
                    
                    //laplacians for diffusion
                    double lap_n = laplacian(n,i,j);
                    double lap_Omega = laplacian(Omega,i,j);
                    
                    // Central difference for derivatives at (i, j)
                    double dPhi_dy = (phi[p(i)][p(j + 1)] - phi[p(i)][p(j - 1)]) / (2 * dy);
                    double dn_dy = (n[p(i)][p(j + 1)] - n[p(i)][p(j - 1)]) / (2 * dy);
                    
                    //caluclate derivatives in time using euler method
                    double dn_dt = drag(n, i, j)
                                   - (n_val * phi[p(i)][p(j)] / Lk)
                                   + (g * n_val * dPhi_dy)
                                   - (g * dn_dy)
                                   + (D * lap_n);

                    double dOmega_dt = drag(Omega, i, j)
                                       + phi[p(i)][p(j)] / Lk
                                       + (nSafe / g * dn_dy)
                                       + (mu * lap_Omega);

                    return make_pair(dn_dt, dOmega_dt);
                };

                // Runge-Kutta coefficients for n and Omega
                double k1_n, k1_Omega, k2_n, k2_Omega, k3_n, k3_Omega, k4_n, k4_Omega;

                // k1
                tie(k1_n, k1_Omega) = compute_der(n[i][j], Omega[i][j], i, j);

                // k2
                tie(k2_n, k2_Omega) = compute_der(
                    n[i][j] + 0.5 * dt * k1_n, Omega[i][j] + 0.5 * dt * k1_Omega, i, j);

                // k3
                tie(k3_n, k3_Omega) = compute_der(
                    n[i][j] + 0.5 * dt * k2_n, Omega[i][j] + 0.5 * dt * k2_Omega, i, j);

                // k4
                tie(k4_n, k4_Omega) = compute_der(
                    n[i][j] + dt * k3_n, Omega[i][j] + dt * k3_Omega, i, j);

                n_next[i][j] = n[i][j] + (dt / 6.0) * (k1_n + 2 * k2_n + 2 * k3_n + k4_n);
                Omega_next[i][j] = Omega[i][j] + (dt / 6.0) * (k1_Omega + 2 * k2_Omega + 2 * k3_Omega + k4_Omega);
            }
        }

        n = n_next;
        Omega = Omega_next;

        updatePotential(t);

        // Save values of n, phi, and Omega in history
        densityHistory[t] = n;
        potentialHistory[t] = phi;
        vorticityHistory[t] = Omega;

        if (t % 1000 == 0) {
            cout << "Step " << t << " - Density center: " << n[Nx / 2][Ny / 2] << endl;
        }
    }
}

#endif // SIMULATION_HPP
