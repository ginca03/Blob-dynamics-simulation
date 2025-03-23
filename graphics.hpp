#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include <TImage.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TPad.h>
#include <TLegend.h>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include "simulation.hpp"

using namespace std;

// Function prototypes for graphics
void plot_profiles(vector<vector<double>>& density, vector<vector<double>>& potential, vector<vector<double>>& vorticity);
void updateGraph(TH2D* hist, TPad* pad, const vector<vector<vector<double>>>& history, int currentTime, const string& variableName);
void createGif(const string& outputGifPath, const string& tempPath);
void animate(const vector<vector<vector<double>>>& densityHistory,
             const vector<vector<vector<double>>>& potentialHistory,
             const vector<vector<vector<double>>>& vorticityHistory);

//Function to plot the profiles of density, potential, and vorticity in the same canvas
inline void plot_profiles(vector<vector<double>>& density,
                   vector<vector<double>>& potential,
                   vector<vector<double>>& vorticity) {

    // Create vectors for x and y profiles
    vector<double> x_coords(Nx), y_coords(Ny);
    for (int i = 0; i < Nx; ++i) x_coords[i] = i * dx;
    for (int j = 0; j < Ny; ++j) y_coords[j] = j * dy;

    // Extract profiles along the middle of the domain
    int mid_y = Ny / 2; // Middle row for x-profile
    int mid_x = Nx / 2; // Middle column for y-profile

    vector<double> density_x(Nx), potential_x(Nx), vorticity_x(Nx);
    vector<double> density_y(Ny), potential_y(Ny), vorticity_y(Ny);

    for (int i = 0; i < Nx; ++i) {
        density_x[i] = density[i][mid_y];
        potential_x[i] = potential[i][mid_y];
        vorticity_x[i] = vorticity[i][mid_y];
    }

    for (int j = 0; j < Ny; ++j) {
        density_y[j] = density[mid_x][j];
        potential_y[j] = potential[mid_x][j];
        vorticity_y[j] = vorticity[mid_x][j];
    }

    // Determine global min and max values for x and y profiles
    double global_min_x = min({*min_element(density_x.begin(), density_x.end()),
                               *min_element(potential_x.begin(), potential_x.end()),
                               *min_element(vorticity_x.begin(), vorticity_x.end())});
    double global_max_x = max({*max_element(density_x.begin(), density_x.end()),
                               *max_element(potential_x.begin(), potential_x.end()),
                               *max_element(vorticity_x.begin(), vorticity_x.end())});

    double global_min_y = min({*min_element(density_y.begin(), density_y.end()),
                               *min_element(potential_y.begin(), potential_y.end()),
                               *min_element(vorticity_y.begin(), vorticity_y.end())});
    double global_max_y = max({*max_element(density_y.begin(), density_y.end()),
                               *max_element(potential_y.begin(), potential_y.end()),
                               *max_element(vorticity_y.begin(), vorticity_y.end())});

    // Create canvas and divide it into two pads
    TCanvas* c = new TCanvas("c", "Profiles", 600, 1200);
    c->Divide(1, 2); // Two pads side by side

    // Plot x-profile in the first pad
    c->cd(1); // Switch to the first pad
    TGraph* g_density_x = new TGraph(Nx, x_coords.data(), density_x.data());
    TGraph* g_potential_x = new TGraph(Nx, x_coords.data(), potential_x.data());
    TGraph* g_vorticity_x = new TGraph(Nx, x_coords.data(), vorticity_x.data());

    g_density_x->SetLineColor(kRed);
    g_density_x->SetLineWidth(2);
    g_density_x->SetTitle("X-Profile;X;Natural Units");
    g_density_x->GetYaxis()->SetRangeUser(global_min_x, global_max_x); // Adjust y-axis range
    g_density_x->Draw("AL");

    g_potential_x->SetLineColor(kBlue);
    g_potential_x->SetLineWidth(2);
    g_potential_x->Draw("L SAME");

    g_vorticity_x->SetLineColor(kGreen + 2);
    g_vorticity_x->SetLineWidth(2);
    g_vorticity_x->Draw("L SAME");

    TLegend* legend_x = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend_x->AddEntry(g_density_x, "Density", "l");
    legend_x->AddEntry(g_potential_x, "Potential", "l");
    legend_x->AddEntry(g_vorticity_x, "Vorticity", "l");
    legend_x->Draw();

    // Plot y-profile in the second pad
    c->cd(2); // Switch to the second pad
    TGraph* g_density_y = new TGraph(Ny, y_coords.data(), density_y.data());
    TGraph* g_potential_y = new TGraph(Ny, y_coords.data(), potential_y.data());
    TGraph* g_vorticity_y = new TGraph(Ny, y_coords.data(), vorticity_y.data());

    g_density_y->SetLineColor(kRed);
    g_density_y->SetLineWidth(2);
    g_density_y->SetTitle("Y-Profile;Y;Natural Units");
    g_density_y->GetYaxis()->SetRangeUser(global_min_y, global_max_y); // Adjust y-axis range
    g_density_y->Draw("AL");

    g_potential_y->SetLineColor(kBlue);
    g_potential_y->SetLineWidth(2);
    g_potential_y->Draw("L SAME");

    g_vorticity_y->SetLineColor(kGreen + 2);
    g_vorticity_y->SetLineWidth(2);
    g_vorticity_y->Draw("L SAME");

    TLegend* legend_y = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend_y->AddEntry(g_density_y, "Density", "l");
    legend_y->AddEntry(g_potential_y, "Potential", "l");
    legend_y->AddEntry(g_vorticity_y, "Vorticity", "l");
    legend_y->Draw();

    // Save the combined canvas
    c->SaveAs("init_profiles.png");
}

// Create a GIF
inline void createGif(const string& outputGifPath, const string& tempPath) {
    string command = "convert -delay 10 -loop 0 " + tempPath + "/frame_*.png " + outputGifPath;
    system(command.c_str());
}

// Update the graph with dynamic title
inline void updateGraph(TH2D* hist, TPad* pad, const vector<vector<vector<double>>>& history, int currentTime, const string& variableName) {
    hist->Reset();
    for (int ix = 0; ix < Nx; ++ix) {
        for (int iy = 0; iy < Ny; ++iy) {
            hist->SetBinContent(ix + 1, iy + 1, history[currentTime][ix][iy]);
        }
    }
    pad->cd(); // Switch to the current pad
    string title = Form("%s at time = %.4f", variableName.c_str(), currentTime * dt);
    hist->SetTitle(title.c_str());
    hist->Draw("COLZ");
}

// Animate with dynamic title
inline void animate(const vector<vector<vector<double>>>& densityHistory,
             const vector<vector<vector<double>>>& potentialHistory,
             const vector<vector<vector<double>>>& vorticityHistory) {
    TApplication app("Animation", nullptr, nullptr);

    TCanvas* canvas = new TCanvas("canvas", "Simulation Results", 600, 1200); // Greater height for vertical layout
    canvas->Divide(1, 3); // Divide the canvas into three rows

    // Histograms for density, potential, and vorticity
    TH2D* densityHist = new TH2D("densityHist", "Density", Nx, 0, Lx, Ny, 0, Ly);
    TH2D* potentialHist = new TH2D("potentialHist", "Potential", Nx, 0, Lx, Ny, 0, Ly);
    TH2D* vorticityHist = new TH2D("vorticityHist", "Vorticity", Nx, 0, Lx, Ny, 0, Ly);

    densityHist->SetStats(false);
    potentialHist->SetStats(false);
    vorticityHist->SetStats(false);

    //Set titles
    densityHist->SetXTitle("x");    densityHist->SetYTitle("y");
    potentialHist->SetXTitle("x");  potentialHist->SetYTitle("y");
    vorticityHist->SetXTitle("x");  vorticityHist->SetYTitle("y");
    
    TImage* img = TImage::Create();

    int numFrames = animationDuration * fps;
    int frameInterval = Nt / numFrames;

    system("mkdir -p frames"); // Create a folder for frames

    for (int i = 0; i < numFrames; ++i) {
        int currentTime = i * frameInterval;
        if (currentTime >= Nt) currentTime = Nt - 1;

        // Update each graph with the current time
        updateGraph(densityHist, (TPad*)canvas->cd(1), densityHistory, currentTime, "Density");
        updateGraph(potentialHist, (TPad*)canvas->cd(2), potentialHistory, currentTime, "Potential");
        updateGraph(vorticityHist, (TPad*)canvas->cd(3), vorticityHistory, currentTime, "Vorticity");

        // Save the frame
        canvas->Update();
        string frameName = Form("frames/frame_%03d.png", i);
        img->FromPad(canvas);
        img->WriteImage(frameName.c_str());
    }

    createGif("output.gif", "frames");
    //system("rm -r frames"); // Remove temporary frames

    gApplication->Terminate();
}

#endif // GRAPHICS_HPP
