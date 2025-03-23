#include <TImage.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TApplication.h>
#include <iostream>
#include "simulation.hpp"
#include "graphics.hpp"

int main() {
    // Initialize the blob
    initialize_blob();
    
    // Alternatively, you can use a wave initialization
    // initialize_wave(4.0, 4.0, 0.1);
    
    // Solve the simulation
    solve();
    
    // Animate the results
    animate(densityHistory, potentialHistory, vorticityHistory);
    
    return 0;
}

void blob() {
    initialize_blob();
    //initialize_wave(4.0, 4.0, 0.1);
    solve();
    animate(densityHistory, potentialHistory, vorticityHistory);
}
