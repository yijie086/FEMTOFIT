#include"newBMKmodel.h"
#include<iostream>

void newBMK(){
    
    double EB = 7.546;

    // Define kinematics: xB, Q2, t, phi
    double xB = 0.16;
    double Q2 = 1.33;
    double t = -0.21;
    double phi = 90.0;  // degrees

    // Beam and target polarizations
    double q_beam = -1;   // electron
    double L_beam = 1;    // +1 or -1 for helicity
    double L_target = 0;  // unpolarized

    // Create the DVCS object
    BMK_DVCS dvcs(q_beam, L_beam, L_target, EB, xB, Q2, t, phi);
    //BMK_DVCS dvcs_FX(q_beam, L_beam, L_target, EB, xB, Q2, t, phi);

    // Optional: enable verbose output
    // dvcs.VERB = true;

    // Compute observables
    double cross_sec = dvcs.CrossSection(); // nb/sr/GeV^4 or so
    double bsa = dvcs.BSA();                // Beam Spin Asymmetry

    //double cross_sec_FX = dvcs_FX.CrossSection(); // nb/sr/GeV^4 or so
    //double bsa_FX= dvcs_FX.BSA();                // Beam Spin Asymmetry


    std::cout << "DVCS Kinematics:\n";
    std::cout << "xB = " << xB << ", Q2 = " << Q2 << " GeV^2, t = " << t << " GeV^2, phi = " << phi << " deg\n";
    std::cout << "Unpolarized Cross Section = " << cross_sec  <<"nb\n";
    std::cout << "Beam Spin Asymmetry (BSA) = " << bsa << "\n";
}