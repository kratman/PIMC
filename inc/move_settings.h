
#pragma once

//Compile-time move options
double step = 0.15; //Initial step size for the beads
const double StepMin = 0.075; //Minimum step size
const double StepMax = 1.0; //Maximum step size
const double Centratio= 10.0; //Scales 'step' for centroids
const int Acc_Check = 5000; //Eq steps before checking accratio
const double RadMin = 0.01; //Minimum electron radius
const double RadMax = 25.0; //Maximum electron radius

//Move Probabilities
//Note: These probabilities allow for multi-particle moves
double BeadProb = 0.55; //Probability to move a single bead
double CentProb = 0.55; //Probability to move a centroid
double ElBeadProb = 0.25; //Probability to move an electron bead
double ElCentProb = 0.25; //Probability to move an electron centroid
double RadProb = 0; //Probability to change electron radius
double SwapProb = 0.05; //Probability to swap spins
double FlipProb = 0.05; //Probability to flip a single spin
double VolProb = 0.05; //Volume change probability
