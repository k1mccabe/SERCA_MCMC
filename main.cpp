/*-----------------------------------------------------------------------------------------------------
//                     University of California, San Diego
//                      Dept. of Chemistry & Biochemistry
//-----------------------------------------------------------------------------------------------------
// Authors: Sophia P. Hirakis & Kimberly J. McCabe
// Year  :  4/2018
//-----------------------------------------------------------------------------------------------------
// This code uses the "Particle Swarm Optimization" algorithm to optimize the rates in a model of the SERCA pump.
//
//
// The SERCA model is based on the published model by Inesi (1988)
//
//          [S1]       [S2]             [S3]                 [S4]                 [S5]
//          E.Ca <==> E'.Ca  + Ca <==> E'.Ca2 (+ ATP) <==> E'.ATP.Ca2  <==>   E'~P.ADP.Ca2
//           /\                                                                  //  \\
//           ||                                                                 //    \\
//           ||                                                          [S6a]  //      \\  [S6]
//     +Ca   ||                                                      *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
//           ||                                                                \\      //
//           ||                                                        (+ ADP)  \\    //
//           \/                                                                  \\  //
//    (Pi +) E <==> *E-Pi <==> *E-P + Ca <==> *E-P.Ca  <==> *E'-P.Ca + Ca <==>  *E'-P.Ca2
//          [S0]    [S11]      [S10]           [S9]          [S8]                 [S7]
//
// State Reaction               State Product          Rate(f)   Rate(r)
//  S0   E + Ca             <==> S1   E.Ca             k_S0_S1   k_S1_S0
//  S1   E.Ca               <==> S2   E'.Ca            k_S1_S2   k_S2_S1
//  S2   E'.Ca + Ca         <==> S3   E'.Ca2           k_S2_S3   k_S3_S2
//  S3   E'.Ca2 (+ ATP)     <==> S4   E'.ATP.Ca2       k_S3_S4   k_S4_S3
//  S4   E'.ATP.Ca2         <==> S5   E'~P.ADP.Ca2     k_S4_S5   k_S5_S4
//  S5   E'~P.ADP.Ca2       <==> S6a *E'-P.ADP.Ca2     k_S5_S6a  k_S6a_S5
//  S6a *E'-P.ADP.Ca2       <==> S7  *E'-P.Ca2 (+ ADP) k_S6a_S7  k_S7_S6a
//  S5   E'~P.ADP.Ca2       <==> S6   E'~P.Ca2 (+ ADP) k_S5_S6   k_S6_S5
//  S6   E'~P.Ca2           <==> S7  *E'-P.Ca2         k_S6_S7   k_S7_S6
//  S7  *E'-P.Ca2           <==> S8  *E-P.Ca   + Ca    k_S7_S8   k_S8_S7
//  S8  *E'-P.Ca            <==> S9  *E-P.Ca           k_S8_S9   k_S9_S8
//  S9  *E-P.Ca             <==> S10 *E-P      + Ca    k_S9_S10  k_S10_S9
//  S10 *E-P                <==> S11 *E-Pi             k_S10_S11 k_S11_S10
//  S11 *E-Pi               <==> S0   E + (Pi)         k_S11_S0  k_S0_S11
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>
//#include </usr/include/openmpi-x86_64/mpi.h>
#include "get_Residual.h"
#include "lastRun.h"

using namespace std;

const int n_particles_PSO = 100;

int   n_s;                 // Number of states
int   n_pCa ;              // Number of simulated pCa or Ca values
int   n_SERCA_Molecules;   // Max number used to repeat the simulation
int   max_tsteps;          // Max number of time stepping
float dt;                  // fixed time step
float residual_cost_func[n_particles_PSO]; // to track the residual between numerics and experiments
float Res_pbest[n_particles_PSO];
int   id, p, ierr, argc; // for parallel    
//---------------------------------------------
// model reference parameters that we need to optimize
//--------------------------------------------
//float Ca_cyt_conc, Ca_cyt_conc_lower , Ca_cyt_conc_upper ; //pCa
//float MgATP_conc , MgATP_conc_lower  ,  MgATP_conc_upper ;
float k_S0_S1    ,  k_S0_S1_lower    ,  k_S0_S1_upper   ; // First Ca Association
float k_S2_S3    ,  k_S2_S3_lower    ,  k_S2_S3_upper   ; // Second Ca Association
float k_S7_S8    ,  k_S7_S8_lower    ,  k_S7_S8_upper  ; // First Ca Dissociation
float k_S9_S10  ,  k_S9_S10_lower  ,  k_S9_S10_upper ; // Second Ca Dissociation
//float k_S5_S6a    ,  k_S5_S6a_lower    ,  k_S5_S6a_upper ; // Forward Phosphorylation
//float k_S6_S7    ,  k_S6_S7_lower    ,  k_S6_S7_upper ; // Forward Phosphorylation
//float k_S0_S11   ,  k_S0_S11_lower   ,  k_S0_S11_upper ; // Reverse Phosphorylation
 

//-----------------------------------------------
// Particle Swarm Optimization (PSO) parameters
//-----------------------------------------------
//float X_Ca_cyt_conc_PSO         [n_particles_PSO] , V_Ca_cyt_conc_PSO         [n_particles_PSO];
float X_k_S0_S1_PSO   [n_particles_PSO] , V_k_S0_S1_PSO    [n_particles_PSO];
float X_k_S2_S3_PSO   [n_particles_PSO] , V_k_S2_S3_PSO    [n_particles_PSO];
float X_k_S7_S8_PSO   [n_particles_PSO] , V_k_S7_S8_PSO    [n_particles_PSO];
float X_k_S9_S10_PSO  [n_particles_PSO] , V_k_S9_S10_PSO   [n_particles_PSO];
//float X_k_S5_S6a_PSO   [n_particles_PSO] , V_k_S5_S6a_PSO    [n_particles_PSO];
//float X_k_S6_S7_PSO   [n_particles_PSO] , V_k_S6_S7_PSO    [n_particles_PSO];
//float X_k_S0_S11_PSO  [n_particles_PSO] , V_k_S0_S11_PSO   [n_particles_PSO];

//-------------------------------------------------------------------
// Parameters used for global best (gbest) and personal best (pbest)
//-------------------------------------------------------------------
//float Ca_cyt_conc_gbest         , Ca_cyt_conc_pbest         [n_particles_PSO];
float k_S0_S1_gbest   , k_S0_S1_pbest   [n_particles_PSO];
float k_S2_S3_gbest   , k_S2_S3_pbest   [n_particles_PSO];
float k_S7_S8_gbest   , k_S7_S8_pbest   [n_particles_PSO];
float k_S9_S10_gbest , k_S9_S10_pbest [n_particles_PSO];
//float k_S5_S6a_gbest   , k_S5_S6a_pbest   [n_particles_PSO];
//float k_S6_S7_gbest   , k_S6_S7_pbest   [n_particles_PSO];
//float k_S0_S11_gbest  , k_S0_S11_pbest  [n_particles_PSO];
// declare all non-changing variables
float    k_S1_S0, k_S1_S2, k_S2_S1, k_S3_S2, k_S3_S4,  k_S4_S3, k_S4_S5, k_S5_S4, k_S5_S6a,  k_S6a_S5, k_S6a_S7, k_S7_S6a, k_S5_S6,  k_S6_S5, k_S6_S7, k_S7_S6,  k_S8_S7, k_S8_S9, k_S9_S8,k_S10_S9, k_S10_S11,k_S11_S10,k_S11_S0,k_S0_S11;
float  Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc;


//-------------------------
// main body code
//------------------------
int main(int argc, char *argv[])
{
    
    long long startTime    = time(NULL);
    n_SERCA_Molecules      = 10000;         // Max number used to repeat the simulation (n_SERCA)
    max_tsteps             = 100001;     // Max number of time stepping
    dt                     = 1e-7;       // fixed time step
    n_s          	   = 12 ;         // Number of states
    n_pCa       	   = 16;         // Number of  pCa or Ca values to be simulated
    // --------------------------------------------------------------------------------------------------
    // Parameters / reference values of the transition rates (will be optimized)
    // Lower and Upper bounds on each parameter. NB: upper = 1.5* lower i.e, 50 % increase of lower value
    //---------------------------------------------------------------------------------------------------
    
    //Ca_cyt_conc_lower       = 1e-8;
    //Ca_cyt_conc_upper       = 1e-4;
    k_S0_S1_lower           = 0.1*4e7;
    k_S0_S1_upper           = 10*4e7;
    k_S2_S3_lower           = 0.1*1e8;
    k_S2_S3_upper           = 10*1e8;
    k_S7_S8_lower           = 0.1*500;
    k_S7_S8_upper           = 10*500;
    k_S9_S10_lower          = 0.1*6e2;
    k_S9_S10_upper          = 10*6e2;
//    k_S5_S6a_lower           = 06*800;
//    k_S5_S6a_upper           = 1.5*800;
//    k_S6_S7_lower           = 0.5*1;
//    k_S6_S7_upper           = 1.5*1;
//    k_S0_S11_lower          = 0.5*1.5e4;
//    k_S0_S11_upper          = 1.5*1.5e4;
    /*----------------------------*/
    /* Assign Model parameters    */
    /*----------------------------*/
    Ca_sr_conc        = 1.3e-3;// needs citation
    MgATP_conc        = 5e-3;  // needs citation
    MgADP_conc        = 36e-6; // needs citation
    Pi_conc           = 1e-3;  // needs citation
    
    //k_S0_S1           = 4e7;   // Transition rate of  E to E.Ca                       Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S1_S0           = 4.5e2;  // Transition rate of  E.Ca to E                       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S1_S2           = 120;   // Transition rate of  E.Ca to E'.Ca                   Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S2_S1           = 25;    // Transition rate of  E'.Ca to E.Ca                   Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    //k_S2_S3           = 1e8;   // Transition rate of  E'.Ca to E'.Ca2                 Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S3_S2           = 16;    // Transition rate of  E'.Ca2 to E'.Ca                 Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S3_S4           = 6e7;   // Transition rate of  E'.Ca2 to E'.ATP.Ca2            Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S4_S3           = 30;    // Transition rate of  E'.ATP.Ca2 to E'.Ca2            Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S4_S5           = 200;   // Transition rate of  E'.ATP.Ca2 to E'~P.ADP.Ca2      Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S5_S4           = 350;   // Transition rate of  E'~P.ADP.Ca2 to E'.ATP.Ca2      Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S5_S6a           = 800;   // Transition rate of  E'~P.ADP.Ca2 to *E'-P.ADP.Ca2   Units (s^-1)      Inesi Methods in Enzymo gy (1988) 157:154-190
    k_S6a_S5           = 200;   // Transition rate of *E'-P.ADP.Ca2 to E'~P.ADP.Ca2    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S6a_S7           = 500;   // Transition rate of *E'-P.ADP.Ca2 to *E'-P.Ca2       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S7_S6a           = 4e6;   // Transition rate of *E'-P.Ca2 to *E'.ADP-P.Ca2       Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S5_S6           = 6;     // Transition rate of *E'~P.ADP.Ca2 to E'~P.Ca2        Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S6_S5           = 1.25e3;// Transition rate of  E'~P.Ca2 to *E'~P.ADP.Ca2       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S6_S7           = 1;     // Transition rate of  E'~P.Ca2 to *E'-P.Ca2           Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S7_S6           = 10;    // Transition rate of *E'-P.Ca2 to E'~P.Ca2            Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    //k_S7_S8           = 500;   // Transition rate of *E'-P.Ca2 to *E-P.Ca2            Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S8_S7           = 5e5;   // Transition rate of *E-P.Ca2 to *E'-P.Ca2            Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S8_S9          = 20;    // Transition rate of *E-P.Ca2 to *E-P.Ca              Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S9_S8          = 20;    // Transition rate of *E-P.Ca to *E-P.Ca2              Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    //k_S9_S10         = 6e2;   // Transition rate of *E-P.Ca to *E-P                  Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S10_S9         = 6e4;   // Transition rate of *E-P to *E-P.Ca                  Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S10_S11         = 60;    // Transition rate of *E-P to *E-Pi                    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S11_S10         = 60;    // Transition rate of *E-Pi to *E-P                    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S11_S0          = 6e2;   // Transition rate of *E-Pi to E                       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S0_S11          = 1.5e4; // Transition rate of  E to *E-Pi                      Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    //end parameter setup
    
    
    
    
  
    srand(time(NULL)); //Random-Seed initialization (must be outside any loop)
    //------------------------------------------------------------------------------------------//
    //                                                                                          //
    //                                                                                          //
    //                                                                                          //
    //                          The Optimization Step using:                                    //
    //                                                                                          //
    //                      Particle Swarm Optimization (PSO) tools                             //
    //                                                                                          //
    //                                                                                          //
    //------------------------------------------------------------------------------------------//
  
  
    //----------------------------------------------
    // Step 1: construct particle-parameter arrays
    //          Particles positions and velocities
    //----------------------------------------------
    
    for (int i = 0; i < n_particles_PSO; i++)
    {
        //-----------
        // positions
        //-----------
        cout << " Particle " << i+1 << " initialized. " << std::endl;
        //X_Ca_cyt_conc_PSO[i] = Ca_cyt_conc_lower  + (Ca_cyt_conc_upper - Ca_cyt_conc_lower) * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	X_k_S0_S1_PSO    [i] = k_S0_S1_lower      + (k_S0_S1_upper    - k_S0_S1_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);   //original Inesi value 4e7
        X_k_S2_S3_PSO    [i] = k_S2_S3_lower      + (k_S2_S3_upper    - k_S2_S3_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);   //original Inesi value 1e8
        X_k_S7_S8_PSO    [i] = k_S7_S8_lower      + (k_S7_S8_upper    - k_S7_S8_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);   //original Inesi value 500
        X_k_S9_S10_PSO  [i] = k_S9_S10_lower    + (k_S9_S10_upper  - k_S9_S10_lower)   * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);    //original Inesi value 6e2
    
//        X_k_S5_S6a_PSO_local    [i] = k_S5_S6a_lower      + (k_S5_S6a_upper    - k_S5_S6a_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//        X_k_S6_S7_PSO_local    [i] = k_S6_S7_lower      + (k_S6_S7_upper    - k_S6_S7_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//        X_k_S0_S11_PSO_local   [i] = k_S0_S11_lower     + (k_S0_S11_upper   - k_S0_S11_lower)    * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);



        //------------------------------------------------------------------
        // Velocities : 0.25*(lower-upper)*rand : this is can be anything
        //--------------------------------------------------------------------
        //V_Ca_cyt_conc_PSO [i] = 0.25* (Ca_cyt_conc_upper - Ca_cyt_conc_lower) * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        V_k_S0_S1_PSO     [i] = 0.25* (k_S0_S1_upper    - k_S0_S1_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        V_k_S2_S3_PSO     [i] = 0.25* (k_S2_S3_upper    - k_S2_S3_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        V_k_S7_S8_PSO     [i] = 0.25* (k_S7_S8_upper    - k_S7_S8_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        V_k_S9_S10_PSO   [i] = 0.25* (k_S9_S10_upper  - k_S9_S10_lower)   * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//        V_k_S5_S6a_PSO     [i] = 0.25* (k_S5_S6a_upper    - k_S5_S6a_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//        V_k_S6_S7_PSO     [i] = 0.25* (k_S6_S7_upper    - k_S6_S7_lower)     * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//        V_k_S0_S11_PSO    [i] = 0.25* (k_S0_S11_upper   - k_S0_S11_lower)    * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    }
    
    //--------------------------------------------------------------------------
    // Step 2: solve for each particle-parameter sets to obtain residual array
    //---------------------------------------------------------------------------
  
    for (int i = 0; i < n_particles_PSO; i++)
    {
        //Ca_cyt_conc = X_Ca_cyt_conc_PSO[i];
        k_S0_S1    = X_k_S0_S1_PSO[i];
        k_S2_S3    = X_k_S2_S3_PSO[i];
        k_S7_S8    = X_k_S7_S8_PSO[i];
        k_S9_S10  = X_k_S9_S10_PSO[i];
//        k_S5_S6a    = X_k_S5_S6a_PSO_local[i];
//        k_S6_S7    = X_k_S6_S7_PSO_local[i];
//        k_S0_S11   = X_k_S0_S11_PSO_local[i];
        //-----------------------------------------------------------------------------------------
        // Call the force_pCa_curve function to get the force as a function of Ca++ concentrations.
        // NB: this function implicitly calls the other functions.
        //-----------------------------------------------------------------------------------------
        
        residual_cost_func[i] = get_Residual  (n_SERCA_Molecules,
                                               max_tsteps,
                                               dt,
                                               n_s,
                                               n_pCa,
                                               k_S0_S1,
                                               k_S2_S3,
                                               k_S7_S8,
                                               k_S9_S10, k_S1_S0, k_S1_S2,  k_S2_S1, k_S3_S2, k_S3_S4,  k_S4_S3, k_S4_S5, k_S5_S4, k_S5_S6a,  k_S6a_S5, k_S6a_S7, k_S7_S6a, k_S5_S6,  k_S6_S5, k_S6_S7, k_S7_S6,  k_S8_S7, k_S8_S9, k_S9_S8,k_S10_S9, k_S10_S11,k_S11_S10,k_S11_S0,k_S0_S11, Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc
                                               );
      
        cout << " Particle # " << i+1 << endl;
        cout << " k_S0_S1 = " << k_S0_S1 << ", k_S2_S3 = " << k_S2_S3 << ", k_S7_S8 = " << k_S7_S8 << ", k_S9_S10 = " << k_S9_S10 << endl;
        cout << " Residual =  " << residual_cost_func[i] << endl;  //**

    } // close loop of particle
     cout << "One iteration runtime: " << (time(NULL)-startTime) << " second(s)" << std::endl;
  
  
    //------------------------------------------
    // Find Min of Residual (i.e., global best)
    //------------------------------------------
    float Res_gbest = residual_cost_func [0];
    int i_Res_gbest = 0;
    for (int i = 0; i < n_particles_PSO; i++)
    {
        if (residual_cost_func [i] < Res_gbest)
        {
            Res_gbest    = residual_cost_func [i];
            i_Res_gbest  = i;

        }
        
        Res_pbest[i] = residual_cost_func [i];  // Residual personal best
    }
    
    //--------------------------------------------------------
    // obtain the parameters that give global (gbest)
    //--------------------------------------------------------
    //Ca_cyt_conc_gbest = X_Ca_cyt_conc_PSO[i_Res_gbest];
    k_S0_S1_gbest     = X_k_S0_S1_PSO[i_Res_gbest];
    k_S2_S3_gbest     = X_k_S2_S3_PSO[i_Res_gbest];
    k_S7_S8_gbest     = X_k_S7_S8_PSO[i_Res_gbest];
    k_S9_S10_gbest   = X_k_S9_S10_PSO[i_Res_gbest];
//    k_S6_S7_gbest     = X_k_S6_S7_PSO[i_Res_gbest];
//    k_S0_S11_gbest    = X_k_S0_S11_PSO[i_Res_gbest];
    //---------------------------------------------------------------------------------------
    // obtain the parameters that give personal best (pbest)
    // Note: this can be combined with one of the other loop but keep it like that for now
    //----------------------------------------------------------------------------------------
    
    for (int i = 0; i < n_particles_PSO; i++)
    {
        //Ca_cyt_conc_pbest[i] = X_Ca_cyt_conc_PSO[i];
        k_S0_S1_pbest[i]     = X_k_S0_S1_PSO[i];
        k_S2_S3_pbest[i]     = X_k_S2_S3_PSO[i];
        k_S7_S8_pbest[i]     = X_k_S7_S8_PSO[i];
        k_S9_S10_pbest[i]   = X_k_S9_S10_PSO[i];
//        k_S5_S6a_pbest[i]     = X_k_S5_S6a_PSO[i];
//        k_S6_S7_pbest[i]     = X_k_S6_S7_PSO[i];
//        k_S0_S11_pbest[i]    = X_k_S0_S11_PSO[i];
    }
    
    //----------------------------------------------------------------------------------
    //
    //
    //                 Swarm Iteration Step over all particles
    //
    //
    //------------------ --------------------------------------------------------------
    const int max_iter = 100;
    float total_gbest[max_iter]; ; 
if (max_iter != 0)
{     
float w_max, w_min, dw, w;
    float c1, c2;
    w_max = 1.0;
    w_min = 0.3;
    dw = (w_max-w_min)/max_iter;
    c1 = 1.05;
    c2 = 1.05;
    for (int it = 0; it < max_iter+1; it++)
    { // begin swarm iteration
        w = w_min +it*dw;
        for (int i = 0; i < n_particles_PSO; i++)
        { // begin loop over all particles
            
            //-----------------
            // Velocity update
            //-----------------
            
            //V_Ca_cyt_conc_PSO[i] = w  * V_Ca_cyt_conc_PSO [i] + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (Ca_cyt_conc_pbest[i]   - X_Ca_cyt_conc_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (Ca_cyt_conc_gbest - X_Ca_cyt_conc_PSO[i]) ;
            
            
            V_k_S0_S1_PSO[i]    = w  * V_k_S0_S1_PSO [i]      + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S0_S1_pbest[i]   - X_k_S0_S1_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S0_S1_gbest      - X_k_S0_S1_PSO[i]) ;
            
            V_k_S2_S3_PSO[i]    = w  * V_k_S2_S3_PSO [i]      + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S2_S3_pbest[i]   - X_k_S2_S3_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S2_S3_gbest      - X_k_S2_S3_PSO[i]) ;
            
            V_k_S7_S8_PSO[i]    = w  * V_k_S7_S8_PSO [i]      + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S7_S8_pbest[i]   - X_k_S7_S8_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S7_S8_gbest      - X_k_S7_S8_PSO[i]) ;
            
            V_k_S9_S10_PSO[i]  = w  * V_k_S9_S10_PSO [i]    + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S9_S10_pbest[i]   - X_k_S9_S10_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S9_S10_gbest      - X_k_S9_S10_PSO[i]) ;

//            V_k_S5_S6a_PSO[i]  = w  * V_k_S5_S6a_PSO [i]    + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S5_S6a_pbest[i]   - X_k_S5_S6a_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S5_S6a_gbest      - X_k_S5_S6a_PSO[i]) ;

//            V_k_S6_S7_PSO[i]  = w  * V_k_S6_S7_PSO [i]    + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S6_S7_pbest[i]   - X_k_S6_S7_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S6_S7_gbest      - X_k_S6_S7_PSO[i]) ;

//            V_k_S0_S11_PSO[i]  = w  * V_k_S0_S11_PSO [i]    + c1 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S0_S11_pbest[i]   - X_k_S0_S11_PSO[i]) + c2 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * (k_S0_S11_gbest      - X_k_S0_S11_PSO[i]) ;
            
            //-----------------
            // position update
            //-----------------
            //X_Ca_cyt_conc_PSO[i] = X_Ca_cyt_conc_PSO[i] + V_Ca_cyt_conc_PSO[i];
            X_k_S0_S1_PSO[i]     = X_k_S0_S1_PSO[i]    + V_k_S0_S1_PSO[i];
            X_k_S2_S3_PSO[i]     = X_k_S2_S3_PSO[i]     + V_k_S2_S3_PSO[i];
            X_k_S7_S8_PSO[i]     = X_k_S7_S8_PSO[i]     + V_k_S7_S8_PSO[i];
            X_k_S9_S10_PSO[i]   = X_k_S9_S10_PSO[i]   + V_k_S9_S10_PSO[i];
            //X_k_S5_S6a_PSO_local[i]     = X_k_S5_S6a_PSO[i]     + V_k_S5_S6a_PSO[i];
 			      //X_k_S6_S7_PSO_local[i]     = X_k_S6_S7_PSO[i]     + V_k_S6_S7_PSO[i];
 			      //X_k_S0_S11_PSO_local[i]    = X_k_S0_S11_PSO[i]    + V_k_S0_S11_PSO[i];
            //---------------------------
            // model parameter update
            //---------------------------
            //Ca_cyt_conc = X_Ca_cyt_conc_PSO[i];
            k_S0_S1     = X_k_S0_S1_PSO[i];
            k_S2_S3     = X_k_S2_S3_PSO[i];
            k_S7_S8     = X_k_S7_S8_PSO[i];
            k_S9_S10   = X_k_S9_S10_PSO[i];
            //k_S5_S6a     = X_k_S5_S6a_PSO_local[i];
            //k_S6_S7     = X_k_S6_S7_PSO_local[i];
            //k_S0_S11    = X_k_S0_S11_PSO_local[i]; 
            //-----------------------------------------------------
            // residual update using the new particles/parameters
            //----------------------------------------------------
            residual_cost_func[i] = get_Residual  (n_SERCA_Molecules,
                                                   max_tsteps,
                                                   dt,
                                                   n_s,
                                                   n_pCa,
                                                   k_S0_S1,
                                                   k_S2_S3,
                                                   k_S7_S8,
                                                   k_S9_S10, k_S1_S0, k_S1_S2,  k_S2_S1, k_S3_S2, k_S3_S4, k_S4_S3, k_S4_S5, k_S5_S4, k_S5_S6a,  k_S6a_S5, k_S6a_S7, k_S7_S6a, k_S5_S6,  k_S6_S5, k_S6_S7, k_S7_S6,  k_S8_S7, k_S8_S9, k_S9_S8,k_S10_S9, k_S10_S11,k_S11_S10,k_S11_S0,k_S0_S11, Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc
                                                   );

        cout << " Particle # " << i+1 << endl;
        cout << "        k_S0_S1 = " << k_S0_S1 << ", k_S2_S3 = " << k_S2_S3 << ", k_S7_S8 = " << k_S7_S8 << ", k_S9_S10 = " << k_S9_S10 << endl;
        cout << " Residual =  " << residual_cost_func[i];            
        }// end looping over all particles to have new Residual vector
        
      
        //--------------------------------------------------
        // Find New Min of Residual (i.e., new global best)
        //--------------------------------------------------
        float min_Res = residual_cost_func [0];
        int i_min_Res = 0;
        for (int i = 0; i < n_particles_PSO; i++)
        {
            if (residual_cost_func [i] < min_Res)
            {
                min_Res    = residual_cost_func [i];
                i_min_Res  = i;   
                 
	    }
	    
            cout << " New Residuals         = " << residual_cost_func [i] << endl;
        }
  
        
        //-------------------------------
        // Check for update min residual
        //-------------------------------
        if (min_Res <= Res_gbest)
        {
            Res_gbest = min_Res;
            //-------------------------------------------------
            // obtain the parameters that give global (gbest)
            //--------------------------------------------------
            //Ca_cyt_conc_gbest = X_Ca_cyt_conc_PSO[i_min_Res];
            k_S0_S1_gbest     = X_k_S0_S1_PSO[i_min_Res];
            k_S2_S3_gbest     = X_k_S2_S3_PSO[i_min_Res];
            k_S7_S8_gbest     = X_k_S7_S8_PSO[i_min_Res];
            k_S9_S10_gbest   = X_k_S9_S10_PSO[i_min_Res];
//            k_S5_S6a_gbest     = X_k_S5_S6a_PSO[i_min_Res];
//            k_S6_S7_gbest     = X_k_S6_S7_PSO[i_min_Res];
//            k_S0_S11_gbest    = X_k_S0_S11_PSO[i_min_Res];
        }

        total_gbest [it] = Res_gbest;
	cout << " " << endl;
	cout << "The total global best is now : " << total_gbest [it]<< endl;
	cout << " " << endl;
        

       // write iteration vs. gbest in a file
       std::string outfilename = "iterations_vs_global_best.csv";
    
       ofstream iter_out(outfilename.c_str()); //opening an output stream for file test.txt
    		if(iter_out.is_open()) //checking whether file could be opened or not.
    		{
        		// create headers for file
        		iter_out << "iteration  vs  gbest" << endl; // write the average force
        		for (int g = 0; g < max_iter; g++)  // time marching
        		{
                        iter_out << g << "  " << total_gbest[g]  << endl; // write the average force
        		}
        		cout << "Iterations and Global best successfully saved into the file " << outfilename << endl;
        	}


        for (int i = 0; i < n_particles_PSO; i++)
        {
            if(residual_cost_func [i]  <= Res_pbest[i] )
            {
                //Ca_cyt_conc_pbest[i] = X_Ca_cyt_conc_PSO[i];
                k_S0_S1_pbest[i]     = X_k_S0_S1_PSO[i];
                k_S2_S3_pbest[i]     = X_k_S2_S3_PSO[i];
                k_S7_S8_pbest[i]     = X_k_S7_S8_PSO[i];
                k_S9_S10_pbest[i]   = X_k_S9_S10_PSO[i];
//                k_S5_S6a_pbest[i]     = X_k_S5_S6a_PSO[i];
//                k_S6_S7_pbest[i]     = X_k_S6_S7_PSO[i];
//                k_S0_S11_pbest[i]    = X_k_S0_S11_PSO[i];
                
                //--------
                Res_pbest[i] = residual_cost_func [i];
                
            }
        }
        
        
    }}// end swarm iteration
  
    cout << "\"Res_gbest\","       << Res_gbest  <<  endl;
    cout << "\"k_S0_S1_gbest  \","  << k_S0_S1_gbest   << " original Inesi value 4e+07" << endl;
    cout << "\"k_S2_S3_gbest  \","  << k_S2_S3_gbest   << " original Inesi value 1e+08" << endl;
    cout << "\"k_S8_S9_gbest \","  << k_S7_S8_gbest   << " original Inesi value 500"   << endl;
    cout << "\"k_S9_S10_gbest\","  << k_S9_S10_gbest << " original Inesi value 600"   << endl;
//    cout << "\"k_S5_S6a_pbest\","   << k_S5_S6a_pbest << endl;
//    cout << "\"k_S6_S7_pbest\","   << k_S6_S7_pbest << endl;
//    cout << "\"k_S0_S11_pbest\","  << k_S0_S11_pbest << endl;
    cout << "Total Optimization Runtime: " << (time(NULL)-startTime) << " second(s)" << endl;
    


    lastRun(n_SERCA_Molecules,
    		max_tsteps,
    		dt,
    		n_s,
    		n_pCa,
    		k_S0_S1_gbest,
    		k_S2_S3_gbest,
    		k_S7_S8_gbest,
    		k_S9_S10_gbest, k_S1_S0, k_S1_S2,  k_S2_S1, k_S3_S2, k_S3_S4,  k_S4_S3, k_S4_S5, k_S5_S4, k_S5_S6a, k_S6a_S5, k_S6a_S7, k_S7_S6a, k_S5_S6,  k_S6_S5, k_S6_S7, k_S7_S6,  k_S8_S7, k_S8_S9, k_S9_S8,k_S10_S9, k_S10_S11, k_S11_S10, k_S11_S0, k_S0_S11, Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc);
    
    return 0;
} // end main function
