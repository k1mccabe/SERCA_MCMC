
//-----------------------------------------------------------------------------------------------------
//                     University of California, San Diego
//                      Dept. of Chemistry & Biochemistry
//-----------------------------------------------------------------------------------------------------
// Authors: Sophia P. Hirakis & Kimberly J. McCabe
// Year  :  4/2018
//-----------------------------------------------------------------------------------------------------
// This code uses the "Markov Chain Monte Carlo" algorithm to solve a model of the SERCA pump.
//
//
// The SERCA MCMC model is based on the published model by Inesi (1988)
//
//          [S1]       [S2]             [S3]                 [S4]                 [S5]
//          E.Ca <==> E'.Ca  + Ca <==> E'.Ca2 (+ ATP) <==> E'.ATP.Ca2  <==>   E'~P.ADP.Ca2
//           /\                                                                  //  \\
//           ||                                                                 //    \\
//           ||                                                         [S6]   //      \\  [S8]
//    + Ca   ||                                                      *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
//           ||                                                                \\      //
//           ||                                                        (+ ADP)  \\    //
//           \/                                                                  \\  //
//           E(+Pi)<==> *E-Pi<==> *E-P + Ca<==> *E-P.Ca + Ca<==> *E-P.Ca2 <==>*E'-P.Ca2
//         [S0]     [S12]      [S11]           [S10]              [S9]            [S7]
//
//
// State Reaction               State Product          Rate(f)   Rate(r)
//  S0   E + Ca             <==> S1   E.Ca             k.S0.S1   k.S1.S0
//  S1   E.Ca               <==> S2   E'.Ca            k.S1.S2   k.S2.S1
//  S2   E'.Ca + Ca         <==> S3   E'.Ca2           k.S2.S3   k.S3.S2
//  S3   E'.Ca2 (+ ATP)     <==> S4   E'.ATP.Ca2       k.S3.S4   k.S4.S3
//  S4   E'.ATP.Ca2         <==> S5   E'~P.ADP.Ca2     k.S4.S5   k.S5.S4
//  S5   E'~P.ADP.Ca2       <==> S6  *E'-P.ADP.Ca2     k.S5.S6   k.S6.S5
//  S6  *E'-P.ADP.Ca2       <==> S7  *E'-P.Ca2 (+ ADP) k.S6.S7   k.S7.S6
//  S5   E'~P.ADP.Ca2       <==> S8   E'~P.Ca2 (+ ADP) k.S5.S8   k.S8.S5
//  S8   E'~P.Ca2           <==> S7  *E'-P.Ca2 (+ ADP) k.S8.S7   k.S7.S8
//  S7  *E'-P.Ca2 (+ ADP)   <==> S9  *E-P.Ca2          k.S7.S9   k.S9.S7
//  S9  *E-P.Ca2            <==> S10 *E-P.Ca + Ca      k.S9.S10  k.S10.S9
//  S10 *E-P.Ca + Ca        <==> S11 *E-P + Ca         k.S10.S11 k.S11.S10
//  S11 *E-P + Ca           <==> S12 *E-Pi             k.S11.12  k.S12.S11
//  S12 *E-Pi               <==> S0   E + (Pi)         k.S12.S0  k.S0.S12


//Libraries defined here
//--------------------------------------------------------------------------
#include <iostream>  // library that contains file input/output functions
#include <fstream>   // library that contains file input/output functions
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include "lastRun.h"

using namespace std;

float last_count_S0, last_count_S1, last_count_S2, last_count_S3, last_count_S4, last_count_S5, last_count_S6, last_count_S7, last_count_S8, last_count_S9, last_count_S10, last_count_S11, last_count_S12;
float S0_SS_last, S1_SS_last, S2_SS_last, S3_SS_last, S4_SS_last, S5_SS_last, S6_SS_last, S7_SS_last, S8_SS_last, S9_SS_last, S10_SS_last, S11_SS_last, S12_SS_last;
float S0_temp_last, S1_temp_last, S2_temp_last, S3_temp_last, S4_temp_last, S5_temp_last, S6_temp_last, S7_temp_last, S8_temp_last, S9_temp_last, S10_temp_last, S11_temp_last, S12_temp_last;
int state_last;
bool open_closed_last; //O is open 1 is closed
float Ca_cyt_conc_last;


//--------------------------------------------------------------------------//


void lastRun  		(int    & n_SERCA_Molecules,
                     int    & tsteps,
                     float  & dt,
                     int    & n_s,
                     int    & n_pCa,
                     float  & k_S0_S1,
                     float  & k_S2_S3,
                     float  & k_S7_S9,
                     float  & k_S10_S11,float  & k_S1_S0, float  & k_S1_S2,  float  & k_S2_S1, float  & k_S3_S2, float  & k_S3_S4,  float  & k_S4_S3, float  & k_S4_S5, float  & k_S5_S4, float  & k_S5_S6,  float  & k_S6_S5, float  & k_S6_S7, float  & k_S7_S6, float  & k_S5_S8,  float  & k_S8_S5, float  & k_S8_S7, float  & k_S7_S8,  float  & k_S9_S7, float  & k_S9_S10, float  & k_S10_S9,float  & k_S11_S10, float  & k_S11_S12,float  & k_S12_S11,float  & k_S12_S0,float  & k_S0_S12, float  & Ca_sr_conc,float  & MgATP_conc,float  & MgADP_conc,float  & Pi_conc
                     )

{
   
   int save_jump = 100; //how many output values should we keep? To minimize memory usage, we will keep every 10 timepoints.
    
    float calConc_Exp[16] = {   1.13465021562703E-07,
        1.48013728928924E-07,
        1.87545047401295E-07,
        2.37746427649773E-07,
        2.86177839072689E-07,
        3.34581558654772E-07,
        3.82579504194903E-07,
        4.40880103529033E-07,
        5.15498018194447E-07,
        6.0268752205741E-07,
        7.04360511231999E-07,
        8.41433890215616E-07,
        9.8310521528177E-07,
        1.209326027507E-06,
        1.46539261994034E-06,
        1.92506766806173E-06};
    float boundSS_max_temp2 = 0; // this will figure out the highest bound Ca for our loop
    float ss_bound_Ca2[n_pCa];
    float norm_ss_bound_Ca2[n_pCa];
    
    
    for (int i = 0; i < n_pCa; i++)
    {
        Ca_cyt_conc_last       = calConc_Exp[i];  // needs citation
    
    //-----------------------
    // SIMULATION FOR SS_last CURVE
    //-----------------------
    //Creating a vector named tstapes (list of elements of States S0-S12) defining the occupancy of each state according to time at each time step
    //CHANGED - this was taking up too much memory. Now only saving every 10 timesteps.
    float S0[(tsteps-1)/10];         // used to find how many S0 molecules in each iteration ("E")
    float S1[(tsteps-1)/10];         //  "E.Ca"
    float S2[(tsteps-1)/10];         //  "E'.Ca"
    float S3[(tsteps-1)/10];         //  "E'.Ca2"
    float S4[(tsteps-1)/10];         //  "E'.ATP.Ca2"
    float S5[(tsteps-1)/10];         //  "E'~P.ADP.Ca2"
    float S6[(tsteps-1)/10];         // "*E'~P.ADP.Ca2"
    float S7[(tsteps-1)/10];         // "*E'-P.Ca2"
    float S8[(tsteps-1)/10];         //  "E'~P.Ca2"
    float S9[(tsteps-1)/10];         // "*E-P.Ca2"
    float S10[(tsteps-1)/10];        // "*E-P.Ca"
    float S11[(tsteps-1)/10];        // "*E-P"
    float S12[(tsteps-1)/10];        // "*E-Pi"

    // create a temp_lastorary last_count (initial last_count of state and set it equal to zero. The last_count will be calculated and set equal to the variables above for each state. This is because C++ will use the last value that was defined, so its a “clear” or reset function
    S0_temp_last = 0;
    S1_temp_last = 0;
    S2_temp_last = 0;
    S3_temp_last = 0;
    S4_temp_last = 0;
    S5_temp_last = 0;
    S6_temp_last = 0;
    S7_temp_last = 0;
    S8_temp_last = 0;
    S9_temp_last = 0;
    S10_temp_last = 0;
    S11_temp_last = 0;
    S12_temp_last = 0;
    
    S0_SS_last = 0;
    S1_SS_last = 0;
    S2_SS_last = 0;
    S3_SS_last = 0;
    S4_SS_last = 0;
    S5_SS_last = 0;
    S6_SS_last = 0;
    S7_SS_last = 0;
    S8_SS_last = 0;
    S9_SS_last = 0;
    S10_SS_last = 0;
    S11_SS_last = 0;
    S12_SS_last = 0;

        
     //-----------------------------------------------------------------------------------------------------------------//
    // start repeat loop i.e., using r-index - REPEAT PROCESS_last FOR EACH SERCA MOLECULE
    //------------------------------------------------------------------------------------------------------------------//
    //start in state 0 which is SERCA in open, then go through in time and ##
    for (int rr = 0; rr < n_SERCA_Molecules; rr++) //repetition of whole simulation to smooth curve
    {
        state_last = 0; // Note: each time we repeat the simulation, we should set all SERCA to state 0 (EiH2)
        open_closed_last = 0; // initially, each molecule is in closed conformation (0)
        //---------------------------------------------------------------------------------------------------------------//
        // start time loop i.e., using n-index
        //----------------------------------------------------------------------------------------------------------------//
        int output_last_count = 10; //every 10 timepoints, we will save the data (starting with data point 0)
        for (int n = 0; n < tsteps; n++)  // time marching
        {  // begin n-loop for time marching
            //
            //setting last_counts of all states equal to zero to initialize which will later be used to find how many SERCA in each iteration
            //
            last_count_S0  = 0.0;   // "E"
            last_count_S1  = 0.0;  //  "E.Ca"
            last_count_S2  = 0.0;  //  "E'.Ca"
            last_count_S3  = 0.0;  //  "E'.Ca2"
            last_count_S4  = 0.0;  //  "E'.ATP.Ca2"
            last_count_S5  = 0.0;  //  "E'~P.ADP.Ca2"
            last_count_S6  = 0.0;  // "*E'~P.ADP.Ca2"
            last_count_S7  = 0.0;  // "*E'-P.Ca2"
            last_count_S8  = 0.0;  //  "E'~P.Ca2"
            last_count_S9  = 0.0;  // "*E-P.Ca2"
            last_count_S10 = 0.0;  // "*E-P.Ca"
            last_count_S11 = 0.0;  // "*E-P"
            last_count_S12 = 0.0;  // "*E-Pi"
            
            
            //last_count << "open_closed_last = " << open_closed_last << endl;
            
            //update all RUs for current timestep
            
            int current_state;
            float p1, p2, p3; //Probabilities of each state
            float randNum = 0;
            //
            //-------------------------------------------------------------------------------------------------------------------------
            current_state = state_last;     // get current state ID number (e.g. 0 = E , etc)
            //
            //setting index equal from 0 to 12 for each state
            //
            //
            // generate random numbers
            randNum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); // Generate a random number between 0 and 1
            //
            //--------------------------------------------------------------------------------__________________________________________________________________________________________________
            if(current_state == 0)
            {
                p1 =       (k_S0_S1  * Ca_cyt_conc_last * dt); // (pseudo first-order) bimolecular  forward transition to  E.Ca          [S1]
                p2 = p1 +  (k_S0_S12 * Pi_conc     * dt); // (pseudo first-order) bimolecular backward transition to *E-Pi          [S12]
                if(randNum < p1)  //If the random number is leSS_last than first probability, then transition to the next state
                {
                    state_last = 1;   // forward transition to  E.Ca [S1]
                }
                else if(randNum < p2) //If the random number is leSS_last than second probability , then transition to the previous state
                {
                    state_last = 12; // backward transition to *E-Pi [S12]
                }
            } //if it is not greater than either probability then stay in the same state
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 1): Then     S0 [E]  <-- S1 {E.Ca} --> S2 [E'.Ca]                  else stay as E.Ca
            //          [S1]      [S2]
            //          E.Ca <==> E'.Ca  + Ca
            //           /\
            //           ||
            //           \/
            //    (Pi +) E
            //          [S0]
            //-----------------------------------------------------------------------------------------------------------------------
            if(current_state == 1)
            {
                p1 =       (k_S1_S2 * dt);           // (first order)        unimolecular forward transition forward E'.Ca          [S2]
                p2 = p1 +  (k_S1_S0 * Pi_conc * dt); // (pseudo first-order) bimolecular backward transition back to E              [S0]
                if(randNum < p1)  //f your random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 2;   // forward transition to  E'.Ca    [S2]
                }
                else if(randNum < p2) //If the random number is leSS_last than second probability , then transition to the previous state
                {
                    state_last = 0;
                }
            }
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 2): Then     S1 [E.Ca]  <-- S2 {E'.Ca} +Ca --> S3 [E'.Ca2]         else stay as E.Ca
            //
            //          [S1]       [S2]             [S3]
            //          E.Ca <==> E'.Ca  + Ca <==> E'.Ca2
            //
            //-----------------------------------------------------------------------------------------------------------------------
            else if(current_state == 2)
            {
                p1 =       (k_S2_S3 * dt);           // (first order)        unimolecular  forward transition to E'.Ca2             [S3]
                p2 = p1 +  (k_S2_S1 * Pi_conc * dt); // (first-order)        unimolecular backward transition to E.Ca               [S0]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state              [S3]
                {
                    state_last = 3;   // forward transition to E'.Ca2    [S3]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state [S0]
                {
                    state_last = 1; // backward transition to  E'.Ca     [S1]
                }
            } //if it is not greater than either probability  then stay in the same state
            //
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 3): Then     S2 [E'.Ca]  <-- S3 {E'.Ca2} --> S4 [E'.ATP.Ca2]       else stay as E'.Ca2
            //
            //
            //           [S2]              [S3]                  [S4]
            //           E'.Ca  + Ca <==> E'.Ca2 (+ ATP) <==> E'.ATP.Ca2
            //
            //
            //-----------------------------------------------------------------------------------------------------------------------
            else if(current_state == 3)
            {
                p1 =       (k_S3_S4 * MgATP_conc * dt);     // (pseudo first-order) bimolecular   forward transition to E'.ATP.Ca2 [S4]
                p2 = p1 +  (k_S3_S2 * dt);                  // (first order)        unimolecular backward transition to E'.Ca      [S2]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 4;   // forward transition to E'.ATP.Ca2 [S4]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 2;  // backward transition to E'.Ca      [S2]
                }
            } //if it is not greater than either probability then stay in the same state
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 4): Then     S3 [E'.Ca2] <-- S4 {E'.ATP.Ca2} --> S5 [E'~P.ADP.Ca2]    else stay as E'.ATP.Ca2
            //
            //
            //   [S3]                 [S4]                 [S5]
            //  E'.Ca2 (+ ATP) <==> E'.ATP.Ca2  <==>   E'~P.ADP.Ca2
            //
            //-----------------------------------------------------------------------------------------------------------------------
            else if(current_state == 4)
            {
                p1 =       (k_S4_S5 * MgATP_conc * dt);     // (pseudo first-order) bimolecular   forward transition to E'~P.ADP.Ca2 [S5]
                p2 = p1 +  (k_S4_S3              * dt);     // (first order)        unimolecular backward transition to E'.Ca2       [S3]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 5;   // forward transition to E'~P.ADP.Ca2 [S5]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 3;  // backward transition to E'.Ca2       [S3]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 5):  Then    S4 [E'.ATP.Ca2]  <--  S5 {E'~P.ADP.Ca2}  --> S6 [*E'-P.ADP.Ca2] else stay as E'~P.ADP.Ca2
            //                                                        ||
            //                                                        \/
            //                                                 S8  [E'~P.Ca2]
            //    [S4]           [S5]
            // E'.ATP.Ca2 <==> E'~P.ADP.Ca2
            //                  //  \\
            //                 //    \\
            //         [S6]   \/      \/ [S8]
            //       *E'-P.ADP.Ca2   E'~P.Ca2 (+ ADP)
            //
            //
            //--------------------------------------------------------------------------------------------------------------------
            else if(current_state == 5)
            {
                p1 =        (k_S5_S6 * dt); // (first order)        unimolecular  forward transition to *E'-P.ADP.Ca2               [S6]
                p2 = p1 +   (k_S5_S4 * dt); // (first order)        unimolecular backward transition to  E'.ATP.Ca2                 [S4]
                p3 = p2 +   (k_S5_S8 * dt); // (first order)        unimolecular  forward transition to  E'~P.Ca2                   [S8]
                if(randNum < p1)
                {
                    state_last = 6;   //forward transition to *E'-P.ADP.Ca2 [S6]
                }
                else if(randNum < p2)
                {
                    state_last = 8;  //forward transition to  E'~P.Ca2      [S8]
                }
                else if(randNum < p3)
                {
                    state_last = 4; //reverse transition to   E'.ATP.Ca2    [S4]
                }
            }
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 6):  Then        S5 [E'~P.ADP.Ca2]  <-- S6 {*E'-P.ADP.Ca2} --> S7 [*E'-P.Ca2]       else stay as *E'-P.ADP.Ca2
            //
            //                     [S5]
            //                  E'~P.ADP.Ca2
            //                  //
            //                 //
            //         [S6]   \/
            //      *E'-P.ADP.Ca2
            //                /\
            //                 \\
            //                  \\
            //                *E'-P.Ca2
            //                   [S7]
            //
            //----------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 6)
            {
                p1 =       (k_S6_S7 * dt);     // (pseudo first-order) bimolecular   forward transition to *E'-P.Ca2                 [S7]
                p2 = p1 +  (k_S6_S5 * dt);     // (first order)        unimolecular backward transition to  E'~P.ADP.Ca2             [S5]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 7;   // forward transition to *E'-P.Ca2    [S7]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 5;  // backward transition to E'~P.ADP.Ca2 [S5]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state = 7) Then          S6 [*E'-P.ADP.Ca2]  <-- S7 {(+ ADP) *E'-P.Ca2} --> S8 [E'~P.Ca2]  else stay as *E'-P.Ca2
            //
            //               [S6]              [S8]
            //            *E'-P.ADP.Ca2      E'~P.Ca2
            //                      \\      //
            //                       \\    //
            //                        \\  //
            // *E-P.Ca2 <==> (+ ADP) *E'-P.Ca2
            //    [S9]                 [S7]
            //
            //
            //-------------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 7)
            {
                p1 =      (k_S7_S9              * dt); // (first order)     unimolecular  forward transition to *E-P.Ca2            [S9]
                p2 = p1 + (k_S7_S6 * MgADP_conc * dt); // (first order)     unimolecular backward transition to *E'-P.ADP.Ca2       [S6]
                p3 = p2 + (k_S7_S8              * dt); // (first order)     unimolecular backward transition to  E'~P.Ca2           [S8]
                if(randNum < p1)
                {
                    state_last = 9; //forward transition to *E-P.Ca2         [S9]
                }
                else if(randNum < p2)
                {
                    state_last = 6; //reverse transition to *E'-P.ADP.Ca2    [S6]
                }
                else if(randNum < p3)
                {
                    state_last = 8; //reverse transition to  E'~P.Ca2        [S8]
                }
            }
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state = 8) Then          S5 [E'~P.ADP.Ca2] <-- S8 {(ADP +) E'~P.Ca2 }  --> S7 [*E'-P.Ca2]
            //
            //        [S5]
            //      E'~P.ADP.Ca2
            //          \\
            //           \\
            //           [S8]
            //        E'~P.Ca2 (+ADP)
            //            //
            //           //
            //          //
            //     *E'-P.Ca2 (+ ADP)
            //        [S7]
            //
            //
            //------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 8)
            {
                p1 =       (k_S8_S7  * dt);     // (first-order) unimolecular  forward transition to *E'-P.Ca2                      [S7]
                p2 = p1 +  (k_S8_S5  * dt);     // (first order) unimolecular backward transition to E'~P.ADP.Ca2                   [S5]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 7;   // forward transition to *E'-P.Ca2      [S7]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 5;  // backward transition to  E'~P.ADP.Ca2  [S5]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 9): Then      S7 [*E'-P.Ca2] <-- S9 {*E-P.Ca2} --> S10 [*E-P.Ca] +Ca   else stay as *E-P.Ca2
            //
            //*E'-P.Ca2 <==  *E-P.Ca2  ==>  *E-P.Ca + Ca
            //  [S7]            [S9]          [S10]
            //
            //------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 9)
            {
                p1 =       (k_S9_S10 * MgATP_conc * dt);     // (pseudo first-order) unimolecular  forward transition to *E-P.Ca    [S10]
                p2 = p1 +  (k_S9_S7               * dt);     // (first order)        unimolecular backward transition to *E'-P.Ca2  [S7]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 10;   // forward transition to *E-P.Ca      [S10]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 7;   // backward transition to *E'-P.Ca2    [S7]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 10): Then   S9 [*E-P.Ca2] <-- S10 {*E-P.Ca} --> S11 [*E-P] +Ca   else stay as *E-P.Ca
            //
            // *E-P.Ca2 <==  Ca + *E-P.Ca   ==> *E-P + Ca
            //   [S9]             [S10]         [S11]
            //
            //-------------------------------------------------------------------------------------------------------------------------
            
            else if(current_state == 10)
            {
                p1 =       (k_S10_S11             * dt);      //(first-order)       unimolecular  forward transition to *E-P        [S11]
                p2 =       (k_S10_S9 * Ca_sr_conc * dt);     // (pseudo first-order) bimolecular backward transition to *E-P.Ca2    [S9]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 11;   // forward transition to *E-P         [S11]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 9; // backward transition to   *E-P.Ca2     [S9]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 11): Then   S10 [*E-P.Ca] <-- S11 {(Ca +) *E-P}--> S12 [*E-Pi]   else stay as *E-P.Ca
            //
            // *E-P.Ca2 <==  Ca + *E-P.Ca    ==> *E-P + Ca
            //   [S9]              [S10]         [S11]
            //
            //-------------------------------------------------------------------------------------------------------------------------
            
            else if(current_state ==11)
            {
                p1 =       (k_S11_S12              * dt);    // (first-order)       unimolecular  forward transition to *E-Pi       [S12]
                p2 =       (k_S11_S10 * Ca_sr_conc * dt);    // (pseudo first-order) bimolecular backward transition to *E-P.Ca     [S10]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 12;   // forward transition to *E-Pi        [S12]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 10; // backward transition to  *E-P.Ca      [S10]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 12): Then   S11 [(Ca +) *E-P] <-- S12 {*E-Pi} --> S0 [E]  else stay as *E-Pi
            //
            //   *E-P  <==  *E-Pi  ==>  E (+ Pi )
            //   [S11]      [S12]     [S0]
            //
            //------------------------------------------------------------------------------------------------------------------------
            else if(current_state ==12)
            {
                p1 =       (k_S12_S0  * dt);    // (first-order)       unimolecular  forward transition to  E                       [S0]
                p2 =       (k_S12_S11 * dt);    // (pseudo first-order) bimolecular backward transition to *E-P                     [S11]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state_last = 0;   // forward transition to  E            [S0]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state_last = 11; // backward transition to *E-P          [S11]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            // ----------------------------------------------------------------------------------------------------------------------
            //                              END UPDATE STEP @ EACH TIME
            //-----------------------------------------------------------------------------------------------------------------------
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            
            //---------------------------------------------------------------
            // Obtaining Force estimate based on the Markov state @ each time (every 10 timesteps)
            //---------------------------------------------------------------
            
            if (output_last_count == save_jump)
            {
                if      (state_last == 0)
                {
                    S0[n/save_jump]   = S0[n/save_jump] + 1.0;
                }
                else if(state_last ==1)
                {
                    S1[n/save_jump]   = S1[n/save_jump] + 1.0;
                }
                else if(state_last ==2)
                {
                    S2[n/save_jump]   = S2[n/save_jump] + 1.0;
                }
                else if(state_last ==3)
                {
                    S3[n/save_jump]   = S3[n/save_jump] + 1.0;
                }
                else if(state_last ==4)
                {
                    S4[n/save_jump]   = S4[n/save_jump] + 1.0;
                }
                else if(state_last ==5)
                {
                    S5[n/save_jump]   = S5[n/save_jump] + 1.0;
                }
                else if(state_last ==6)
                {
                    S6[n/save_jump]   = S6[n/save_jump] + 1.0;
                }
                else if(state_last ==7)
                {
                    S7[n/save_jump]   = S7[n/save_jump] + 1.0;
                }
                else if(state_last ==8)
                {
                    S8[n/save_jump]   = S8[n/save_jump] + 1.0;
                }
                else if(state_last ==9)
                {
                    S9[n/save_jump]   = S9[n/save_jump] + 1.0;
                }
                else if(state_last ==10)
                {
                    S10[n/save_jump]  = S10[n/save_jump] + 1.0;
                }
                else if(state_last ==11)
                {
                    S11[n/save_jump]  = S11[n/save_jump] + 1.0;
                }
                else if(state_last ==12)
                {
                    S12[n/save_jump]  = S12[n/save_jump] + 1.0;
                }
            }
            if (output_last_count == save_jump)
            {
                output_last_count = 1;
            }
            else
            {
                output_last_count+=1;
            }
        } // end the (n-loop) of the time marching
        
        
    } // end the (r-loop) of the repeat step
    
    
    //last_counting total number of molecules in each state at every time step
    //--------------------------------------------------------------------------------------
    // Calculate The Steady-State Force using Impluse using data from the last 0.5 sec
    // (i.e., just 10000 time steps, considering we only saved every 10 timesteps) only using numerical trapaziodal integration
    //--------------------------------------------------------------------------------------
    
    for (int n = tsteps-10000; n < tsteps-1; n++)  // time marching
    {
        S0_temp_last  = S0_temp_last   +(S0[n/save_jump]  /n_SERCA_Molecules);
        S1_temp_last  = S1_temp_last   +(S1[n/save_jump]  /n_SERCA_Molecules);
        S2_temp_last  = S2_temp_last   +(S2[n/save_jump]  /n_SERCA_Molecules);
        S3_temp_last  = S3_temp_last   +(S3[n/save_jump]  /n_SERCA_Molecules);
        S4_temp_last  = S4_temp_last   +(S4[n/save_jump]  /n_SERCA_Molecules);
        S5_temp_last  = S5_temp_last   +(S5[n/save_jump]  /n_SERCA_Molecules);
        S6_temp_last  = S6_temp_last   +(S6[n/save_jump]  /n_SERCA_Molecules);
        S7_temp_last  = S7_temp_last   +(S7[n/save_jump]  /n_SERCA_Molecules);
        S8_temp_last  = S8_temp_last   +(S8[n/save_jump]  /n_SERCA_Molecules);
        S9_temp_last  = S9_temp_last   +(S9[n/save_jump]  /n_SERCA_Molecules);
        S10_temp_last = S10_temp_last  +(S10[n/save_jump] /n_SERCA_Molecules);
        S11_temp_last = S11_temp_last  +(S11[n/save_jump] /n_SERCA_Molecules);
        S12_temp_last = S12_temp_last  +(S12[n/save_jump] /n_SERCA_Molecules);
        
    }
    
    S0_SS_last  =   S0_temp_last  / 10000;
    S1_SS_last  =   S1_temp_last  / 10000;
    S2_SS_last  =   S2_temp_last  / 10000;
    S3_SS_last  =   S3_temp_last  / 10000;
    S4_SS_last  =   S4_temp_last  / 10000;
    S5_SS_last  =   S5_temp_last  / 10000;
    S6_SS_last  =   S6_temp_last  / 10000;
    S7_SS_last  =   S7_temp_last  / 10000;
    S8_SS_last  =   S8_temp_last  / 10000;
    S9_SS_last  =   S9_temp_last  / 10000;
    S10_SS_last =   S10_temp_last / 10000;
    S11_SS_last =   S11_temp_last / 10000;
    S12_SS_last =   S12_temp_last / 10000;
    
    
    //----------------------------------------------------------
    // Write time-state data into a text file so we can plot it
    //----------------------------------------------------------

    

        
        ss_bound_Ca2[i] = S1_SS_last + S2_SS_last + S10_SS_last + 2* (S3_SS_last + S4_SS_last + S5_SS_last + S6_SS_last + S7_SS_last + S8_SS_last + S9_SS_last);
        
        
        if (ss_bound_Ca2[i] > boundSS_max_temp2)
        {
            boundSS_max_temp2 = ss_bound_Ca2[i];
        }
    }
    
    for (int i2 = 0; i2 < n_pCa; i2++)
    {
        norm_ss_bound_Ca2[i2] = ss_bound_Ca2[i2] / boundSS_max_temp2;
	cout << "Using norm_ss_bound_Ca2" << endl;  
    }
    
    // write steady state info in a file
    std::string filename = "best_residual_SSpCa_Curve.csv";
    
    ofstream time_states_out(filename.c_str()); //opening an output stream for file test.txt
    if(time_states_out.is_open()) //checking whether file could be opened or not.
    {
        // create headers for file
        time_states_out << "pCa,Bound_Ca" << endl; // write the average force
        for (int m = 0; m < n_pCa; m++)  // time marching
        {
            
            time_states_out << calConc_Exp[m] << "," << norm_ss_bound_Ca2[m] << endl; // write the average force
        }
        cout << "Array data succeSS_lastfully saved into the file " << filename << endl;
    }
    
    return;
    
} // end main function
