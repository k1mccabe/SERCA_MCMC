
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


using namespace std;

float Ca_cyt_conc, Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc;

/*  You can't use '.' in variable names. I changed them all to '_'*/
float   k_S1_S0, k_S1_S2,  k_S2_S1,
        k_S3_S2, k_S3_S4,  k_S4_S3,
        k_S4_S5, k_S5_S4,  k_S5_S6,  k_S6_S5,
        k_S6_S7, k_S7_S6,  k_S5_S8,  k_S8_S5,
        k_S8_S7, k_S7_S8,  k_S9_S7,
        k_S9_S10,k_S10_S9, k_S11_S10,
        k_S11_S12,k_S12_S11,k_S12_S0,k_S0_S12;   // Transition rates

float count_S0, count_S1, count_S2, count_S3, count_S4, count_S5, count_S6, count_S7, count_S8, count_S9, count_S10, count_S11, count_S12;
float S0_SS, S1_SS, S2_SS, S3_SS, S4_SS, S5_SS, S6_SS, S7_SS, S8_SS, S9_SS, S10_SS, S11_SS, S12_SS;
float S0_temp, S1_temp, S2_temp, S3_temp, S4_temp, S5_temp, S6_temp, S7_temp, S8_temp, S9_temp, S10_temp, S11_temp, S12_temp;
int state;
bool open_closed; //O is open 1 is closed



//--------------------------------------------------------------------------//


void lastRun  		(int    & n_SERCA_Molecules,
                     int    & tsteps,
                     float  & dt,
                     int    & n_s,
                     float  & k_S0_S1,
                     float  & k_S2_S3,
                     float  & k_S7_S9,
                     float  & k_S10_S11
                     )
{
   
   int save_jump = 100; //how many output values should we keep? To minimize memory usage, we will keep every 10 timepoints.
     
     /*----------------------------*/
    /* Assign Model parameters    */
    /*----------------------------*/
    Ca_cyt_conc       = 1e-6;  // needs citation
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
    k_S5_S6           = 800;   // Transition rate of  E'~P.ADP.Ca2 to *E'-P.ADP.Ca2   Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S6_S5           = 200;   // Transition rate of *E'-P.ADP.Ca2 to E'~P.ADP.Ca2    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S6_S7           = 500;   // Transition rate of *E'-P.ADP.Ca2 to *E'-P.Ca2       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S7_S6           = 4e6;   // Transition rate of *E'-P.Ca2 to *E'.ADP-P.Ca2       Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S5_S8           = 6;     // Transition rate of *E'~P.ADP.Ca2 to E'~P.Ca2        Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S8_S5           = 1.25e3;// Transition rate of  E'~P.Ca2 to *E'~P.ADP.Ca2       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S8_S7           = 1;     // Transition rate of  E'~P.Ca2 to *E'-P.Ca2           Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S7_S8           = 10;    // Transition rate of *E'-P.Ca2 to E'~P.Ca2            Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    //k_S7_S9           = 500;   // Transition rate of *E'-P.Ca2 to *E-P.Ca2            Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S9_S7           = 5e5;   // Transition rate of *E-P.Ca2 to *E'-P.Ca2            Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S9_S10          = 20;    // Transition rate of *E-P.Ca2 to *E-P.Ca              Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S10_S9          = 20;    // Transition rate of *E-P.Ca to *E-P.Ca2              Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    //k_S10_S11         = 6e2;   // Transition rate of *E-P.Ca to *E-P                  Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S11_S10         = 6e4;   // Transition rate of *E-P to *E-P.Ca                  Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    k_S11_S12         = 60;    // Transition rate of *E-P to *E-Pi                    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S12_S11         = 60;    // Transition rate of *E-Pi to *E-P                    Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S12_S0          = 6e2;   // Transition rate of *E-Pi to E                       Units (s^-1)      Inesi Methods in Enzymology (1988) 157:154-190
    k_S0_S12          = 1.5e4; // Transition rate of  E to *E-Pi                      Units (M^-1 s^-1) Inesi Methods in Enzymology (1988) 157:154-190
    //end parameter setup
    
    //-------------------------------------------------
    
    

    
    //-----------------------
    // SIMULATION FOR SS CURVE
    //-----------------------
    //Creating a vector named tstapes (list of elements of States S0-S12) defining the occupancy of each state according to time at each time step
    //CHANGED - this was taking up too much memory. Now only saving every 10 timesteps.
    float Time[(tsteps-1)/10];
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

    // create a temporary count (initial count of state and set it equal to zero. The count will be calculated and set equal to the variables above for each state. This is because C++ will use the last value that was defined, so its a “clear” or reset function
    S0_temp = 0;
    S1_temp = 0;
    S2_temp = 0;
    S3_temp = 0;
    S4_temp = 0;
    S5_temp = 0;
    S6_temp = 0;
    S7_temp = 0;
    S8_temp = 0;
    S9_temp = 0;
    S10_temp = 0;
    S11_temp = 0;
    S12_temp = 0;
    
    S0_SS = 0;
    S1_SS = 0;
    S2_SS = 0;
    S3_SS = 0;
    S4_SS = 0;
    S5_SS = 0;
    S6_SS = 0;
    S7_SS = 0;
    S8_SS = 0;
    S9_SS = 0;
    S10_SS = 0;
    S11_SS = 0;
    S12_SS = 0;

        
     //-----------------------------------------------------------------------------------------------------------------//
    // start repeat loop i.e., using r-index - REPEAT PROCESS FOR EACH SERCA MOLECULE
    //------------------------------------------------------------------------------------------------------------------//
    //start in state 0 which is SERCA in open, then go through in time and ##
    for (int rr = 0; rr < n_SERCA_Molecules; rr++) //repetition of whole simulation to smooth curve
    {
        state = 0; // Note: each time we repeat the simulation, we should set all SERCA to state 0 (EiH2)
        open_closed = 0; // initially, each molecule is in closed conformation (0)
        //---------------------------------------------------------------------------------------------------------------//
        // start time loop i.e., using n-index
        //----------------------------------------------------------------------------------------------------------------//
        int output_count = 10; //every 10 timepoints, we will save the data (starting with data point 0)
        for (int n = 0; n < tsteps; n++)  // time marching
        {  // begin n-loop for time marching
            //
            //setting counts of all states equal to zero to initialize which will later be used to find how many SERCA in each iteration
            //
            count_S0  = 0.0;   // "E"
            count_S1  = 0.0;  //  "E.Ca"
            count_S2  = 0.0;  //  "E'.Ca"
            count_S3  = 0.0;  //  "E'.Ca2"
            count_S4  = 0.0;  //  "E'.ATP.Ca2"
            count_S5  = 0.0;  //  "E'~P.ADP.Ca2"
            count_S6  = 0.0;  // "*E'~P.ADP.Ca2"
            count_S7  = 0.0;  // "*E'-P.Ca2"
            count_S8  = 0.0;  //  "E'~P.Ca2"
            count_S9  = 0.0;  // "*E-P.Ca2"
            count_S10 = 0.0;  // "*E-P.Ca"
            count_S11 = 0.0;  // "*E-P"
            count_S12 = 0.0;  // "*E-Pi"
            
            
            //count << "open_closed = " << open_closed << endl;
            
            //update all RUs for current timestep
            
            int current_state;
            float p1, p2, p3; //Probabilities of each state
            float randNum = 0;
            //
            //-------------------------------------------------------------------------------------------------------------------------
            current_state = state;     // get current state ID number (e.g. 0 = E , etc)
            //
            //setting index equal from 0 to 12 for each state
            //
            //
            // generate random numbers
            randNum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); // Generate a random number between 0 and 1
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // This portion of the code is used to update the states of each RUs based on the Markov step
            //
            //-----------------------------------------------------------------------------------------------------------------------
            // EXAMPLE COMMENT :
            // if(state  = $STARTING_STATE): Then   $REVERSE_STATE {$REVERSE_STATE_NAME}  <-- $STARTING_STATE_NAME{$GIVEN_NAME} --> S2 {$FORWARD_STATE_NAME}    else stay as GIVEN_NAME
            //
            //          [S1]       [S2]             [S3]                 [S4]                 [S5]
            //          E.Ca <==> E'.Ca  + Ca <==> E'.Ca2 (+ ATP) <==> E'.ATP.Ca2  <==>   E'~P.ADP.Ca2
            //           /\                                                                  //  \\
            //           ||                                                                 //    \\
            //           ||                                                          [S6]  //      \\  [S8]
            //           ||                                                      *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
            //           ||                                                                \\      //
            //           ||                                                        (+ ADP)  \\    //
            //           \/                                                                  \\  //
            //    (Pi +) E <==> *E-Pi <==> *E-P + Ca <==> *E-P.Ca + Ca <==> *E-P.Ca2 <==>  *E'-P.Ca2
            //          [S0]    [S12]      [S11]           [S10]              [S9]            [S7]
            //
            //----------------------------------------------------------------------------------------------------------------------
            //
            //EXAMPLE CODE FOR TWO STATE DECAY:
            //
            //           if(current_state == $STARTING_STATE_INDEX)
            //            {
            //                p1 =       ($k.forward                  * dt); // (first order)        unimolecular forward transition to $FORWARD_STATE
            //                p2 = p1 +  ($k.reverse * $species_concentration * dt); // (pseudo first-order) bimolecular backward transition to $REVERSE_STATE
            //
            //                if(randNum < p1)  // If your random number is less than first probability , then transition to next state
            //                {
            //                    state = $FORWARD_STATE_INDEX;   //  transition to $FORWARD_STATE_NAME
            //                }
            //                else if(randNum < p2) //If the random number is less than second probability , then transition to previous state
            //                {
            //                    state = $REVERSE_STATE_INDEX; //    transition to $REVERSE_STATE_NAME
            //                }
            //            } //if it is not greater than either probability then stay in the same state
            //
            //EXAMPLE CODE FOR THREE STATE DECAY
            //--------------------------------------------------------------------------------------------------------
            // if(state  = 2): Then                                   $REVERSE_STATE_NAME-1
            //                                                                 /\
            //                                                                 ||
            //                    $REVERSE_STATE_NAME-0 + $species <-- [$STARTING_STATE_NAME] --> $FORWARD_STATE_NAME     else stay as $STARTING_STATE_NAME]
            //--------------------------------------------------------------------------------------------------------
            //            if(current_state == 2)
            //            {
            //                p1 =        ($k.forward                            * dt); // (first order)        unimolecular forward  transition to $FORWARD_STATE_NAME   $FORWARD_STATE_INDEX
            //                p2 = p1 +   ($k.reverse-0                          * dt); // (pseudo first-order) unimolecular backward transition to $REVERSE_STATE_NAME-0 $REVERSE_STATE_INDEX-0
            //                p3 = p2 +   ($k.reverse-1 * $species_concentration * dt); // (first order)        bimolecular forward   transition to $FORWARD_STATE_NAME-1   $FORWARD_STATE_INDEX-1
            //              if(randNum < p1)
            //              {
            //                    state = FORWARD_STATE_INDEX;   //  transition to $FORWARD_STATE_NAME
            //                }
            //                else if(randNum < p2)
            //                {
            //                    state = REVERSE_STATE_INDEX-0; //transition to $REVERSE_STATE_NAME-0
            //                }
            //                else if(randNum < p3)
            //                {
            //                    state = REVERSE_STATE_INDEX-1; //transition to $REVERSE_STATE_NAME-1
            //               }
            //            }
            //
            //EXAMPLE CODE FOR ONE STATE DECAY or "DEAD-END STATE" outside of the loop
            //---------------------------------------------------------------------------------------------------------
            // if(state  = 13): Then       [$STARTING_STATE_NAME]          else stay as $FORWARD_STATE_NAME
            //                                      \/
            //                             $FORWARD_STATE_NAME
            //---------------------------------------------------------------------------------------------------------
            //if(current_state == 13)
            //            {
            //                p1 =        ($k.forward * dt); // transition to $FORWARD_STATE_NAME $
            //    SPH         p2 = 1 - p1
            //                if(randNum < p1)
            //                {
            //                    state = $FORWARD_STATE_INDEX;  //transition to $FORWARD_STATE_NAME $FORWARD_STATE_INDEX
            //                }
            //    SPH         else if(randNum < p2)
            //    SPH           {
            //    SPH             state = $REVERSE_STATE_INDEX; // transition to $REVERSE_STATE_NAME $REVERSE_STATE_INDEX
            //               }
            //            }
            //    SPH for KJM       How does it come out of this state if not with the new lines that I put in?
            //      4/10/18             I think thats why all the molecules were ending up in the dead end state.
            //                              Not sure though..
            //
            //_________________________________________________________________________________________________
            //
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //_________________________________________________________________________________________________
            //
            //
            //-------------------------------------------------------------------------------------------------
            //
            //                      SERCA MONTE CARLO MARKOV CHAIN LOOP INTEGRATION STARTS HERE:
            //
            //-------------------------------------------------------------------------------------------------
            //
            // if(state  = 0): Then   S12 [*E-Pi]  <-- S0 {E} --> S1 [E.Ca]                 else stay as E
            //          [S1]
            //          E.Ca
            //           /\
            //           ||
            //           \/
            //    (Pi +) E <==> *E-Pi
            //          [S0]    [S12]
            //
            //__________________________________________________________________________________________________
            //
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //__________________________________________________________________________________________________
            if(current_state == 0)
            {
                p1 =       (k_S0_S1  * Ca_cyt_conc * dt); // (pseudo first-order) bimolecular  forward transition to  E.Ca          [S1]
                p2 = p1 +  (k_S0_S12 * Pi_conc     * dt); // (pseudo first-order) bimolecular backward transition to *E-Pi          [S12]
                if(randNum < p1)  //If the random number is less than first probability, then transition to the next state
                {
                    state = 1;   // forward transition to  E.Ca [S1]
                }
                else if(randNum < p2) //If the random number is less than second probability , then transition to the previous state
                {
                    state = 12; // backward transition to *E-Pi [S12]
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
                if(randNum < p1)  //f your random number is less than first probability , then transition to next state
                {
                    state = 2;   // forward transition to  E'.Ca    [S2]
                }
                else if(randNum < p2) //If the random number is less than second probability , then transition to the previous state
                {
                    state = 0;
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state              [S3]
                {
                    state = 3;   // forward transition to E'.Ca2    [S3]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state [S0]
                {
                    state = 1; // backward transition to  E'.Ca     [S1]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 4;   // forward transition to E'.ATP.Ca2 [S4]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 2;  // backward transition to E'.Ca      [S2]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 5;   // forward transition to E'~P.ADP.Ca2 [S5]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 3;  // backward transition to E'.Ca2       [S3]
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
                    state = 6;   //forward transition to *E'-P.ADP.Ca2 [S6]
                }
                else if(randNum < p2)
                {
                    state = 8;  //forward transition to  E'~P.Ca2      [S8]
                }
                else if(randNum < p3)
                {
                    state = 4; //reverse transition to   E'.ATP.Ca2    [S4]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 7;   // forward transition to *E'-P.Ca2    [S7]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 5;  // backward transition to E'~P.ADP.Ca2 [S5]
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
                    state = 9; //forward transition to *E-P.Ca2         [S9]
                }
                else if(randNum < p2)
                {
                    state = 6; //reverse transition to *E'-P.ADP.Ca2    [S6]
                }
                else if(randNum < p3)
                {
                    state = 8; //reverse transition to  E'~P.Ca2        [S8]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 7;   // forward transition to *E'-P.Ca2      [S7]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 5;  // backward transition to  E'~P.ADP.Ca2  [S5]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 10;   // forward transition to *E-P.Ca      [S10]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 7;   // backward transition to *E'-P.Ca2    [S7]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 11;   // forward transition to *E-P         [S11]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 9; // backward transition to   *E-P.Ca2     [S9]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 12;   // forward transition to *E-Pi        [S12]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 10; // backward transition to  *E-P.Ca      [S10]
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
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 0;   // forward transition to  E            [S0]
                }
                else if(randNum < p2)//If the random number is less than second probability , then transition to previous state
                {
                    state = 11; // backward transition to *E-P          [S11]
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
            
            if (output_count == save_jump)
            {
                if      (state == 0)
                {
                    S0[n/save_jump]   = S0[n/save_jump] + 1.0;
                }
                else if(state ==1)
                {
                    S1[n/save_jump]   = S1[n/save_jump] + 1.0;
                }
                else if(state ==2)
                {
                    S2[n/save_jump]   = S2[n/save_jump] + 1.0;
                }
                else if(state ==3)
                {
                    S3[n/save_jump]   = S3[n/save_jump] + 1.0;
                }
                else if(state ==4)
                {
                    S4[n/save_jump]   = S4[n/save_jump] + 1.0;
                }
                else if(state ==5)
                {
                    S5[n/save_jump]   = S5[n/save_jump] + 1.0;
                }
                else if(state ==6)
                {
                    S6[n/save_jump]   = S6[n/save_jump] + 1.0;
                }
                else if(state ==7)
                {
                    S7[n/save_jump]   = S7[n/save_jump] + 1.0;
                }
                else if(state ==8)
                {
                    S8[n/save_jump]   = S8[n/save_jump] + 1.0;
                }
                else if(state ==9)
                {
                    S9[n/save_jump]   = S9[n/save_jump] + 1.0;
                }
                else if(state ==10)
                {
                    S10[n/save_jump]  = S10[n/save_jump] + 1.0;
                }
                else if(state ==11)
                {
                    S11[n/save_jump]  = S11[n/save_jump] + 1.0;
                }
                else if(state ==12)
                {
                    S12[n/save_jump]  = S12[n/save_jump] + 1.0;
                }
            }
            if (output_count == save_jump)
            {
                output_count = 1;
            }
            else
            {
                output_count+=1;
            }
        } // end the (n-loop) of the time marching
        
        
    } // end the (r-loop) of the repeat step
    
    
    //counting total number of molecules in each state at every time step
    //--------------------------------------------------------------------------------------
    // Calculate The Steady-State Force using Impluse using data from the last 0.5 sec
    // (i.e., just 10000 time steps, considering we only saved every 10 timesteps) only using numerical trapaziodal integration
    //--------------------------------------------------------------------------------------
    
    for (int n = tsteps-10000; n < tsteps-1; n++)  // time marching
    {
        S0_temp  = S0_temp   +(S0[n/save_jump]  /n_SERCA_Molecules);
        S1_temp  = S1_temp   +(S1[n/save_jump]  /n_SERCA_Molecules);
        S2_temp  = S2_temp   +(S2[n/save_jump]  /n_SERCA_Molecules);
        S3_temp  = S3_temp   +(S3[n/save_jump]  /n_SERCA_Molecules);
        S4_temp  = S4_temp   +(S4[n/save_jump]  /n_SERCA_Molecules);
        S5_temp  = S5_temp   +(S5[n/save_jump]  /n_SERCA_Molecules);
        S6_temp  = S6_temp   +(S6[n/save_jump]  /n_SERCA_Molecules);
        S7_temp  = S7_temp   +(S7[n/save_jump]  /n_SERCA_Molecules);
        S8_temp  = S8_temp   +(S8[n/save_jump]  /n_SERCA_Molecules);
        S9_temp  = S9_temp   +(S9[n/save_jump]  /n_SERCA_Molecules);
        S10_temp = S10_temp  +(S10[n/save_jump] /n_SERCA_Molecules);
        S11_temp = S11_temp  +(S11[n/save_jump] /n_SERCA_Molecules);
        S12_temp = S12_temp  +(S12[n/save_jump] /n_SERCA_Molecules);
        
    }
    
    S0_SS  =   S0_temp  / 10000;
    S1_SS  =   S1_temp  / 10000;
    S2_SS  =   S2_temp  / 10000;
    S3_SS  =   S3_temp  / 10000;
    S4_SS  =   S4_temp  / 10000;
    S5_SS  =   S5_temp  / 10000;
    S6_SS  =   S6_temp  / 10000;
    S7_SS  =   S7_temp  / 10000;
    S8_SS  =   S8_temp  / 10000;
    S9_SS  =   S9_temp  / 10000;
    S10_SS =   S10_temp / 10000;
    S11_SS =   S11_temp / 10000;
    S12_SS =   S12_temp / 10000;
    
    
    //----------------------------------------------------------
    // Write time-state data into a text file so we can plot it
    //----------------------------------------------------------
    std::string filename = "Time_Data_gbest.csv";
    
    ofstream time_states_out(filename); //opening an output stream for file test.txt
    if(time_states_out.is_open()) //checking whether file could be opened or not.
    {
        // create headers for file
        time_states_out << "Time,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12" << endl; // write the average force
        for (int n = 0; n < (tsteps-1)/save_jump; n++)  // time marching
        {
            
            Time [n] = (n*save_jump*dt);
            time_states_out << Time [n] << "," << S0[n] << "," << S1[n] << "," << S2[n] << "," << S3[n] << "," << S4[n] << "," << S5[n] << "," << S6[n] << "," << S7[n] << "," << S8[n] << "," << S9[n] << "," << S10[n] << "," << S11[n] << "," << S12[n] << endl; // write the average force
        }
        cout << "Array data successfully saved into the file " << filename << endl;
    }
    
    //----------------------------------------------------------
    // Write SS data into a text file
    //----------------------------------------------------------
    std::string filename2 = "SS_Data_gbest.csv";
    
    ofstream ss_out(filename2); //opening an output stream for file test.txt
    if(ss_out.is_open()) //checking whether file could be opened or not.
    {
        
        ss_out << "S0,"   << S0_SS     << endl;
        ss_out << "S1,"   << S1_SS     << endl;
        ss_out << "S2,"   << S2_SS     << endl;
        ss_out << "S3,"   << S3_SS     << endl;
        ss_out << "S4,"   << S4_SS     << endl;
        ss_out << "S5,"   << S5_SS     << endl;
        ss_out << "S6,"   << S6_SS     << endl;
        ss_out << "S7,"   << S7_SS     << endl;
        ss_out << "S8,"   << S8_SS     << endl;
        ss_out << "S9,"   << S9_SS     << endl;
        ss_out << "S10,"  << S10_SS    << endl;
        ss_out << "S11,"  << S11_SS    << endl;
        ss_out << "S12,"  << S12_SS    << endl;
        cout << "Array data successfully saved into the file " << filename2 << endl;
    }
    return;
    
} // end main function