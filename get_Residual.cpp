/*
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
//           ||                                                          [S6a]  //      \\  [S6]
//     +Ca   ||                                                      *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
//           ||                                                                \\      //
//           ||                                                        (+ ADP)  \\    //
//           \/                                                                  \\  //
//    (Pi +) E <==> *E-Pi <==> *E-P + Ca <==> *E-P.Ca  <==> *E'-P.Ca + Ca <==>  *E'-P.Ca2
//          [S0]    [S11]      [S10]           [S9]          [S8]                 [S7]
//
// State Reaction               State Product          Rate(f)   Rate(r)
//  S0   E + Ca             <==> S1   E.Ca             k.S0.S1   k.S1.S0
//  S1   E.Ca               <==> S2   E'.Ca            k.S1.S2   k.S2.S1
//  S2   E'.Ca + Ca         <==> S3   E'.Ca2           k.S2.S3   k.S3.S2
//  S3   E'.Ca2 (+ ATP)     <==> S4   E'.ATP.Ca2       k.S3.S4   k.S4.S3
//  S4   E'.ATP.Ca2         <==> S5   E'~P.ADP.Ca2     k.S4.S5   k.S5.S4
//  S5   E'~P.ADP.Ca2       <==> S6a  *E'-P.ADP.Ca2    k.S5.S6a  k.S6a.S5
//  S6a  *E'-P.ADP.Ca2      <==> S7  *E'-P.Ca2 (+ ADP) k.S6a.S7  k.S7.S6a
//  S5   E'~P.ADP.Ca2       <==> S6   E'~P.Ca2 (+ ADP) k.S5.S6   k.S6.S5
//  S6   E'~P.Ca2           <==> S7  *E'-P.Ca2 (+ ADP) k.S6.S7   k.S7.S6
//  S7  *E'-P.Ca2 (+ ADP)   <==> S8  *E-P.Ca2          k.S7.S8   k.S8.S7
//  S8  *E'-P.Ca + Ca       <==> S9  *E-P.Ca + Ca      k.S8.S9   k.S9.S8
//  S9  *E-P.Ca             <==> S10 *E-P + Ca         k.S9.S10  k.S10.S9
//  S10 *E-P + Ca           <==> S11 *E-Pi             k.S10.12  k.S11.S10
//  S11 *E-Pi               <==> S0   E + (Pi)         k.S11.S0  k.S0.S11
*/

//Libraries defined here
//--------------------------------------------------------------------------
#include <iostream>  // library that contains file input/output functions
#include <fstream>   // library that contains file input/output functions
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include "get_Residual.h"

//const int   max_tsteps              = 1000001;
const int save_jump = 1000; //how many output values should we keep? To minimize memory usage, we will keep every 10 timepoints.

using namespace std;
float Ca_cyt_conc;

float count_S0, count_S1, count_S2, count_S3, count_S4, count_S5, count_S6a, count_S7, count_S6, count_S8, count_S9, count_S10, count_S11;
float S0_SS, S1_SS, S2_SS, S3_SS, S4_SS, S5_SS, S6a_SS, S7_SS, S6_SS, S8_SS, S9_SS, S10_SS, S11_SS;
float S0_temp, S1_temp, S2_temp, S3_temp, S4_temp, S5_temp, S6a_temp, S7_temp, S6_temp, S8_temp, S9_temp, S10_temp, S11_temp;
int state;
bool open_closed; //O is open 1 is closed
float residual;

//functions to be called
void update_States(int &state, float &dt,
                   float &k_S0_S1, float &k_S0_S11,
                   float &Ca_cyt_conc, float &Ca_sr_conc,
                   float &Pi_conc, float &MgATP_conc, float &MgADP_conc,
                   float &k_S1_S2, float &k_S1_S0,
                   float &k_S2_S3, float &k_S2_S1,
                   float &k_S3_S4, float &k_S3_S2,
                   float &k_S4_S5, float &k_S4_S3,
                   float &k_S5_S6a, float &k_S5_S4,
                   float &k_S5_S6, float &k_S6_S5,
                   float &k_S6a_S7, float &k_S6a_S5,
                   float &k_S7_S8, float &k_S7_S6a,
                   float &k_S7_S6, float &k_S6_S7,
                   float &k_S8_S9, float &k_S8_S7,
                   float &k_S9_S10, float &k_S9_S8,
                   float &k_S10_S11, float &k_S10_S9,
                   float &k_S11_S0, float &k_S11_S10
                   );



//--------------------------------------------------------------------------//

float get_Residual(int    & n_SERCA_Molecules,
                   int    & max_tsteps,
                   float  & dt,
                   int    & n_s,
                   int    & n_pCa,
                   float  & k_S0_S1,
                   float  & k_S2_S3,
                   float  & k_S7_S8,
                   float  & k_S9_S10,float  & k_S1_S0, float  & k_S1_S2,  float  & k_S2_S1, float  & k_S3_S2, float  & k_S3_S4,  float  & k_S4_S3, float  & k_S4_S5, float  & k_S5_S4, float  & k_S5_S6a,  float  & k_S6a_S5, float  & k_S6a_S7, float  & k_S7_S6a, float  & k_S5_S6,  float  & k_S6_S5, float  & k_S6_S7, float  & k_S7_S6,  float  & k_S8_S7, float  & k_S8_S9, float  & k_S9_S8,float  & k_S10_S9, float  & k_S10_S11,float  & k_S11_S10,float  & k_S11_S0,float  & k_S0_S11, float  & Ca_sr_conc,float  & MgATP_conc,float  & MgADP_conc,float  & Pi_conc
                   )

{
    float boundSS_max_temp = 0; // this will figure out the highest bound Ca for our loop
    float calConc[16] = {   1.13465021562703E-07,
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
    
    float norm_bound_Ca_Exp[16] = { 0.056698042688369,
                                    0.100474048769127,
                                    0.159057553309407,
                                    0.23871761522272,
                                    0.30582603399111,
                                    0.385634231598201,
                                    0.459159901847406,
                                    0.551640962566692,
                                    0.63566454650982,
                                    0.715472730890509,
                                    0.778419740266281,
                                    0.835003303161443,
                                    0.885304272566714,
                                    0.935510990636819,
                                    0.970991050011721,
                                    1};

    
    //CHANGED - this was taking up too much memory. Now only saving every 10 timesteps.
    //float Time[(max_tsteps-1)/10];
    float S0[(max_tsteps-1)/10];         // used to find how many S0 molecules in each iteration ("E")
    float S1[(max_tsteps-1)/10];         //  "E.Ca"
    float S2[(max_tsteps-1)/10];         //  "E'.Ca"
    float S3[(max_tsteps-1)/10];         //  "E'.Ca2"
    float S4[(max_tsteps-1)/10];         //  "E'.ATP.Ca2"
    float S5[(max_tsteps-1)/10];         //  "E'~P.ADP.Ca2"
    float S6a[(max_tsteps-1)/10];         // "*E'~P.ADP.Ca2"
    float S7[(max_tsteps-1)/10];         // "*E'-P.Ca2"
    float S6[(max_tsteps-1)/10];         //  "E'~P.Ca2"
    float S8[(max_tsteps-1)/10];         // "*E-P.Ca2"
    float S9[(max_tsteps-1)/10];        // "*E-P.Ca"
    float S10[(max_tsteps-1)/10];        // "*E-P"
    float S11[(max_tsteps-1)/10];        // "*E-Pi"
    
    
    
    residual = 0;
    float ss_bound_Ca[n_pCa];
    float norm_ss_bound_Ca[n_pCa];
    
    for (int cal = 0; cal < n_pCa; cal++)
    {
            Ca_cyt_conc       = calConc[cal];  // needs citation
        //-----------------------
        // SIMULATION FOR SS CURVE
        //-----------------------
        //Creating a vector named tstapes (list of elements of States S0-S11) defining the occupancy of each state according to time at each time step

        
        // create a temporary count (initial count of state and set it equal to zero. The count will be calculated and set equal to the variables above for each state. This is because C++ will use the last value that was defined, so its a “clear” or reset function
        S0_temp = 0;
        S1_temp = 0;
        S2_temp = 0;
        S3_temp = 0;
        S4_temp = 0;
        S5_temp = 0;
        S6a_temp = 0;
        S7_temp = 0;
        S6_temp = 0;
        S8_temp = 0;
        S9_temp = 0;
        S10_temp = 0;
        S11_temp = 0;
        
        S0_SS = 0;
        S1_SS = 0;
        S2_SS = 0;
        S3_SS = 0;
        S4_SS = 0;
        S5_SS = 0;
        S6a_SS = 0;
        S7_SS = 0;
        S6_SS = 0;
        S8_SS = 0;
        S9_SS = 0;
        S10_SS = 0;
        S11_SS = 0;
        
        
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
            for (int n = 0; n < max_tsteps; n++)  // time marching
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
                count_S6a  = 0.0;  // "*E'~P.ADP.Ca2"
                count_S7  = 0.0;  // "*E'-P.Ca2"
                count_S6  = 0.0;  //  "E'~P.Ca2"
                count_S8  = 0.0;  // "*E-P.Ca2"
                count_S9 = 0.0;  // "*E-P.Ca"
                count_S10 = 0.0;  // "*E-P"
                count_S11 = 0.0;  // "*E-Pi"
                
                
                //count << "open_closed = " << open_closed << endl;
                
                //update all RUs for current timestep
                
                //
                //-----------------------------------------------------------------------------------------------------------------------
                
                update_States(state, dt,
                              k_S0_S1, k_S0_S11,
                              Ca_cyt_conc,  Ca_sr_conc,
                              Pi_conc, MgATP_conc, MgADP_conc,
                              k_S1_S2, k_S1_S0,
                              k_S2_S3, k_S2_S1,
                              k_S3_S4, k_S3_S2,
                              k_S4_S5, k_S4_S3,
                              k_S5_S6a, k_S5_S4,
                              k_S5_S6, k_S6_S5,
                              k_S6a_S7, k_S6a_S5,
                              k_S7_S8, k_S7_S6a,
                              k_S7_S6,  k_S6_S7,
                              k_S8_S9,  k_S8_S7,
                              k_S9_S10,  k_S9_S8,
                              k_S10_S11, k_S10_S9,
                              k_S11_S0, k_S11_S10);
                
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
                        S6a[n/save_jump]   = S6a[n/save_jump] + 1.0;
                    }
                    else if(state ==7)
                    {
                        S7[n/save_jump]   = S7[n/save_jump] + 1.0;
                    }
                    else if(state ==8)
                    {
                        S6[n/save_jump]   = S6[n/save_jump] + 1.0;
                    }
                    else if(state ==9)
                    {
                        S8[n/save_jump]   = S8[n/save_jump] + 1.0;
                    }
                    else if(state ==10)
                    {
                        S9[n/save_jump]  = S9[n/save_jump] + 1.0;
                    }
                    else if(state ==11)
                    {
                        S10[n/save_jump]  = S10[n/save_jump] + 1.0;
                    }
                    else if(state ==12)
                    {
                        S11[n/save_jump]  = S11[n/save_jump] + 1.0;
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
        
        for (int n = max_tsteps-10000; n < max_tsteps-1; n++)  // time marching
        {
            S0_temp  = S0_temp   +(S0[n/save_jump]  /n_SERCA_Molecules);
            S1_temp  = S1_temp   +(S1[n/save_jump]  /n_SERCA_Molecules);
            S2_temp  = S2_temp   +(S2[n/save_jump]  /n_SERCA_Molecules);
            S3_temp  = S3_temp   +(S3[n/save_jump]  /n_SERCA_Molecules);
            S4_temp  = S4_temp   +(S4[n/save_jump]  /n_SERCA_Molecules);
            S5_temp  = S5_temp   +(S5[n/save_jump]  /n_SERCA_Molecules);
            S6a_temp  = S6a_temp   +(S6a[n/save_jump]  /n_SERCA_Molecules);
            S7_temp  = S7_temp   +(S7[n/save_jump]  /n_SERCA_Molecules);
            S6_temp  = S6_temp   +(S6[n/save_jump]  /n_SERCA_Molecules);
            S8_temp  = S8_temp   +(S8[n/save_jump]  /n_SERCA_Molecules);
            S9_temp = S9_temp  +(S9[n/save_jump] /n_SERCA_Molecules);
            S10_temp = S10_temp  +(S10[n/save_jump] /n_SERCA_Molecules);
            S11_temp = S11_temp  +(S11[n/save_jump] /n_SERCA_Molecules);
            
        }
        
        S0_SS  =   S0_temp  / 10000;
        S1_SS  =   S1_temp  / 10000;
        S2_SS  =   S2_temp  / 10000;
        S3_SS  =   S3_temp  / 10000;
        S4_SS  =   S4_temp  / 10000;
        S5_SS  =   S5_temp  / 10000;
        S6a_SS  =   S6a_temp  / 10000;
        S7_SS  =   S7_temp  / 10000;
        S6_SS  =   S6_temp  / 10000;
        S8_SS  =   S8_temp  / 10000;
        S9_SS =   S9_temp / 10000;
        S10_SS =   S10_temp / 10000;
        S11_SS =   S11_temp / 10000;
        
        
        //----------------------------------------------------------
        // Write time-state data into a text file so we can plot it
        //----------------------------------------------------------
       /* std::string filename = "Time_Data.csv";
        
        ofstream time_states_out(filename); //opening an output stream for file test.txt
        if(time_states_out.is_open()) //checking whether file could be opened or not.
        {
            // create headers for file
            time_states_out << "Time,S0,S1,S2,S3,S4,S5,S6a,S7,S6,S8,S9,S10,S11" << endl; // write the average force
            for (int n = 0; n < (max_tsteps-1)/save_jump; n++)  // time marching
            {
                
                Time [n] = (n*save_jump*dt);
                time_states_out << Time [n] << "," << S0[n] << "," << S1[n] << "," << S2[n] << "," << S3[n] << "," << S4[n] << "," << S5[n] << "," << S6a[n] << "," << S7[n] << "," << S6[n] << "," << S8[n] << "," << S9[n] << "," << S10[n] << "," << S11[n] << endl; // write the average force
            }
            cout << "Array data successfully saved into the file " << filename << endl;
        }*/
        
        //----------------------------------------------------------
        // Write SS data into a text file
        //----------------------------------------------------------
        /*std::string filename2 = "SS_Data.csv";
        
        ofstream ss_out(filename2); //opening an output stream for file test.txt
        if(ss_out.is_open()) //checking whether file could be opened or not.
        {
            
            ss_out << "S0,"   << S0_SS     << endl;
            ss_out << "S1,"   << S1_SS     << endl;
            ss_out << "S2,"   << S2_SS     << endl;
            ss_out << "S3,"   << S3_SS     << endl;
            ss_out << "S4,"   << S4_SS     << endl;
            ss_out << "S5,"   << S5_SS     << endl;
            ss_out << "S6a,"   << S6a_SS     << endl;
            ss_out << "S7,"   << S7_SS     << endl;
            ss_out << "S6,"   << S6_SS     << endl;
            ss_out << "S8,"   << S8_SS     << endl;
            ss_out << "S9,"  << S9_SS    << endl;
            ss_out << "S10,"  << S10_SS    << endl;
            ss_out << "S11,"  << S11_SS    << endl;
            cout << "Array data successfully saved into the file " << filename2 << endl;
        }
        */
        
        ss_bound_Ca[cal] = S1_SS + S2_SS + S9_SS + S8_SS + 2* (S3_SS+ S4_SS + S5_SS + S6a_SS + S7_SS + S6_SS);
    
    
            if (ss_bound_Ca[cal] > boundSS_max_temp)
            {
                boundSS_max_temp = ss_bound_Ca[cal];
            }
    }// end loop through pCa cal
    
    
       //another loop to normalize - divide by largest number
    for (int cal2 = 0; cal2 < n_pCa; cal2++)
    {
        norm_ss_bound_Ca[cal2] = ss_bound_Ca[cal2] / boundSS_max_temp;
    }
    
    //-------------------------------------
    // Formulate Residual/Cost Function :
    //--------------------------------------
    double residual_temp = 0.0;
    for (int cc = 0; cc < n_pCa; cc++)  // Ca-loop
    {
        residual_temp = residual_temp + pow ((norm_bound_Ca_Exp[cc] - norm_ss_bound_Ca[cc]),2); // normalized force
// residual_temp = residual_temp + pow (((norm_bound_Ca_Exp[cc] - norm_ss_bound_Ca[cc])/boundSS_max_temp),2); // normalized force
    }
    residual = pow(residual_temp,0.5);
    cout << " " << std::endl;
    cout << "       Residual being passed on : " << residual << std::endl;
    
    return residual;
    
} // end function
