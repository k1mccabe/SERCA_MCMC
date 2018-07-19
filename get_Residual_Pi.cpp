/*
//-----------------------------------------------------------------------------------------------------
//                     University of piifornia, San Diego
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
//           ||                                                          [S6a] //      \\  [S6]
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
//  S5   E'~P.ADP.Ca2       <==> S6a  *E'-P.ADP.Ca2    k.S5.S6a   k.S6a.S5
//  S6a  *E'-P.ADP.Ca2      <==> S7  *E'-P.Ca2 (+ ADP) k.S6a.S7   k.S7.S6a
//  S5   E'~P.ADP.Ca2       <==> S6   E'~P.Ca2 (+ ADP) k.S5.S6   k.S6.S5
//  S6   E'~P.Ca2           <==> S7  *E'-P.Ca2 (+ ADP) k.S6.S7   k.S7.S6
//  S7  *E'-P.Ca2 (+ ADP)   <==> S8  *E-P.Ca2          k.S7.S8   k.S8.S7
//  S8  *E'-P.Ca + Ca       <==> S9 *E-P.Ca + Ca       k.S8.S9  k.S9.S8
//  S9 *E-P.Ca              <==> S10 *E-P + Ca         k.S9.S10 k.S10.S9
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
float Pi_conc_Pi;

float count_S0_Pi, count_S1_Pi, count_S2_Pi, count_S3_Pi, count_S4_Pi, count_S5_Pi, count_S6a_Pi, count_S7_Pi, count_S6_Pi, count_S8_Pi, count_S9_Pi, count_S10_Pi, count_S11_Pi;
float S0_SS_Pi, S1_SS_Pi, S2_SS_Pi, S3_SS_Pi, S4_SS_Pi, S5_SS_Pi, S6a_SS_Pi, S7_SS_Pi, S6_SS_Pi, S8_SS_Pi, S9_SS_Pi, S10_SS_Pi, S11_SS_Pi;
float S0_temp_Pi, S1_temp_Pi, S2_temp_Pi, S3_temp_Pi, S4_temp_Pi, S5_temp_Pi, S6a_temp_Pi, S7_temp_Pi, S6_temp_Pi, S8_temp_Pi, S9_temp_Pi, S10_temp_Pi, S11_temp_Pi;
int state_Pi;
bool open_closed_Pi; //O is open 1 is closed
float residual_pi;

//functions to be piled
void update_States(int &state_Pi, float &dt,
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

float get_Residual_Pi  (int    & n_SERCA_Molecules,
                     int    & max_tsteps,
                     float  & dt,
                     int    & n_s,
                     int    & n_pPi,
                     float  & k_S0_S1,
                     float  & k_S2_S3,
                     float  & k_S7_S8,
                     float  & k_S9_S10, 
		     float  & k_S5_S6a,
		     float  & k_S6_S7,
		     float  & k_S0_S11,		
		     float  & k_S1_S0, 
		     float  & k_S1_S2,    float  & k_S2_S1,  
		     float  & k_S3_S2,  
		     float  & k_S3_S4,    float  & k_S4_S3, 
		     float  & k_S4_S5,    float  & k_S5_S4, 
		     float  & k_S5_S6,    float  & k_S6_S5,  
	             float  & k_S6a_S5, 
		     float  & k_S6a_S7,   float  & k_S7_S6a, 
		     float  & k_S7_S6, 
	             float  & k_S8_S7,  
	             float  & k_S8_S9,    float  & k_S9_S8, 
	             float  & k_S10_S9, 
		     float  & k_S10_S11,  float  & k_S11_S10, 
	             float  & k_S11_S0, 
		     float  & Ca_sr_conc, float  & MgATP_conc, float  & MgADP_conc, float & Ca_cyt_conc
                     )

{
    float boundSS_max_temp = 0; // this will figure out the highest bound Pi for our loop
    float piConc[13] = { 1.03376779868233E-6,
			2.05352502645715E-6,
			3.37944980307518E-6,
			6.35162720217447E-6,
			1.06867586159251E-5,
			1.44096688378901E-5,
			2.34526541941682E-5,
			3.94597200689256E-5,
			8.37677788439386E-5,
			0.0001778279,
			0.0004079219,
			0.0006493816,
			0.0009357374};

    
    float norm_phosphorylated_Exp[13] = {0.08470588,
					0.17735294,
					0.26794118,
					0.40382353,
					0.54382353,
					0.60764706,
					0.72294118,
					0.82176471,
					0.90617647,
					0.95970588,
					0.98029412,
					0.99264706,
					1};

    
    //CHANGED - this was taking up too much memory. Now only saving every 10 timesteps.
    //float Time[(max_tsteps-1)/10];
    float S0[(max_tsteps-1)/10];         // used to find how many S0 molecules in each iteration ("E")
    float S1[(max_tsteps-1)/10];         //  "E.Ca"
    float S2[(max_tsteps-1)/10];         //  "E'.Ca"
    float S3[(max_tsteps-1)/10];         //  "E'.Ca2"
    float S4[(max_tsteps-1)/10];         //  "E'.ATP.Ca2"
    float S5[(max_tsteps-1)/10];         //  "E'~P.ADP.Ca2"
    float S6a[(max_tsteps-1)/10];        // "*E'~P.ADP.Ca2"
    float S7[(max_tsteps-1)/10];         // "*E'-P.Ca2"
    float S6[(max_tsteps-1)/10];         //  "E'~P.Ca2"
    float S8[(max_tsteps-1)/10];         // "*E-P.Ca2"
    float S9[(max_tsteps-1)/10];         // "*E-P.Ca"
    float S10[(max_tsteps-1)/10];        // "*E-P"
    float S11[(max_tsteps-1)/10];        // "*E-Pi"
    
    
    
    residual_pi = 0;
    float ss_bound_Pi[n_pPi];
    float norm_ss_bound_Pi[n_pPi];
    
    for (int pi = 0; pi < n_pPi; pi++)
    {
            Pi_conc_Pi       = piConc[pi];  // needs citation
        //-----------------------
        // SIMULATION FOR SS CURVE
        //-----------------------
        //Creating a vector named tstapes (list of elements of States S0-S11) defining the occupancy of each state according to time at each time step

        
        // create a temporary count (initial count of state and set it equal to zero. The count will be piculated and set equal to the variables above for each state. This is because C++ will use the last value that was defined, so its a “clear” or reset function
        S0_temp_Pi  = 0;
        S1_temp_Pi  = 0;
        S2_temp_Pi  = 0;
        S3_temp_Pi  = 0;
        S4_temp_Pi  = 0;
        S5_temp_Pi  = 0;
        S6a_temp_Pi = 0;
        S7_temp_Pi  = 0;
        S6_temp_Pi  = 0;
        S8_temp_Pi  = 0;
        S9_temp_Pi  = 0;
        S10_temp_Pi = 0;
        S11_temp_Pi = 0;
        
        S0_SS_Pi  = 0;
        S1_SS_Pi  = 0;
        S2_SS_Pi  = 0;
        S3_SS_Pi  = 0;
        S4_SS_Pi  = 0;
        S5_SS_Pi  = 0;
        S6a_SS_Pi = 0;
        S7_SS_Pi  = 0;
        S6_SS_Pi  = 0;
        S8_SS_Pi  = 0;
        S9_SS_Pi  = 0;
        S10_SS_Pi = 0;
        S11_SS_Pi = 0;
        
        
        //-----------------------------------------------------------------------------------------------------------------//
        // start repeat loop i.e., using r-index - REPEAT PROCESS FOR EACH SERCA MOLECULE
        //------------------------------------------------------------------------------------------------------------------//
        //start in state 0 which is SERCA in open, then go through in time and ##
        for (int rr = 0; rr < n_SERCA_Molecules; rr++) //repetition of whole simulation to smooth curve
        {
            state_Pi = 0; // Note: each time we repeat the simulation, we should set all SERCA to state 0 (EiH2)
            open_closed_Pi = 0; // initially, each molecule is in closed conformation (0)
            //---------------------------------------------------------------------------------------------------------------//
            // start time loop i.e., using n-index
            //----------------------------------------------------------------------------------------------------------------//
            int output_count = 10; //every 10 timepoints, we will save the data (starting with data point 0)
            for (int n = 0; n < max_tsteps; n++)  // time marching
            {  // begin n-loop for time marching
                //
                //setting counts of all states equal to zero to initialize which will later be used to find how many SERCA in each iteration
                //
                count_S0_Pi  = 0.0;  // "E"
                count_S1_Pi  = 0.0;  //  "E.Ca"
                count_S2_Pi  = 0.0;  //  "E'.Ca"
                count_S3_Pi  = 0.0;  //  "E'.Ca2"
                count_S4_Pi  = 0.0;  //  "E'.ATP.Ca2"
                count_S5_Pi  = 0.0;  //  "E'~P.ADP.Ca2"
                count_S6a_Pi = 0.0;  // "*E'~P.ADP.Ca2"
                count_S7_Pi  = 0.0;  // "*E'-P.Ca2"
                count_S6_Pi  = 0.0;  //  "E'~P.Ca2"
                count_S8_Pi  = 0.0;  // "*E-P.Ca2"
                count_S9_Pi  = 0.0;  // "*E-P.Ca"
                count_S10_Pi = 0.0;  // "*E-P"
                count_S11_Pi = 0.0;  // "*E-Pi"
                
                
                //count << "open_closed = " << open_closed << endl;
                
                //update all RUs for current timestep
                
                //
                //-----------------------------------------------------------------------------------------------------------------------
                
                update_States(state_Pi, dt,
                              k_S0_S1, k_S0_S11,
                              Ca_cyt_conc,  Ca_sr_conc,
                              Pi_conc_Pi, MgATP_conc, MgADP_conc,
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
                    if      (state_Pi == 0)
                    {
                        S0[n/save_jump]   = S0[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==1)
                    {
                        S1[n/save_jump]   = S1[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==2)
                    {
                        S2[n/save_jump]   = S2[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==3)
                    {
                        S3[n/save_jump]   = S3[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==4)
                    {
                        S4[n/save_jump]   = S4[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==5)
                    {
                        S5[n/save_jump]   = S5[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==6)
                    {
                        S6a[n/save_jump]   = S6a[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==7)
                    {
                        S7[n/save_jump]   = S7[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==8)
                    {
                        S6[n/save_jump]   = S6[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==9)
                    {
                        S8[n/save_jump]   = S8[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==10)
                    {
                        S9[n/save_jump]  = S9[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==11)
                    {
                        S10[n/save_jump]  = S10[n/save_jump] + 1.0;
                    }
                    else if(state_Pi ==12)
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
        // piculate The Steady-State Force using Impluse using data from the last 0.5 sec
        // (i.e., just 10000 time steps, considering we only saved every 10 timesteps) only using numeripi trapaziodal integration
        //--------------------------------------------------------------------------------------
        
        for (int n = max_tsteps-10000; n < max_tsteps-1; n++)  // time marching
        {
            S0_temp_Pi   = S0_temp_Pi   +(S0[n/save_jump]  /n_SERCA_Molecules);
            S1_temp_Pi   = S1_temp_Pi   +(S1[n/save_jump]  /n_SERCA_Molecules);
            S2_temp_Pi   = S2_temp_Pi   +(S2[n/save_jump]  /n_SERCA_Molecules);
            S3_temp_Pi   = S3_temp_Pi   +(S3[n/save_jump]  /n_SERCA_Molecules);
            S4_temp_Pi   = S4_temp_Pi   +(S4[n/save_jump]  /n_SERCA_Molecules);
            S5_temp_Pi   = S5_temp_Pi   +(S5[n/save_jump]  /n_SERCA_Molecules);
            S6a_temp_Pi  = S6a_temp_Pi  +(S6a[n/save_jump] /n_SERCA_Molecules);
            S7_temp_Pi   = S7_temp_Pi   +(S7[n/save_jump]  /n_SERCA_Molecules);
            S6_temp_Pi   = S6_temp_Pi   +(S6[n/save_jump]  /n_SERCA_Molecules);
            S8_temp_Pi   = S8_temp_Pi   +(S8[n/save_jump]  /n_SERCA_Molecules);
            S9_temp_Pi   = S9_temp_Pi   +(S9[n/save_jump]  /n_SERCA_Molecules);
            S10_temp_Pi  = S10_temp_Pi  +(S10[n/save_jump] /n_SERCA_Molecules);
            S11_temp_Pi  = S11_temp_Pi  +(S11[n/save_jump] /n_SERCA_Molecules);
            
        }
        
        S0_SS_Pi  =   S0_temp_Pi  / 10000;
        S1_SS_Pi  =   S1_temp_Pi  / 10000;
        S2_SS_Pi  =   S2_temp_Pi  / 10000;
        S3_SS_Pi  =   S3_temp_Pi  / 10000;
        S4_SS_Pi  =   S4_temp_Pi  / 10000;
        S5_SS_Pi  =   S5_temp_Pi  / 10000;
        S6a_SS_Pi =   S6a_temp_Pi / 10000;
        S7_SS_Pi  =   S7_temp_Pi  / 10000;
        S6_SS_Pi  =   S6_temp_Pi  / 10000;
        S8_SS_Pi  =   S8_temp_Pi  / 10000;
        S9_SS_Pi  =   S9_temp_Pi  / 10000;
        S10_SS_Pi =   S10_temp_Pi / 10000;
        S11_SS_Pi =   S11_temp_Pi / 10000;
        
        
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
            ss_out << "S6a,"  << S6a_SS    << endl;
            ss_out << "S7,"   << S7_SS     << endl;
            ss_out << "S6,"   << S6_SS     << endl;
            ss_out << "S8,"   << S8_SS     << endl;
            ss_out << "S9,"   << S9_SS     << endl;
            ss_out << "S10,"  << S10_SS    << endl;
            ss_out << "S11,"  << S11_SS    << endl;
            cout << "Array data successfully saved into the file " << filename2 << endl;
        }
        */
        
        ss_bound_Pi[pi] = S5_SS_Pi + S6a_SS_Pi + S6_SS_Pi + S7_SS_Pi + S8_SS_Pi+ S9_SS_Pi + S10_SS_Pi + S11_SS_Pi ;
    
    
            if (ss_bound_Pi[pi] > boundSS_max_temp)
            {
                boundSS_max_temp = ss_bound_Pi[pi];
            }
    }// end loop through Pi
    
    
       //another loop to normalize - divide by largest number
    for (int pi2 = 0; pi2 < n_pPi; pi2++)
    {
        norm_ss_bound_Pi[pi2] = ss_bound_Pi[pi2] / boundSS_max_temp;
    }
    
    //-------------------------------------
    // Formulate Residual/Cost Function :
    //--------------------------------------
    double residual_temp = 0.0;
    for (int cc = 0; cc < n_pPi; cc++)  // Pi-loop
    {
        residual_temp = residual_temp + pow ((norm_phosphorylated_Exp[cc] - norm_ss_bound_Pi[cc]/boundSS_max_temp),2); // normalized force
    }
    residual_pi = pow(residual_temp,0.5);

    cout << " residual from Phosphate : " << residual_pi << std::endl;
    
    return residual_pi;
    
} // end function
