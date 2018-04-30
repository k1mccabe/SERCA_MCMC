
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
//           ||                                                          [S6]  //      \\  [S8]
//     +Ca   ||                                                      *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
//           ||                                                                \\      //
//           ||                                                        (+ ADP)  \\    //
//           \/                                                                  \\  //
//    (Pi +) E <==> *E-Pi <==> *E-P + Ca <==> *E-P.Ca  <==> *E'-P.Ca + Ca <==>  *E'-P.Ca2
//          [S0]    [S12]      [S11]           [S10]          [S9]                 [S7]
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
//  S9  *E'-P.Ca + Ca       <==> S10 *E-P.Ca + Ca      k.S9.S10  k.S10.S9
//  S10 *E-P.Ca             <==> S11 *E-P + Ca         k.S10.S11 k.S11.S10
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


//const int   tsteps              = 1000001;
const int save_jump = 1000; //how many output values should we keep? To minimize memory usage, we will keep every 10 timepoints.

using namespace std;

float Ca_cyt_conc, Ca_sr_conc, MgATP_conc, MgADP_conc, Pi_conc;

/*  You can't use '.' in variable names. I changed them all to '_'*/
float    k_S1_S0, k_S1_S2,  k_S2_S1, k_S3_S2, k_S3_S4,  k_S4_S3, k_S4_S5, k_S5_S4, k_S5_S6,  k_S6_S5, k_S6_S7, k_S7_S6, k_S5_S8,  k_S8_S5, k_S8_S7, k_S7_S8,  k_S9_S7, k_S9_S10, k_S10_S9,k_S11_S10, k_S11_S12,k_S12_S11,k_S12_S0,k_S0_S12;   // Transition rates
//float k_S0_S1, k_S2_S3,k_S7_S9 , k_S10_S11;

float count_S0, count_S1, count_S2, count_S3, count_S4, count_S5, count_S6, count_S7, count_S8, count_S9, count_S10, count_S11, count_S12;
float S0_SS, S1_SS, S2_SS, S3_SS, S4_SS, S5_SS, S6_SS, S7_SS, S8_SS, S9_SS, S10_SS, S11_SS, S12_SS;
float S0_temp, S1_temp, S2_temp, S3_temp, S4_temp, S5_temp, S6_temp, S7_temp, S8_temp, S9_temp, S10_temp, S11_temp, S12_temp;
int state;
bool open_closed; //O is open 1 is closed
float residual;

//functions to be called
void update_States(int &state, float &dt,
                   float &k_S0_S1, float &k_S0_S12,
                   float &Ca_cyt_conc, float &Ca_sr_conc,
                   float &Pi_conc, float &MgATP_conc, float &MgADP_conc,
                   float &k_S1_S2, float &k_S1_S0,
                   float &k_S2_S3, float &k_S2_S1,
                   float &k_S3_S4, float &k_S3_S2,
                   float &k_S4_S5, float &k_S4_S3,
                   float &k_S5_S6, float &k_S5_S4,
                   float &k_S5_S8, float &k_S8_S5,
                   float &k_S6_S7, float &k_S6_S5,
                   float &k_S7_S9, float &k_S7_S6,
                   float &k_S7_S8, float &k_S8_S7,
                   float &k_S9_S10, float &k_S9_S7,
                   float &k_S10_S11, float &k_S10_S9,
                   float &k_S11_S12, float &k_S11_S10,
                   float &k_S12_S0, float &k_S12_S11
                   );



//--------------------------------------------------------------------------//


float get_Residual(int    & n_SERCA_Molecules,
                   int    & tsteps,
                   float  & dt,
                   int    & n_s,
                   int    & n_pCa,
                   float  & k_S0_S1,
                   float  & k_S2_S3,
                   float  & k_S7_S9,
                   float  & k_S10_S11)
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

    
    
    
    
    
    residual = 0;
    clock_t t1;
    t1=clock();
    float ss_bound_Ca[n_pCa];
    float norm_ss_bound_Ca[n_pCa];
    
    for (int cal = 0; cal < n_pCa; cal++)
    {
        /*----------------------------*/
        /* Assign Model parameters    */
        /*----------------------------*/
        Ca_cyt_conc       = calConc[cal];  // needs citation
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
        srand(time(NULL)); //Random-Seed initialization (must be outside any loop)
        
        
        
        
        //-----------------------
        // SIMULATION FOR SS CURVE
        //-----------------------
        //Creating a vector named tstapes (list of elements of States S0-S12) defining the occupancy of each state according to time at each time step
        //CHANGED - this was taking up too much memory. Now only saving every 10 timesteps.
        //float Time[(tsteps-1)/10];
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
                
                //
                //-----------------------------------------------------------------------------------------------------------------------
                
                update_States(state, dt,
                              k_S0_S1, k_S0_S12,
                              Ca_cyt_conc,  Ca_sr_conc,
                              Pi_conc, MgATP_conc, MgADP_conc,
                              k_S1_S2, k_S1_S0,
                              k_S2_S3, k_S2_S1,
                              k_S3_S4, k_S3_S2,
                              k_S4_S5, k_S4_S3,
                              k_S5_S6, k_S5_S4,
                              k_S5_S8, k_S8_S5,
                              k_S6_S7, k_S6_S5,
                              k_S7_S9, k_S7_S6,
                              k_S7_S8,  k_S8_S7,
                              k_S9_S10,  k_S9_S7,
                              k_S10_S11,  k_S10_S9,
                              k_S11_S12, k_S11_S10,
                              k_S12_S0, k_S12_S11);
                
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
       /* std::string filename = "Time_Data.csv";
        
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
            ss_out << "S6,"   << S6_SS     << endl;
            ss_out << "S7,"   << S7_SS     << endl;
            ss_out << "S8,"   << S8_SS     << endl;
            ss_out << "S9,"   << S9_SS     << endl;
            ss_out << "S10,"  << S10_SS    << endl;
            ss_out << "S11,"  << S11_SS    << endl;
            ss_out << "S12,"  << S12_SS    << endl;
            cout << "Array data successfully saved into the file " << filename2 << endl;
        }
        */
        
        ss_bound_Ca[cal] = S1_SS + S2_SS + S10_SS + 2* (S3_SS+ S4_SS + S5_SS + S6_SS + S7_SS + S8_SS + S9_SS);
    
    
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
        residual_temp = residual_temp + pow ((norm_bound_Ca_Exp[cc] - norm_ss_bound_Ca[cc]/boundSS_max_temp),2); // normalized force
    }
    residual = pow(residual_temp,0.5);

    cout << " residual : " << residual << std::endl;
    
    return residual;
    
} // end function
