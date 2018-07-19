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
//           /\  \\     ||                                                       //  \\
//           ||   \\     =====================>  E'.ATP.Ca                      //    \\
//           ||    ===========> E.ATP.Ca    //    [S3a]                 [S6a] //      \\  [S6]
//     Ca -> ||             //   [S2a]                               *E'-P.ADP.Ca2      E'~P.Ca2 (+ ADP)
//           ||     E.ATP                                                      \\      //
//           ||     [S1a]                                              (+ ADP)  \\    //
//           \/                                                                  \\  //
//    (Pi +) E <==> *E-Pi <==> *E-P + Ca <==> *E-P.Ca  <==> *E'-P.Ca + Ca <==>  *E'-P.Ca2
//          [S0]    [S11]      [S10]           [S9]          [S8]                 [S7]
//
// State Reaction               State Product          Rate(f)   Rate(r)
//  S0   E + Ca             <==> S1   E.Ca             k_S0_S1   k_S1_S0
//  S0   E(+ATP)            <==> S1a  E.ATP            k_S0_S1a  k_S1a_S0
//  S1   E.Ca               <==> S2   E'.Ca            k_S1_S2   k_S2_S1
//  S1a  E.ATP + Ca         <==> S2a  E.ATP.Ca         k_S1a_S2a k_S2a_S1a
//  S2   E'.Ca + Ca         <==> S3   E'.Ca2           k_S2_S3   k_S3_S2
//  S2a  E.ATP.Ca           <==> S3a  E'.ATP.Ca        k_S2a_S3a k_S3a_S2a
//  S2   E'.Ca  (+ ATP)     <==> S3   E'.ATP.Ca        k_S2_S3a  k_S3a_S2
//  S3   E'.Ca2 (+ ATP)     <==> S4   E'.ATP.Ca2       k_S3_S4   k_S4_S3
//  S3a  E.ATP.Ca           <==> S4   E'.ATP.Ca2       k_S3a_S4  k_S4_S3a
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
//----------------------------------------------------------------------------------------------------------------------
//
//EXAMPLE CODE FOR TWO STATE DECAY:
//
//           if(current_state == $STARTING_STATE_INDEX)
//            {
//                p1 =       ($k.forward                         E * dt); // (first order)        unimolecular forward transition to $FORWARD_STATE
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
// if(state  = 0): Then   S11 [*E-Pi]  <-- S0 {E} --> S1 [E.Ca]                 else stay as E
//          [S1]
//          E.Ca
//           /\
//           ||
//           \/
//    (Pi +) E <==> *E-Pi
//          [S0]    [S11]
//
//__________________________________________________________________________________________________
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//__________________________________________________________________________________________________
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>


void update_States( int &state,       float &dt,
                   float &Ca_cyt_conc,float &Ca_sr_conc,
                   float &Pi_conc,    float &MgATP_conc,
                   float &k_S0_S1,    float &k_S1_S0,
                   float &k_S0_S1a,   float &k_S1a_S0,
                   float &k_S1_S2,    float &k_S2_S1,
                   float &k_S1_S2a,   float &k_S2a_S1,
                   float &k_S1a_S2a,  float &k_S2a_S1a,
                   float &k_S2_S3,    float &k_S3_S2,
                   float &k_S2_S3a,   float &k_S3a_S2,
                   float &k_S2a_S3a,  float &k_S3a_S2a,
                   float &k_S3_S4,    float &k_S4_S3,
                   float &k_S3a_S4,   float &k_S4_S3a,
                   float &k_S4_S5,    float &k_S5_S4,
                   float &k_S5_S6a,   float &k_S6a_S5,
                   float &k_S5_S6,    float &k_S6_S5,
                   float &k_S6a_S7,   float &k_S7_S6a,
                   float &k_S6_S7,    float &k_S7_S6,
                   float &k_S7_S8,    float &k_S8_S7,
                   float &k_S8_S9,    float &k_S9_S8,
                   float &k_S9_S10,   float &k_S10_S9,
                   float &k_S10_S11,  float &k_S11_S10,
                   float &k_S11_S0,   float &k_S0_S11 ,
                   float &MgADP_conc
                   )

{
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
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 0): Then      S1 {E.Ca} <--- S0 [E] --> S11 [E*.Pi]                  else stay as E
            //          [S1]     						    ||
            //          E.Ca 							    \/
            //           /\        E.ATP 			 S1a [E.ATP]	
            //           ||	(+ ATP) [S1a]		
            //           \/ //=======//
            //    (Pi +) E    <==> E*.Pi  
            //          [S0]	   [S11]
            //-----------------------------------------------------------------------------------------------------------------------
            if(current_state == 0)
            {
                p1 =        (k_S0_S1  * Ca_cyt_conc * dt); // (pseudo-first order) bimolecular  forward transition to  E.Ca       [S1]
                p2 = p1 +   (k_S0_S11 * Pi_conc     * dt); //  (pseudo-firstorder) bimolecular backward transition to  E*.Pi      [S11]
                p3 = p2 +   (k_S0_S1a * MgATP_conc  * dt); // (first order)        unimolecular  forward transition to E'~P.Ca2   [S1a] 
                if(randNum < p1)
                {
                    state = 1;   //forward transition to  *E'-P.ADP.Ca2 [S1]
                }
                else if(randNum < p2)
                {
                    state = 12;  //forward transition to  *E-Pi         [S11]
                }
                else if(randNum < p3)
                {
                    state = 13; //reverse transition to    E.ATP        [S1a]
                }
            } //if it is not greater than either probability then stay in the same state
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 1): Then     S0 [E]  <-- S1 {E.Ca} --> S2 [E'.Ca]                  else stay as E.Ca
            //          [S1]      [S2]					|
            //          E.Ca <==> E'.Ca  + Ca 		    \/
            //			 /\							S2a[E.ATP.Ca]
            //           ||
            //           ||  (+ATP)   E.ATP.Ca
            //			 ||           [S2a]	
            //           \/
            //     (Pi +) E
            //          [S0]
            //-----------------------------------------------------------------------------------------------------------------------
            if(current_state == 1)
            {
                p1 =       (k_S1_S2  			  * dt); // (first order)        unimolecular forward transition forward to E'.Ca          [S2]
                p2 = p1 +  (k_S1_S0  * Pi_conc    * dt); // (pseudo first-order) bimolecular backward transition back    to E*.Pi          [S0]
                p3 = p2 +  (k_S1_S2a * MgATP_conc * dt); // (pseudo first order) bimolecular  forward transition forward to E.ATP.Ca       [S2a] 
                if(randNum < p1)  //f your random number is less than first probability , then transition to next state
                {
                    state = 2;   // forward transition to  E'.Ca    [S2]
                }
                else if(randNum < p2) //If the random number is leSS_last than second probability , then transition to the previous state
                {
                    state = 0;  // backward transition to  E        [S0]
                }
                else if(randNum < p3)
                {
                    state = 14; // forward transition to   E.ATP.Ca [S2a]
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
                p1 =       (k_S2_S3 * Ca_cyt_conc * dt); //(pseudo-first order)  bimolecular   forward transition to E'.Ca2             [S3]
                p2 = p1 +  (k_S2_S1 			  * dt); // (first order) 		unimolecular backward transition to E.Ca               [S1]
                p3 = p2 +  (k_S2_S3a * MgATP_conc * dt); //(pseudo-first order)  bimolecular   forward transition to E'.ATP.Ca    	   [S3a]
                if(randNum < p1) //If the random number is less than first probability , then transition to next state         [S3]
                {
                    state = 3;   // forward transition to E'.Ca2    [S3]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state [S0]
                {
                    state = 1; // backward transition to  E.Ca     [S1]
                }
                else if(randNum < p3)
                {
                    state = 13; // forward transition to  E'.ATP.Ca [S3a]
                }

            } //if it is not greater than either probability  then stay in the same state
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
                p1 =       (k_S3_S4 * MgATP_conc * dt);   // (pseudo first-order) bimolecular   forward transition to E'.ATP.Ca2 [S4]
                p2 = p1 +  (k_S3_S2 			 * dt);   // (first order)        unimolecular backward transition to E'.Ca      [S2]
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 4;   // forward transition to E'.ATP.Ca2 [S4]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 2;  // backward transition to E'.Ca      [S2]
                }
            } //if it is not greater than either probability then stay in the same state
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 4): Then     S3 [E'.Ca2] <-- S4 {E'.ATP.Ca2} --> S5 [E'~P.ADP.Ca2]    else stay as E'.ATP.Ca2
            //													|
            //												    \/
            //					
            //
            //
            //   [S3]                 [S4]                 [S5]
            //  E'.Ca2 (+ ATP) <==> E'.ATP.Ca2  <==>   E'~P.ADP.Ca2
            //						    ||				
            //						 E'.ATP.Ca
            //						   [S3a]
            //-----------------------------------------------------------------------------------------------------------------------
            else if(current_state == 4)
            {
                p1 =       (k_S4_S5 * MgATP_conc * dt); // (pseudo first-order) bimolecular   forward transition to E'~P.ADP.Ca2 [S5]
                p2 = p1 +  (k_S4_S3              * dt); // (first order)        unimolecular backward transition to E'.Ca2       [S3]
                p3 = p2 +  (k_S4_S3a 			 * dt); // (first order) 		unimolecular backward transition to E'.ATP.Ca    [S3a] 
                if(randNum < p1) //If the random number is less than first probability , then transition to next state
                {
                    state = 5;   // forward transition to E'~P.ADP.Ca2 [S5]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 3;  // backward transition to E'.Ca2       [S3]
                }
                else if(randNum < p3)
                {
                    state = 15; // backward transition to E'.ATP.Ca    [S3a]
                }
            } //if it is not greater than either probability  then stay in the same state
            



            /*-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 5):  Then    S4 [E'.ATP.Ca2]  <--  S5 {E'~P.ADP.Ca2}  --> S6a [*E'-P.ADP.Ca2] else stay as E'~P.ADP.Ca2
            //                                                        ||
            //                                                        \/
            //                                                 S6  [E'~P.Ca2]
            //    [S4]           [S5]
            // E'.ATP.Ca2 <==> E'~P.ADP.Ca2
            //                  //  \\
            //                 //    \\
            //         [S6a]   \/      \/ [S6]
            //       *E'-P.ADP.Ca2   E'~P.Ca2 (+ ADP)
            //
            //
            //--------------------------------------------------------------------------------------------------------------------
*/
            else if(current_state == 5)
            {
                p1 =        (k_S5_S6a * dt); // (first order)        unimolecular  forward transition to *E'-P.ADP.Ca2               [S6a]
                p2 = p1 +   (k_S5_S4  * dt); // (first order)        unimolecular backward transition to  E'.ATP.Ca2                 [S4]
                p3 = p2 +   (k_S5_S6  * dt); // (first order)        unimolecular  forward transition to  E'~P.Ca2                   [S6]
                if(randNum < p1)
                {
                    state = 6;   //forward transition to *E'-P.ADP.Ca2 [S6a]
                }
                else if(randNum < p2)
                {
                    state = 8;  //forward transition to  E'~P.Ca2      [S6]
                }
                else if(randNum < p3)
                {
                    state = 4; //reverse transition to   E'.ATP.Ca2    [S4]
                }
            }
            /*-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 6):  Then        S5 [E'~P.ADP.Ca2]  <-- S6a {*E'-P.ADP.Ca2} --> S7 [*E'-P.Ca2]       else stay as *E'-P.ADP.Ca2
            //
            //                     [S5]
            //                  E'~P.ADP.Ca2
            //                  //
            //                 //
            //         [S6a]   \/
            //      *E'-P.ADP.Ca2
            //                /\
            //                 \\
            //                  \\
            //                *E'-P.Ca2
            //                   [S7]
            //
            //----------------------------------------------------------------------------------------------------------------------------
            */
	    else if(current_state == 6)
            {
                p1 =       (k_S6a_S7 * dt);     // (pseudo first-order) bimolecular   forward transition to *E'-P.Ca2                 [S7]
                p2 = p1 +  (k_S6a_S5 * dt);     // (first order)        unimolecular backward transition to  E'~P.ADP.Ca2             [S5]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 7;   // forward transition to *E'-P.Ca2    [S7]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 5;  // backward transition to E'~P.ADP.Ca2 [S5]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state = 7) Then          S6a [*E'-P.ADP.Ca2]  <-- S7 {(+ ADP) *E'-P.Ca2} --> S6 [E'~P.Ca2]  else stay as *E'-P.Ca2
            //
            //               [S6a]              [S6]
            //            *E'-P.ADP.Ca2      E'~P.Ca2
            //                      \\      //
            //                       \\    //
            //                        \\  //
            // *E-P.Ca2 <==> (+ ADP) *E'-P.Ca2
            //    [S8]                 [S7]
            //
            //
            //-------------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 7)
            {
                p1 =      (k_S7_S8              * dt); // (first order)     unimolecular  forward transition to *E-P.Ca2            [S8]
                p2 = p1 + (k_S7_S6a * MgADP_conc * dt); // (first order)     unimolecular backward transition to *E'-P.ADP.Ca2       [S6a]
                p3 = p2 + (k_S7_S6              * dt); // (first order)     unimolecular backward transition to  E'~P.Ca2           [S6]
                if(randNum < p1)
                {
                    state = 9; //forward transition to *E-P.Ca2         [S8]
                }
                else if(randNum < p2)
                {
                    state = 6; //reverse transition to *E'-P.ADP.Ca2    [S6a]
                }
                else if(randNum < p3)
                {
                    state = 8; //reverse transition to  E'~P.Ca2        [S6]
                }
            }
            /*-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state = 8) Then     S5 [E'~P.ADP.Ca2] <-- S6 {(ADP +) E'~P.Ca2 }  --> S7 [*E'-P.Ca2]  else stay as   E'~P.Ca2 
            //
            //        [S5]
            //      E'~P.ADP.Ca2
            //          \\
            //           \\
            //           [S6]
            //        E'~P.Ca2 (+ADP)
            //            //
            //           //
            //          //
            //     *E'-P.Ca2 (+ ADP)
            //        [S7]
            //
            //
            //------------------------------------------------------------------------------------------------------------------------
            */
	    else if(current_state == 8)
            {
                p1 =       (k_S6_S7  * dt);     // (first-order) unimolecular  forward transition to *E'-P.Ca2                      [S7]
                p2 = p1 +  (k_S6_S5  * dt);     // (first order) unimolecular backward transition to E'~P.ADP.Ca2                   [S5]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 7;   // forward transition to *E'-P.Ca2      [S7]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 5;  // backward transition to  E'~P.ADP.Ca2  [S5]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 9): Then      S7 [*E'-P.Ca2] <-- S8 {*E-P.Ca2} --> S9 [*E-P.Ca] +Ca   else stay as 		*E-P.Ca2
            //
            //*E'-P.Ca2 <==  *E-P.Ca2  ==>  *E-P.Ca + Ca
            //  [S7]            [S8]          [S9]
            //
            //------------------------------------------------------------------------------------------------------------------------
            else if(current_state == 9)
            {
                p1 =       (k_S8_S9 * MgATP_conc * dt);     // (pseudo first-order) unimolecular  forward transition to *E-P.Ca    [S9]
                p2 = p1 +  (k_S8_S7               * dt);     // (first order)        unimolecular backward transition to *E'-P.Ca2  [S7]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 10;   // forward transition to *E-P.Ca      [S9]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 7;   // backward transition to *E'-P.Ca2    [S7]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 10): Then   S8 [*E-P.Ca2] <-- S9 {*E-P.Ca} --> S10 [*E-P] +Ca   				else stay as *E-P.Ca
            //
            // *E-P.Ca2 <==  Ca + *E-P.Ca   ==> *E-P + Ca
            //   [S8]             [S9]         [S10]
            //
            //-------------------------------------------------------------------------------------------------------------------------
            
            else if(current_state == 10)
            {
                p1 =       (k_S9_S10             * dt);      //(first-order)       unimolecular  forward transition to *E-P        [S10]
                p2 =       (k_S9_S8 * Ca_sr_conc * dt);     // (pseudo first-order) bimolecular backward transition to *E-P.Ca2    [S8]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 11;   // forward transition to *E-P         [S10]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 9; // backward transition to   *E-P.Ca2     [S8]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            /*-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 11): Then   S9 [*E-P.Ca] <-- S10 {(Ca +) *E-P}--> S11 [*E-Pi]   else stay as *E-P.Ca
            //
            // *E-P.Ca2 <==  Ca + *E-P.Ca    ==> *E-P + Ca
            //   [S8]              [S9]         [S10]
            //
            //-------------------------------------------------------------------------------------------------------------------------*/
            
            else if(current_state ==11)
            {
                p1 =       (k_S10_S11             * dt);    // (first-order)       unimolecular  forward transition to *E-Pi       [S11]
                p2 =       (k_S10_S9 * Ca_sr_conc * dt);    // (pseudo first-order) bimolecular backward transition to *E-P.Ca     [S9]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 12;   // forward transition to *E-Pi        [S11]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 10; // backward transition to  *E-P.Ca      [S9]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            /*-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 12): Then   S10 [(Ca +) *E-P] <-- S11 {*E-Pi} --> S0 [E]  else stay as *E-Pi
            //
            //   *E-P  <==  *E-Pi  ==>  E (+ Pi )
            //   [S10]      [S11]     [S0]
            //
            //------------------------------------------------------------------------------------------------------------------------*/
            else if(current_state ==12)
            {
                p1 =       (k_S11_S0  * dt);    // (first-order)       unimolecular  forward transition to  E                       [S0]
                p2 =       (k_S11_S10 * dt);    // (pseudo first-order) bimolecular backward transition to *E-P                     [S10]
                if(randNum < p1) //If the random number is leSS_last than first probability , then transition to next state
                {
                    state = 0;   // forward transition to  E            [S0]
                }
                else if(randNum < p2)//If the random number is leSS_last than second probability , then transition to previous state
                {
                    state = 11; // backward transition to *E-P          [S10]
                }
            } //if it is not greater than either probability  then stay in the same state
            
            /*----------------------------------------------------------------------------------------------------------------------
            //                              END UPDATE STEP @ EACH TIME
            //-----------------------------------------------------------------------------------------------------------------------
            //
            //
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 13): Then     S0 [E]  <-- S1a {E.ATP} --> S2a [E.ATP.Ca]                  else stay as E.ATP
            //          
            //           					||==>E.ATP.Ca
            //                 ||==>E.ATP<==||   [S2a]
            //             <===||   [S1a]
            //            E
            //          [S0]
            //-----------------------------------------------------------------------------------------------------------------------*/
            if(current_state == 13)
            {
                p1 =       (k_S1a_S2a * Ca_cyt_conc * dt);  // ((pseudo-first order)  bimolecular forward transition forward to E.ATP.Ca  [S2a]
                p2 = p1 +  (k_S1a_S0  * dt); 				// (first-order)          bimolecular backward transition back to   E         [S0]
                if(randNum < p1)  //f your random number is leSS_last than first probability , then transition to next state
                {
                    state = 14;   // forward transition to  E.ATP    [S2a]
                }
                else if(randNum < p2) //If the random number is leSS_last than second probability , then transition to the previous state
                {
                    state = 0;  // backward transition to  E         [S0]
                }
            }
            /*
            //
            //
            //-----------------------------------------------------------------------------------------------------------------------
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 14): Then     S1 [E.Ca]  <-- S2a {E.ATP.Ca} --> S3a [E'.ATP.Ca]                  else stay as E.ATP.Ca
            //                				    	    	    ||
            //           					                    \/
            //		    [S1]							  S1a[E.ATP]
            //          E.Ca<====E.ATP.Ca===> E'.ATP.Ca
            //			          [S2a]			[S3a]
            //           		   ||
            //   	   E.ATP <=====||
            //         [S1a]
            //          
            //-----------------------------------------------------------------------------------------------------------------------*/
            if(current_state == 14)
            {
                p1 =       (k_S2a_S3a  * dt);  // (first order)   unimolecular  forward transition to E'.Ca    [S3a]
                p2 = p1 +  (k_S2a_S1a  * dt);  // (first order)   unimolecular backward transition to E.ATP    [S1a]
                p3 = p2 +  (k_S2a_S1   * dt);  // (first order)   unimolecular backward transition to E.Ca     [S1] 
                if(randNum < p1)  //f your random number is leSS_last than first probability , then transition to next state
                {
                    state = 15;   // forward transition to  E'.ATP.Ca [S3a]
                }
                else if(randNum < p2) //If the random number is leSS_last than second probability , then transition to the previous state
                {
                    state = 13;  // backward transition to  E.ATP     [S1a]
                }
                else if(randNum < p3)
                {
                    state = 1; // backward transition to    E.Ca      [S1]
                }
            }
	
            /*-----------------------------------------------------------------------------------------------------------------------
            //
            // if(state  = 15): Then     S2a [E.ATP.Ca] <-- S3a {E'.ATP.Ca} --> S4 [E'.ATP.Ca2]         else stay as E'.ATP.Ca
            //                		    	    	    	        ||
            //           		     	   	   [S4]		            \/
            //		   					    E'.ATP.Ca2   		S2[E'.Ca]
            //									||
            //          		 E'.Ca<===>E'.ATP.Ca<=====> E.ATP.Ca
            //			          [S2]		   [S3a]         [S2a]
            //         
            //          
            //-----------------------------------------------------------------------------------------------------------------------*/
            if(current_state == 15)
            {
                p1 =       (k_S3a_S4  * Ca_cyt_conc * dt);  // (first order)   unimolecular  forward transition to E'.ATP.Ca2 [S4]
                p2 = p1 +  (k_S3a_S2a 		    * dt);  // (first order)   unimolecular backward transition to E.ATP.Ca   [S2a]
                p3 = p2 +  (k_S3a_S2  		    * dt);  // (first order)   unimolecular backward transition to E'.Ca      [S2] 
                if(randNum < p1)  //if your random number is less than first probability , then transition to next state
                {
                    state = 4;   // forward transition to  E'.ATP.Ca2 [S4]
                }
                else if(randNum < p2) //If the random number is less than second probability , then transition to the previous state
                {
                    state = 14;  // backward transition to E.ATP.Ca   [S2a]
                }
                else if(randNum < p3)
                {
                    state = 2;  // backward transition to  E'.Ca      [S2]
                }
            }//if it is not greater than either probability then stay in the same state
			//
    return;
}
