//genesis  genesis2.2 script

/************************************************************************** 
* This file defines important simulation parameters for GnRH cell  model  *
* with synapses.  Everything uses SI units.                               *
***************************************************************************/

echo loading  constants file, Crank Nicholson method.

/* variables controlling hsolve integration */
float dt = 1.0e-5

int tab_xdivs = 149; int tab_xfills = 2999
/* This line sets up the size of the tabchannels */

/* The model is quite sensitive to these values in NO_INTERP (caclmode=0) */
float tab_xmin = -0.10; float tab_xmax = 0.05; float Ca_tab_max = 0.300

/* make filename for input of compartments to save Vm */
str outputcompsfname = "ksg0423a_outputcomps.txt"

/* Make filename for storing a snapshot of the cell for RESTORE */
str snapshotname = "ksg0423a_snapshot_singlecomp.save"
						 
/*Make filenames for output of Membrane Voltage and Synaptic Conductances*/
str data_fname = "./output/Vm_traces.asc" 
str gsyn_AMPA_fname = "./output/gsyn_AMPA_traces.asc"
str gsyn_GABA_fname = "./output/gsyn_GABA_traces.asc"
str current_data_fname = "./output/Current_traces.asc" 
/*Set up the name for the cell, which will be the hsolve element */
str cellpath = "/gnrh"

/*Set up the name for the Morphology file to be read in, defining the cell */
//str dotp = "variable.p" 
str dotp = "single_comp.p"
/********************  Cell Specific Parameters for KSG0425a **************************/
//float CM =   0.008179 //0.01	// Result of G.A. for Ksg0425a passive short-cip model
//float RMs =  1.754 		// Result of G.A. for Ksg0425a passive short-cip model
//float RMd = {RMs}		// keep RM's the same for now
//float RA =  2.362		// Result of G.A. for Ksg0425a passive short-cip model
/***** preset constants *****/
//float ELEAK = -0.0477 //-0.0513// -0.070	// Ek value used for the leak conductance
//float EREST_ACT = -0.05524  //{ELEAK} // -0.056 //-0.0513 //-0.070	// Vm value used for the //RESET
/************************************************************************************
**   End of cell-specific parameters for KSG0425a:  Uncomment for work on this cell 
************************************************************************************/


/****************************** Cell Specific Parameters for KSG0423a ***********/
/***** cable parameters *****/
float CM =   0.0084075 //0.01	// Result of G.A. for Ksg0423a passive short-cip model
float RMs =  3.06381 		// Result of G.A. for Ksg0423a passive short-cip model
float RMd = {RMs}		// keep RM's the same for now
float RA = 0.966306		// Result of G.A. for Ksg0423a passive short-cip model
/***** preset constants *****/
float ELEAK =-0.062 //-0.0675  //-0.04985 gives RMP = -50 mV RMP
 ///Resting Membrane Potential(a,b) = -50 mV for ELEAK =-0.06976
         //-0.036 from cip response tuning with 20 mV shifted Na and Kdr activations (-38.6 mV RMP)    
                     //-0.048 //-0.048  from model tuning
              // -0.05 //-0.0513
               // -0.070  // Ek value used for the leak conductance
float EREST_ACT = - 0.050 //-0.05549 // {ELEAK}  // -0.059  //-0.070	// Vm value used for the RESET
/************************************************************************************
**   End of cell-specific parameters for KSG0423a:  Uncomment for work on this cell 
************************************************************************************/
float soma_NaF_scale = 9.6 //4 //36 //40 //25 //120
     //Resting Membrane Potential(a) = -50 mV for soma_NaF_scale =25
     //Resting Membrane Potential(b) = -50 mV for soma_NaF_scale =37
     //float soma_NaF_scale = 8
float soma_Kdr_scale =3.2 // 0.56//3.8
 //7.44 //7.0// 5.78 //4
     //Resting Membrane Potential(a) = -50 mV for soma_Kdr_scale =7.14
     //Resting Membrane Potential(b) = -50 mV for soma_Kdr_scale =7.40
     //float soma_Kdr_scale = 45
float axon_NaF_scale = 14 //2
     //float axon_NaF_scale = 25
float axon_Kdr_scale = 50 //50
     //float axon_Kdr_scale = 125
float CaL_scale = 0.400
 //2.88 //2.467 //2.44 // 2.4
     //Resting Membrane Potential(a) = -50 mV for CaL_scale =2.469
     //Resting Membrane Potential(b) = -50 mV for CaL_scale =2.48

/*     Conductance Scaling Values for Ksg0423a with GnRH synapses 
(May not be good for Ksg0425a)*/
float cond_scale_s =1.0
float cond_scale_a = 1
//float cond_scale_s = 1
//float cond_scale_a = 0.41


/*   Conductance Scaling Values for Ksg0425a with GnRH synapses  */
//float cond_scale_s = 0.1
//float cond_scale_a = 1.0


float GNa = 28.87               // 28.87
float GKdr = 17.32		// 17.32
float GCaL = 0.9384             // 0.9384
/***** Active Channels ******/
// only used for proto channels
float ECa = 0.1
float GCa = 1, GK = 1, Gh = 1    //GNa = 1
float ENa = 0.060
float GNaFs = {soma_NaF_scale}*{GNa}*{cond_scale_s}	//433.1  //48043.0 //75000.0
float GNaFd = 0				//48043.0 //27740.0 //75000.0
float GNaFa = {axon_NaF_scale}*{GNa}*{cond_scale_a} 	//48043.0 //75000.0
float EK = -0.085
float GKdrs = {soma_Kdr_scale}*{GKdr}*{cond_scale_s}	//144.36 //6000.0
float GKdrd = 0
float GKdra = {axon_Kdr_scale}*{GKdr}*{cond_scale_a}
float GCaLs = {GCaL}*{CaL_scale}
float GCaLa = {GCaL}*{CaL_scale}

str syn_compt = "soma"			//The compartment with an AMPA synapse
str syn_compt_GABA = "soma"            // The compartment with a GABA synapse
float run_time = 0.5                    // number of seconds to run simulation
//Synaptic conductances
/**********************************************************************
** The parameters for these synapses were copied from J. Edgerton 
** and J. Hanson's GP scripts, specifically "GP1_defaults.g"
** AMPA Synapse parameters revised to reflect data from K. Suter (04)
**********************************************************************/

// Excitatory inputs
float num_AMPA_syns 	= 0 //5		// Numnber of synapses
float syn_strength_AMPA = 2.5 //2.5          //Magnitude of unitary conductance 
float G_AMPA_gnrh    	= 0.5e-9	// base conductance for gnrh AMPA synapse //(0.25e-9 for GP)
float G_AMPA		={syn_strength_AMPA}*{G_AMPA_gnrh}
float tauRise_AMPA 	= 0.0005
float tauFall_AMPA	= 0.0012
float AMPA_freq     = 10  // AMPA-mediated excitation is always 1 Hz unless title indicates other  

float G_NMDA    	= {{G_AMPA}*.05}  //not used yet!
float tauRise_NMDA	= 0.01
float tauFall_NMDA	= 0.03

// Inhibitory inputs
float num_GABA_syns     = 0
float syn_strength_GABA = 1
float G_GABA_gnrh    	= 0.5e-9
float GABA_freq         = 5

float G_GABA            = {syn_strength_GABA}*{G_GABA_gnrh}
float tauRise_GABA	= 0.000764
float tauFall_GABA	= 0.01210

// Reversal potentials
float E_AMPA 		= 0
float E_NMDA 		= 0
float E_GABA 		= -0.0365
