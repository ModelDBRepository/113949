// genesis

int i,n
str hstr // string for findsolvefield operations


/********************************************************************
** Modified by Carson Roberts 09/2004 From Dieter Jaeger's GP simulations
** (05/2000) for use in gnrh  neurons.  
** In its present form, it  uses the scripts:
** gnrh_const.g, gnrh_actcomps.g,  gnrh_chansave.g, gnrh_syns.g,
** and gnrh_current_inject.g
*******************************************************************/

setenv SIMPATH {getenv SIMPATH} ./prototypes  

/* always include these default definitions! */
include defaults
//include userprefs

/* simulation constants */
include gnrh_const.g
/*****************************************************************************
** The above file, "gnrh_const.g" holds settings for ELEAK, EREST_ACT,
** The Tabchannel xmin and xmax etc, default values CM,RM,RA (which may be over-
** written by values read in from a "parameter.asc" file).  It is also the place 
** to set the densities of active channels and synapses
******************************************************************************/
echo data_fname {data_fname}
echo gsyn_AMPA_fname {gsyn_AMPA_fname}
echo gsyn_GABA_fname {gsyn_GABA_fname}
/* Skip input of parameters, take values from "gnrh_const" */

/************************************************************/
echo Cell Parameters are: RMs {RMs}, CM {CM}, RA {RA} RMd {RMd}

/* scripts to create the prototypes */

include gnrh_actcomps.g
/******************************************************************************** 
** Sets up prototype Soma, Dendrite and Axon compartments. Differences in channel
** Types and Densities should be set up in this script.
*********************************************************************************/

//include gnrh_chanload.g	// This line for already saved channels
include gnrh_chansave.g		// This line to initially generate Table Channels
/******************************************************************************
** The two lines above will contain the actual kinetic equations for whatever 
** channels are included in the model. 
********************************************************************************/

include gnrh_syns.g		// This script contains the definitions for various synapses

include gnrh_current_inject.g  // This script contains codes for various current injections

/* Create the output elements */

//Create Output element for Membrane Potential 
create asc_file /out
setfield /out filename {data_fname} initialize 1 flush 0 append 1 leave_open 1
useclock /out 9
silent -1

/* To ensure that all subsequent elements are made in the library */
ce /library

//These make the prototypes of channels and compartments that can be
//  invoked in .p files 
make_Kdr_gnrh   // function in gnrh_chansave.g
make_NaF_gnrh   // function in gnrh_chansave.g
make_CaL_gnrh   // function in gnrh_chansave.g
make_gnrh_comps   // function in gnrh_actcomps.g

/* create the model and set up the run cell mode */
// read cell data from .p file AND set the cell itself up as the hsolve element
readcell {dotp} {cellpath} -hsolve
echo done reading cell


/* Set the clocks */
for (i=0; i<=8; i=i+1)
    setclock {i} {dt}
end
setclock 9 1.0e-4


// Add Synapses

//create input element tree outside of the cell path
if (!{exists /inputs})
	create neutral /inputs
end
create neutral /inputs/AMPAs
create neutral /inputs/GABAs

/******************************************************************/
/*   Block of code to create AMPA Synapses, and set up timetables */
/* to drive them.                                                 */
/******************************************************************/

/* Generate the prototype AMPA synapse, as defined in gnrh_syns */
make_gnrh_AMPA

randseed = 1234567

/* create a large number of independant  synapses,  */
/*  receiving activation info from timetables */

/*AMPA synapses*/
create neutral /inputs/AMPAs/{syn_compt}
for (i=0; i<{num_AMPA_syns}; i=i+1)  
     copy /library/AMPA {cellpath}/{syn_compt}/AMPAs{i}

     addmsg  {cellpath}/{syn_compt}/AMPAs{i} {cellpath}/{syn_compt} CHANNEL Gk Ek
     addmsg  {cellpath}/{syn_compt} {cellpath}/{syn_compt}/AMPAs{i} VOLTAGE Vm
     
    create timetable /inputs/AMPAs/{syn_compt}/AMPAtt{i}
    // set properties of spike timing
     setfield /inputs/AMPAs/{syn_compt}/AMPAtt{i}   \
           maxtime 10 \                   // Max. Time up to which table is filled
	   act_val 1.0 \                 // Value of Activation when set to ON
	   method 2  \                   // Gamma Dist. for Inter Time Intervals
	     meth_desc1 {1/{AMPA_freq}} \  //mean Inter Time Interval (method 2)
	     meth_desc2 0.003 \            //absolute refractory period of 3 ms (method 2) 
	     meth_desc3 2                  //Order of Gamma Dist. (method 2)

  call /inputs/AMPAs/{syn_compt}/AMPAtt{i} TABFILL
  // addmsg /inputs/AMPAs/{syn_compt}/AMPAtt{i} {cellpath}/{syn_compt}/AMPAs{i} ACTIVATION activation


     //set up spikegen
       create spikegen /inputs/AMPAs/{syn_compt}/spiker{i}
       setfield /inputs/AMPAs/{syn_compt}/spiker{i} \
	     output_amp 1 \
             thresh 0.5
      //connect timetables to AMPA synapses
        addmsg /inputs/AMPAs/{syn_compt}/AMPAtt{i} \
               /inputs/AMPAs/{syn_compt}/spiker{i} INPUT activation
        	addmsg /inputs/AMPAs/{syn_compt}/spiker{i} \
                            {cellpath}/{syn_compt}/AMPAs{i} SPIKE

end

echo we have made {num_AMPA_syns} AMPA syns

/******************************************************************/
/*   Block of code to create GABA Synapses, and set up timetables */
/* to drive them.                                                 */
/******************************************************************/

/* Generate the prototype GABA synapse, as defined in gnrh_syns */
make_gnrh_GABA

randseed = 2345671

/* create a large number of independant  synapses,  */
/*  receiving activation info from timetables */

/*GABA synapses*/
create neutral /inputs/GABAs/{syn_compt_GABA}
for (i=0; i<{num_GABA_syns}; i=i+1)  
     copy /library/GABA {cellpath}/{syn_compt_GABA}/GABAs{i}

     addmsg  {cellpath}/{syn_compt_GABA}/GABAs{i} {cellpath}/{syn_compt_GABA} CHANNEL Gk Ek
     addmsg  {cellpath}/{syn_compt_GABA} {cellpath}/{syn_compt_GABA}/GABAs{i} VOLTAGE Vm
     
    create timetable /inputs/GABAs/{syn_compt_GABA}/GABAtt{i}
    // set properties of spike timing
     setfield /inputs/GABAs/{syn_compt_GABA}/GABAtt{i}   \
           maxtime 10 \                   // Max. Time up to which table is filled
	   act_val 1.0 \                 // Value of Activation when set to ON
	   method 2  \                   // Gamma Dist. for Inter Time Intervals
	     meth_desc1 {1/{GABA_freq}} \  //mean Inter Time Interval (method 2)
	     meth_desc2 0.003 \            //absolute refractory period of 3 ms (method 2) 
	     meth_desc3 3                 //Order of Gamma Dist. (method 2) 
 // for the last variable, use "meth_desc3 3" for inputs independant of the AMPA synapses, and
 // use "meth_desc3 2" for GABA inputs coincident with AMPA ones.
  call /inputs/GABAs/{syn_compt_GABA}/GABAtt{i} TABFILL
  // addmsg /inputs/GABAs/{syn_compt_GABA}/GABAtt{i} \
         // {cellpath}/{syn_compt_GABA}/GABAs{i} ACTIVATION activation

     //set up spikegen
       create spikegen /inputs/GABAs/{syn_compt_GABA}/spiker_GABA{i}
       setfield /inputs/GABAs/{syn_compt_GABA}/spiker_GABA{i} \
	     output_amp 1 \
             thresh 0.5
      //connect timetables to GABA synapses
        addmsg /inputs/GABAs/{syn_compt_GABA}/GABAtt{i} \
               /inputs/GABAs/{syn_compt_GABA}/spiker_GABA{i} INPUT activation
        	addmsg /inputs/GABAs/{syn_compt_GABA}/spiker_GABA{i} \
                            {cellpath}/{syn_compt_GABA}/GABAs{i} SPIKE
end

echo we have made {num_GABA_syns} GABA syns

//	make_gnrh_GABA 
//	copy GABA {cellpath}/soma/GABA 
//	addmsg  {cellpath}/soma/GABA {cellpath}/soma CHANNEL Gk Ek
//	addmsg  {cellpath}/soma {cellpath}/soma/GABA VOLTAGE Vm
//	addmsg /input/Apulse/spiketrain {cellpath}/soma/GABA SPIKE


//set up current injection
/****************************************************************************
* This is a bit of a hack, to get around a Genesis bug involving injecting  *
* current into a cell that has been taken over by the Hines solver.  The    *
* Manual says to use code like:                                             *
*                           setfield {cellpath}/soma inject {curr}          *
*			     call {cellpath} HPUT {cellpath}/soma           *
* But, it turns out that there is a disturbance in the membrane potential   *
* every time the HPUT command is excecuted.  The solution (as implemented   *
* here) is to set up a "pulse" that is actually a DC current injection, and *
* in the run code for current injections to use commands like:              *
*            setfield /pulse level1 {curr}                                  *
*            step 0.05 -time                                                *
*            setfield /pulse level1 0                                       *
* to generate a pulse.  The injections are set up and run from the script   *
* "gnrh_current_inject.g                                                    *
****************************************************************************/

create pulsegen /pulse
setfield /pulse \
        level1          0 \
        width1          0.1      \
        delay1          0      	\
        delay2          0      \
        baselevel       0 	\
        trig_mode       0

addmsg /pulse {cellpath}/soma INJECT output

/****************************
 *  SETUP THE HINES SOLVER  *
 ****************************/
echo preparing hines solver...
//create hsolve {cellpath}
// The hsolve element has already been created.... It is the cell itself
setfield {cellpath} \
        path {cellpath}/##[][TYPE=compartment] \
	comptmode	1 \   // uses less memory, and a bit slower than comptmode 1
        chanmode        4 \   // Mode for efficient saving of multiple compartment values
        calcmode 	0 \   // no linear interpolation of values from lookup tables
	outclock	9 \   // clock to be used for element updates
	storemode   	2     // Total Conductances are stored
        call {cellpath} SETUP
        setmethod 11  // Crank-Nicholson Method

/********************************
 *             OUTPUT           *
 ********************************/

/* Output of Soma trace only:    */
hstr ={findsolvefield {cellpath} {cellpath}/soma Vm}
     addmsg {cellpath} /out SAVE {hstr}
hstr ={findsolvefield {cellpath} {cellpath}/p0[196] Vm}
     addmsg {cellpath} /out SAVE {hstr}

/*  uncomment above line to save only one Vm trace (for simple graphs) */

/************************************************************************/
/* Save Vm of all compartments to a file for use in making movie frames */
/*             Uncomment block below for saving stuff for making movies */
/************************************************************************/
//int n
//str readcompartment
//openfile {outputcompsfname} r
//readcompartment = {readfile {outputcompsfname}}
//while (! {eof {outputcompsfname}})
//     hstr ={findsolvefield {cellpath} {cellpath}/{readcompartment} Vm}
//     addmsg {cellpath} /out SAVE {hstr}
//     readcompartment = {readfile {outputcompsfname}}
//end
//closefile {outputcompsfname}
//

/* ****************************************************************************/
/* Block of code to output synaptic conductances: leave uncommented if not    */
/* desired, to save time in repeated simulations (actually, it does not seem  */
/* to save much time at all......                                             */
/******************************************************************************/
/* Create Output element for AMPA synaptic conductances */
 
if ({num_AMPA_syns} >0)
create asc_file /out2
  setfield /out2 filename {gsyn_AMPA_fname} initialize 1 flush 0 append 1 leave_open 1
  useclock /out2 9
  silent -1
/* Loop through and save each AMPA conductance in the output file */
for (i=0; i<{num_AMPA_syns}; i=i+1)
   hstr ={findsolvefield {cellpath} {cellpath}/{syn_compt}/AMPAs{i} Gk}
   addmsg {cellpath} /out2 SAVE {hstr}
 end
end
  
/* Create Output element for GABA synaptic conductances */
 if ({num_GABA_syns}>0) 
create asc_file /out3
  setfield /out3 filename {gsyn_GABA_fname} initialize 1 flush 0 append 1 leave_open 1
  useclock /out3 9
  silent -1
/* Loop through and save each GABA conductance in the output file */
for (i=0; i<{num_GABA_syns}; i=i+1)
   hstr ={findsolvefield {cellpath} {cellpath}/{syn_compt_GABA}/GABAs{i} Gk}
   addmsg {cellpath} /out3 SAVE {hstr}
 end
end

/* Create Output element for Injected Currents */
create asc_file /out4
  setfield /out4 filename {current_data_fname} initialize 1 flush 0 append 1 leave_open 1
  useclock /out4 9
  silent -1

   hstr ={findsolvefield {cellpath} {cellpath}/soma inject}
   addmsg {cellpath} /out4 SAVE {hstr}


/****************************************************************************/
/*   End of block for saving all voltages and synapses for movie            */
/****************************************************************************/


/*************************************/
/*   Start running the simulation    */
/*************************************/

reset

/*********************************************/
/* Restore a saved snapshot of the cell      */
/* See code below for saving the snapshot    */
/* Uncomment the two lines below to restore  */
/*********************************************/
 restore {snapshotname}
  call {cellpath} HRESTORE


/**********************/
/*  Run for some time */
/**********************/

/*********************************************************************************/
/*  Use the line below to run for some time and with synapses, etc, as set up in */
/* "gnrh_const.g"...  Comment it out for current injection simulations.          */		  
/*********************************************************************************/
		  // step {run_time} -time

/*********************************************************************************/
/*  Use the line below to run for current injection protocols to match           */
/* experiments  as set up in  "gnrh_current_inject.g....                         */ 
/* Comment it out for synaptic activation  simulations.                          */		  
/*********************************************************************************/
				   do_current_injections

/**************************************************/
/*  Take a snapshot of the cell for RESTORE       */
/* Uncomment the lines below to take the snapshot */
/**************************************************/

		  // call {cellpath} HSAVE
/* This updates all compartments from the hsolve element */

		  // save {cellpath}/##[] {snapshotname}
/* This saves all values for all elements in the file "test.save" */
/* NOTE:  Just using the wildcard /## will only save elements of the type "p0[XX]" */
/* but WILL NOT save elements of the type "p1[XX]" */


  quit             // exits Genesis
