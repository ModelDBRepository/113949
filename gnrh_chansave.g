//genesis - gnrh cell  genesis2 script

/*********************************************************************
** The shell of this script is based on the Purkinje cell script of
** De Schutter, as distributed with Genesis.
** This version actually computes all the tables and saves them to files
** E. De Schutter, Caltech, 1991-1992
**********************************************************************/
 

/************************************************************************
** Script to calculate the various parameters used in 
** LeBeau, Van Goor, et al's modification of the Hodgkin-Huxley 
** formalism for various voltage-activated channels in GT1 Neurons.
**
** These calculations are based on: 
** LeBeau, Van Goor, Stojilkovic and Sherman,"Modeling of Membrane 
** Excitability in Gonadotropin-Releasing Hormone-Secreting Neurons Regulated 
** by Ca2+-Mobilizing and Adenylyl Cyclase-Coupled Receptors", 
** J. Neurosci., 20(24):9290-9297 (2000). (Referred to as LBVG(00))
**
** Additional equations come from:
** Van Goor, LeBeau, Krsmanovic, Sherman, Catt, and Stojilkovic, 
** "Amplitude-Dependant Spike-Broadening and Enhanced Ca2+ Signaling in 
** GnRH-Secreting Neurons", Biophysical Journal 79(3):1310-1323 (2000)
** (Referred to as VGLB(00)).
**
**
** These equations have been de-constructed and re-worked by Carson Roberts
** in the essay "GT1 Channels from LBVG(00)" (2004).(unpublished)
**
** Written by Carson Roberts, 11/2004-05/2005
**************************************************************************/



/* Explanation of Temperature Scaling from De Schutter:               */
/**************************************************************************
** The data on which these currents are based were collected at room
**   temperature (anything in the range 20 to 25C).  All rate factors were
**   adapted for body temperature (37C), assuming a Q10 of about 3;
**   in practice all rate factors were multiplied by a factor of 5 
**************************************************************************/

// CONSTANTS
/* should be defined by calling routine (all correctly scaled):
**	ENa, GNa
**	ECa, GCa
**	EK, GK 
**  	Eh, Gh */

int include_gnrh_chansave

if ( {include_gnrh_chansave} == 0 )   /* Higest-Level Loop in Program */

	include_gnrh_chansave = 1

	int i
	float x, dx, y

/***********************************************************************************
** For this model, the Potassium Delayed Rectifier is taken from LBVG(00).  They
** have only "n" gates for Their Kdr, and a constant "h" gate, if used at all. 
*************************************************************************************/


  function make_Kdr_gnrh
        if (({exists Kdr_gnrh}))
                return
        end

        create tabchannel Kdr_gnrh

	float npower =  4 // 4  

   	

        setfield Kdr_gnrh Ek {EK} Gbar {GKdrs} Ik 0 Gk 0\
                Xpower {npower} Ypower 0 Zpower 0
/************************************************************************
** The names of the variables used here are copied from the XPP scripts
** "gt1v2.041.ms1b.ode" and "gt1v3.031.0de" from Artie Sherman, one of the
** co-authors on LBVG(00).
*************************************************************************/

        float snkdr = -1	// This is s (sign) on p.8 of GT1 Channels...
        
/*******************************************************************************/
/* In order to adjust "Spike Threshold" to match more nearly the observed      */
/* values in the gfp-gnrh neurons, it may be useful to shift the activation of */
/* the Kdr to higher voltages.  This will be attempted by changing the         */
/* half-activation value, "vnkdr" in LBVG's model.  Their value (translated     */
/* into Genesis, was vnkdr = 0.025.....                                        */
/*******************************************************************************/
	float vkdr_shift =0.0  //0.0149 //0.022 //0.0031        
	float vkdr_LBVG = 0.025
	float vnkdr = {vkdr_LBVG} -{vkdr_shift}


	//float vnkdr = 0.025	// This is b (V) on p.8 of GT1 Channels....
        float knkdr = 0.015	// This is c (V) on P.8 of GT1 Channels..
	float taunkdrss = 0.015 // This is a (sec) on p.7 of GT1 Channels...
	float vtnkdr = 0.030	// This is b (V) on p.7...
	float ktnkdr = 0.015	// This is c (V) on p.7....
	float dtnkdr = 0.001	// This is d (sec) on p.7 .....
	float rate_factor = 5 //5	
	/* The above is a factor to compensate for LBVG(00) using rate factors measured at Room Temp, and */
	/* this model representing slices measured at 32 degrees (see discussion in De Schutter-Bower (94)) */              

	float  ninf, taun,taun2, alpha, beta

	call Kdr_gnrh TABCREATE X {tab_xdivs} {tab_xmin} {tab_xmax}
	x = {tab_xmin}
	dx = ({tab_xmax} - {tab_xmin})/{tab_xdivs}
	for (i = 0; i <= {tab_xdivs}; i = i + 1)

		ninf = 1 / (1 + {exp { {snkdr}*(x+{vnkdr}) / {knkdr} } } )
		taun = {taunkdrss} / ({exp { (x+{vtnkdr} ) / {ktnkdr} } } + {exp {-1*(x+{vtnkdr}) / {ktnkdr} }}) + {dtnkdr}
		taun2 = {taun}/{rate_factor}
		/* Decrease Time Constant (increase rate) for higher temperature */
		setfield Kdr_gnrh X_A->table[{i}] {taun2}
		setfield Kdr_gnrh X_B->table[{i}] {ninf}

		x = x + dx
	end
	tweaktau Kdr_gnrh X
	call Kdr_gnrh TABFILL X {tab_xfills + 1} 0
	setfield Kdr_gnrh X_A->calc_mode {NO_INTERP}
	setfield Kdr_gnrh X_B->calc_mode {NO_INTERP}
	/****************************************************************************
	** End of equations for Kdr "n" gate
	***************************************************************************/


	call Kdr_gnrh TABSAVE Kdr_gnrh.tab

   end	/* of function "make_Kdr_gnrh":  This channel creation will have to be explicitly called in setup*/

/******************************************************************************************************
**  Next channel: Fast Sodium.  Note that this is VGLB(00)'s "Hodgkin-Huxley-like" fast Na
**  Model, and in later work they moved to the Markov Chain model of Kuo and Bean.
*******************************************************************************************************/


function make_NaF_gnrh
        if (({exists NaF_gnrh}))
                return
        end

        create tabchannel NaF_gnrh

	float mpower =  3 // 3  

        setfield NaF_gnrh Ek {ENa} Gbar {GNaFs} Ik 0 Gk 0\
                Xpower {mpower} Ypower 1 Zpower 0
/******************************************************************************
** The names of the variables used here are copied from the XPP script
** "gt1v2.041.ms1b.ode"  from Artie Sherman, one of the  co-authors on LBVG(00).
********************************************************************************/

        float smna = -1		// This is s (sign) on p.8 of GT1 Channels...
/*******************************************************************************/
/* In order to adjust "Spike Threshold" to match more nearly the observed      */
/* values in the gfp-gnrh neurons, it may be useful to shift the activation of */
/* the fast sodium to higher voltages.  This will be attempted by changing the */
/* half-activation value, "vmna" in LBVG's model.  Their value (translated     */
/* into Genesis, was vmna = 0.0421.....                                        */
/*******************************************************************************/
	float vna_shift = 0.0146 //0.01976//0.0031        
	float vmna_LBVG = 0.0421
	float vmna = {vmna_LBVG} -{vna_shift}
	//float vmna = 0.0421	// This is b (V) on p.8 of GT1 Channels....
        float kmna = 0.0043	// This is c (V) on P.8 of GT1 Channels..

	float taumnass = 0.0043  // This is a (sec) on p.7 of GT1 Channels...
	float vmnat = 0.047	// This is b (V) on p.7...
	float kmnat = 0.011	// This is c (V) on p.7....
	float dmnat = 0.0001	// This is d (sec) on p.7 .....
	float rate_factor = 5	
	/* The above is a factor to compensate for LBVG(00) using rate factors measured at Room Temp, and */
	/* this model representing slices measured at 32 degrees (see discussion in De Schutter-Bower (94)) */              

	float  minf, taum,taum2, alpha, beta

	call NaF_gnrh TABCREATE X {tab_xdivs} {tab_xmin} {tab_xmax}
	x = {tab_xmin}
	dx = ({tab_xmax} - {tab_xmin})/{tab_xdivs}
	for (i = 0; i <= {tab_xdivs}; i = i + 1)

		minf = 1 / (1 + {exp { {smna}*(x+{vmna}) / {kmna} } } )
		taum = {taumnass} / ({exp { (x+{vmnat} ) / {kmnat} } } + 2*{exp {-2*(x+{vmnat}) / {kmnat} }}) + {dmnat}
		taum2 = {taum}/{rate_factor}
		/* Decrease Time Constant (increase rate) for higher temperature */
		setfield NaF_gnrh X_A->table[{i}] {taum2}
		setfield NaF_gnrh X_B->table[{i}] {minf}

		x = x + dx
	end
	tweaktau NaF_gnrh X
	call NaF_gnrh TABFILL X {tab_xfills + 1} 0
	setfield NaF_gnrh X_A->calc_mode {NO_INTERP}
	setfield NaF_gnrh X_B->calc_mode {NO_INTERP}
	/****************************************************************************
	** End of equations for NaF "m" gate
	***************************************************************************/
	
	/**********************************************************************************
	** Equations for NaF "h" gate: Again, variable names are copied from the XPP script
	** "gt1v2.041.ms1b.ode"  from Artie Sherman, one of the  co-authors on LBVG(00).
	***********************************************************************************/
 	float shna = 1		// This is s (sign) on p.8 of GT1 Channels...
       
/*******************************************************************************/
/* In order to adjust "Spike Threshold" to match more nearly the observed      */
/* values in the gfp-gnrh neurons, it may be useful to shift the inactivation of */
/* the fast sodium to higher voltages.  This will be attempted by changing the */
/* half-activation value, "vmha" in LBVG's model.  Their value (translated     */
/* into Genesis, was vmha = 0.0682.....                                        */
/*******************************************************************************/
        float vhna_LBVG = 0.0682
	float vhna = {vhna_LBVG} -{vna_shift}




        //float vhna = 0.0682	// This is b (V) on p.8 of GT1 Channels....
	//float vhna = 0.0682
        float khna = 0.0108	// This is c (V) on P.8 of GT1 Channels..

	float tauhnass = 0.150  // This is a (sec) on p.7 of GT1 Channels...
	float vhnat = 0.080	// This is b (V) on p.7...
	float khnat = 0.019	// This is c (V) on p.7....
	float dhnat = 0.0	// This is d (sec) on p.7 .....
	float rate_factor = 5	
	/* The above is a factor to compensate for LBVG(00) using rate factors measured at Room Temp, and */
	/* this model representing slices measured at 32 degrees (see discussion in De Schutter-Bower (94)) */              

	float  hinf, tauh,tauh2, alpha, beta

	call NaF_gnrh TABCREATE Y {tab_xdivs} {tab_xmin} {tab_xmax}
	x = {tab_xmin}
	dx = ({tab_xmax} - {tab_xmin})/{tab_xdivs}
	for (i = 0; i <= {tab_xdivs}; i = i + 1)

		hinf = 1 / (1 + {exp { {shna}*(x+{vhna}) / {khna} } } )
         
		tauh = {tauhnass} / ({exp { (x+{vhnat} ) / {khnat} } } + 2*{exp {-2*(x+{vhnat}) / {khnat} }}) + {dhnat}              
		tauh2 = {tauh}/{rate_factor}
		/* Decrease Time Constant (increase rate) for higher temperature */
		setfield NaF_gnrh Y_A->table[{i}] {tauh2}
		setfield NaF_gnrh Y_B->table[{i}] {hinf}

		x = x + dx
	end
	tweaktau NaF_gnrh Y
	call NaF_gnrh TABFILL Y {tab_xfills + 1} 0
	setfield NaF_gnrh Y_A->calc_mode {NO_INTERP}
	setfield NaF_gnrh Y_B->calc_mode {NO_INTERP}
	/****************************************************************************
	** End of equations for NaF "h" gate
	***************************************************************************/

	call NaF_gnrh TABSAVE NaF_gnrh.tab

   end	/* of function "make_NaF_gnrh":  This channel creation will have to be explicitly called in setup*/
/**************************************************************************************************************/

/******************************************************************************************************/
/* Set up Calcium "L" channels, after kinetics in LBVG  */

function make_CaL_gnrh
        if (({exists CaL_gnrh}))
                return
        end

        create tabchannel CaL_gnrh

	float mpower =  2 // 2  

   	

        setfield CaL_gnrh Ek {ECa} Gbar {GCaLs} Ik 0 Gk 0\
                Xpower {mpower} Ypower 0 Zpower 0
/************************************************************************
** The names of the variables used here are copied from the XPP scripts
** "gt1v2.041.ms1b.ode" and "gt1v3.031.0de" from Artie Sherman, one of the
** co-authors on LBVG(00).
*************************************************************************/

        float smcal = -1	// This is s (sign) on p.8 of GT1 Channels...
        
/*******************************************************************************/
/* In order to adjust "Spike Threshold" to match more nearly the observed      */
/* values in the gfp-gnrh neurons, it may be useful to shift the activation of */
/* the CaL to higher voltages.  This will be attempted by changing the         */
/* half-activation value, "vmcal" in LBVG's model.  Their value (translated     */
/* into Genesis, was vmcal = 0.040.....                                        */
/*******************************************************************************/
	float vmcal_shift =0.0      
	float vmcal_LBVG = 0.040   //This is b (V) on p.8 of GT1 Channels....
	float vmcal = {vmcal_LBVG} -{vmcal_shift}

        float kmcal = 0.012	// This is c (V) on P.8 of GT1 Channels..
	float tmlss = 0.005     // This is a (sec) on p.7 of GT1 Channels...
	float vtml = 0.015	// This is b (V) on p.7...
	float ktml = 0.025	// This is c (V) on p.7....
	float dtml = 0.0	// This is d (sec) on p.7 .....
	float rate_factor = 5 //5	
	/* The above is a factor to compensate for LBVG(00) using rate factors measured at Room Temp, and */
	/* this model representing slices measured at 32 degrees (see discussion in De Schutter-Bower (94)) */              

	float  minf, taum,taum2, alpha, beta

	call CaL_gnrh TABCREATE X {tab_xdivs} {tab_xmin} {tab_xmax}
	x = {tab_xmin}
	dx = ({tab_xmax} - {tab_xmin})/{tab_xdivs}
	for (i = 0; i <= {tab_xdivs}; i = i + 1)

		minf = 1 / (1 + {exp { {smcal}*(x+{vmcal}) / {kmcal} } } )
		taum = {tmlss} / ({exp { (x+{vtml} ) / {ktml} } } + {exp {-1*(x+{vtml}) / {ktml} }}) + {dtml}
		taum2 = {taum}/{rate_factor}
		/* Decrease Time Constant (increase rate) for higher temperature */
		setfield CaL_gnrh X_A->table[{i}] {taum2}
		setfield CaL_gnrh X_B->table[{i}] {minf}

		x = x + dx
	end
	tweaktau CaL_gnrh X
	call CaL_gnrh TABFILL X {tab_xfills + 1} 0
	setfield CaL_gnrh X_A->calc_mode {NO_INTERP}
	setfield CaL_gnrh X_B->calc_mode {NO_INTERP}
	/****************************************************************************
	** End of equations for CaL "n" gate
	***************************************************************************/


	call CaL_gnrh TABSAVE CaL_gnrh.tab

   end	/* of function "make_CaL_gnrh":  This channel creation will have to be explicitly called in setup*/
/***********************************************************************************************************/





/*****************************************************************/
end	/* End of Highest-level "If" loop */
















