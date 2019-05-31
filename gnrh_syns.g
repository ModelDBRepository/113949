// GENESIS script: gnrh model synaptic properties
/*******************************************************************
** This script is to set up properties for synapses in gnrh 
** neurons.  The basic script is copied from "GPsyns.g", written
** by Jesse Hanson for Globus Pallidus neurons.
*******************************************************************/
// Created by Jesse Hanson
// Modified by Carson Roberts for use in gnrh modeling

/*******************************************************************
** Each of these prototype functions will need to be explicitly 
** called in the run script (gnrh_autofit_act.g).  The constants used
** in the synapse definitions are defined in the script "gnrh_const.g"
********************************************************************/ 


function make_gnrh_AMPA
	if (!({exists AMPA}))
         create synchan AMPA
	end
	 
	setfield AMPA \
		  Ek {E_AMPA} \
                  tau1 {tauRise_AMPA} \
                  tau2 {tauFall_AMPA} \
                  gmax {G_AMPA}  \
                  frequency 0
end

function make_gnrh_GABA
	if (!({exists GABA}))
	       	create synchan GABA
	end
	setfield GABA \
	          Ek {E_GABA} \  
                  tau1 {tauRise_GABA} \ 
                  tau2 {tauFall_GABA}  \
	          gmax {G_GABA}  \
                  frequency 0
end

function make_GP_GABA_pallidum
	if (!({exists GABA_GP}))
	       	create synchan GABA_GP
	end
	setfield GABA_GP Ek {E_GABA} tau1 {tauRise_GABA_GP} \
		tau2 {tauFall_GABA_GP} gmax {G_GABA_GP} frequency 0
end

function make_gnrh_NMDA
	if (!({exists NMDA}))
                create synchan NMDA
        end
	setfield NMDA Ek {E_NMDA} tau1 {tauRise_NMDA} tau2 {tauFall_NMDA} \
		gmax {G_NMDA} frequency 0
end
	
