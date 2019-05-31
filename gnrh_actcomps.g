// genesis - gnrh compartment definitions
// Steve Van Hooser 9/98 (modified from Dieter Jaeger's gpkit)
// A. Tobin 2/01 (modified from D. Jaeger's gpkit)
// C. B. Roberts 8/04 (modified from above for gnrh neurons)

/*************************************************************
* Sets up passive membrane gnrh cell compartment prototypes 
* Compartment types:
* soma   dendrite   synaptic contacts  neurite   neurax1,2,3   axon
*
***************************************************************/

function make_gnrh_comps
/* separate function so we can have local variables */

	float len, dia, surf
	int i

	echo Making gnrh active compartment library...

/******************************SOMA*****************************/

	len = 0.00e-6
	dia = 20.0e-6
	surf = dia*dia*{PI}
	if (!({exists gnrh_soma}))
		create compartment gnrh_soma
	end
	setfield gnrh_soma Cm {{CM}*surf} Ra {8.0*{RA}/(dia*{PI})}  \
	    Em {ELEAK} initVm {EREST_ACT} Rm {{RMs}/surf} \
            inject 0.0  dia {dia} len {len}

	copy NaF_gnrh gnrh_soma/NaF
	addmsg gnrh_soma gnrh_soma/NaF VOLTAGE Vm
	addmsg gnrh_soma/NaF gnrh_soma CHANNEL Gk Ek
	setfield gnrh_soma/NaF Gbar {{GNaFs}*surf}

        copy Kdr_gnrh gnrh_soma/Kdr
	addmsg gnrh_soma gnrh_soma/Kdr VOLTAGE Vm
	addmsg gnrh_soma/Kdr gnrh_soma CHANNEL Gk Ek
	setfield gnrh_soma/Kdr Gbar {{GKdrs}*surf}

        copy CaL_gnrh gnrh_soma/CaL
	addmsg gnrh_soma gnrh_soma/CaL VOLTAGE Vm
	addmsg gnrh_soma/CaL gnrh_soma CHANNEL Gk Ek
	setfield gnrh_soma/CaL Gbar {{GCaLs}*surf}



/***************************DENDRITE***************************/

    len = 200.00e-6
    dia = 2.00e-6
    surf = len*dia*{PI}
    if (!({exists gnrh_dendrite}))
        create compartment gnrh_dendrite
    end
    setfield gnrh_dendrite Cm {{CM}*surf} Ra {4.0*{RA}*len/(dia*dia*{PI})}  \
        Em {ELEAK} initVm {EREST_ACT} Rm {{RMs}/surf} inject 0.0  \
        dia {dia} len {len}
/*
	copy NaF_gnrh gnrh_dendrite/NaF
	addmsg gnrh_dendrite gnrh_dendrite/NaF VOLTAGE Vm
	addmsg gnrh_dendrite/NaF gnrh_dendrite CHANNEL Gk Ek
	setfield gnrh_dendrite/NaF Gbar {{GNaFd}*surf}

	copy Kdr_gnrh gnrh_dendrite/Kdr
	addmsg gnrh_dendrite gnrh_dendrite/Kdr VOLTAGE Vm
	addmsg gnrh_dendrite/Kdr gnrh_dendrite CHANNEL Gk Ek
	setfield gnrh_dendrite/Kdr Gbar {{GKdrd}*surf}

*/
/****************************AXON *****************************/

    len = 200.00e-6
    dia = 2.00e-6
    surf = len*dia*{PI}
    if (!({exists gnrh_axon}))
        create compartment gnrh_axon
    end
    setfield gnrh_axon Cm {{CM}*surf} Ra {4.0*{RA}*len/(dia*dia*{PI})}  \
        Em {ELEAK} initVm {EREST_ACT} Rm {{RMs}/surf} inject 0.0  \
        dia {dia} len {len}

	copy NaF_gnrh gnrh_axon/NaF
	addmsg gnrh_axon gnrh_axon/NaF VOLTAGE Vm
	addmsg gnrh_axon/NaF gnrh_axon CHANNEL Gk Ek
	setfield gnrh_axon/NaF Gbar {{GNaFa}*surf}

	copy Kdr_gnrh gnrh_axon/Kdr
	addmsg gnrh_axon gnrh_axon/Kdr VOLTAGE Vm
	addmsg gnrh_axon/Kdr gnrh_axon CHANNEL Gk Ek
	setfield gnrh_axon/Kdr Gbar {{GKdra}*surf}

        copy CaL_gnrh gnrh_axon/CaL
	addmsg gnrh_axon gnrh_axon/CaL VOLTAGE Vm
	addmsg gnrh_axon/CaL gnrh_axon CHANNEL Gk Ek
	setfield gnrh_axon/CaL Gbar {{GCaLa}*surf}





end

// actual code

pushe library
//make_gnrh_comps
pope
