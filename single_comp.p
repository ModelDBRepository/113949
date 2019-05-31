// genesis
//
// Cell morphology file for GENESIS.
// Written by cvapp (http://www.neuro.soton.ac.uk/cells/#software).

*absolute
*asymmetric
*cartesian


// End of cvapp-generated header file.

// Settings suggested by D. Jaeger and AE Tobin to make sure that proper parameters
// are passed to all compartments
*set_compt_param RM {RMs} // This is the one sent in from the GA
*set_compt_param RA {RA}
*set_compt_param CM {CM} // F/m^2
*set_compt_param ELEAK {ELEAK} // Volts
*set_compt_param EREST_ACT {EREST_ACT} // Volt

*compt /library/gnrh_soma
//soma  none    0  0  0  16.169
soma none 0 0 0 21
//soma  none    0  0  0  1.063 


