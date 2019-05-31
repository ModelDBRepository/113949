//genesis - gnrh cell  genesis2 script

function do_current_injections

/* do a negative current injection to inhibit spiking */
setfield /pulse level1 0
//setfield /pulse level1 +0.0003e-9 //  For Single Spike Sim
//setfield /pulse level1 -0.0008e-9  // For Double Spike Sim	 
//setfield /pulse level1 -0.0037e-9   // For Latch-up Sim 
step 0.4 -time

/* do a short, strong positive current injection (hope to elicit ONE spike) */
//setfield /pulse level1 0.33e-9       //For Single Compartment Model
setfield /pulse level1 0.46e-9	       //For Full model
step 0.0005 -time

/* do a negative current injection to inhibit spiking */
//setfield /pulse level1 0	
//step 0.1 -time

/* do a positive current injection  (hope to elicit TWO spikes)*/
//setfield /pulse level1 0.0165e-9	  // For Single compartment
//setfield /pulse level1 0.027e-9 //0.018e-9           // For Full Model
//step 0.05 -time

/* do a negative current injection to inhibit spiking  */
setfield /pulse level1 0
//setfield /pulse level1 +0.0001e-9        // For Single Spike sim
//setfield /pulse level1 -0.0004e-9      // For Double Spike Sim	
  step 0.8 -time

/***********  Additional Current Injections *************************/



/* do a long, strong positive current injection (hope to elicit ONE spike and latch-up) */	
//setfield /pulse level1 2e-9	
//step 0.2 -time

/* do a negative current injection to inhibit spiking  */
//setfield /pulse level1 +0.0001e-9        // For Single Spike sim
//setfield /pulse level1 -0.0004e-9      // For Double Spike Sim
//setfield /pulse level1 -0.001e-9	 // For Bigspike Sim
//  step 0.8 -time

//setfield /pulse level1 0
//step 0.05 -time

//setfield /pulse level1 0.0169e-9
//step 0.05 -time

setfield /pulse level1 0
step 0.1 -time


float curr = 0.0e-9
 while ({curr}<0.45e-9)
  echo  
  echo injected current  {{curr}*1e9} nA				     				     
  setfield /pulse level1 {curr}
  step 0.05 -time
  curr={curr}+0.005e-9
end


float curr = 0.45e-9
 while ({curr}<2e-9)  
  step 0.05 -time 
  echo injected current  {{curr}*1e9} nA
  echo
  setfield /pulse level1 {curr}				     
  curr={curr}+0.05e-9  
end


end     // of function do_current_injections
