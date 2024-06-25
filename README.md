# Models 
Hi there. Welcome to my documents. Here, I have placed my Python scripts on Residual Strength. 
These scripts are based on the methodologies described in:

(1)   Mourik, Prediction of the erosion velocity of a slope of clay due to wave attack, (2020)
(2)   Klein Breteler, Erosie van kleibekleding met gras op boventalud van Waddenzeedijken (2022)

What can be changed, when applying these models to other sites, are the following:

For (1)   Mourik, Prediction of the erosion velocity of a slope of clay due to wave attack, (2020)

A. Hydraulic conditions:
  line 12,13
  Duration of the storm -> default set at 35 hours 
  Water level min./max. -> default set at -0.1 and 2.72 m+NAP, respectively
  Wave height min./max. -> default set at 0.1 and 2.65 meters, respectively
  Peak period max. -> default set at 5.75 seconds

B. Dike related constants:
   line 64
   Ce = erosion coefficient -> default set at 0.55
   alpha = slope of original profile
   s_op = wave steepness -> default set at 0.05 (-)
   

For (2)   Klein Breteler, Erosie van kleibekleding met gras op boventalud van Waddenzeedijken (2022)

A. Hydraulic conditions:
  Same as for (1) Mourik above.
  
B. Dike related constants
  line 47-55
  m2 = model constant phase 2 -> default set at 1.0 (-)
  s_op = wave steepness -> default set at 0.05 (-)
  h_berm = height point where berm intersects upper slope -> default set at 2.52 m+NAP
  B_berm = width berm -> default set at 3.0 m
  tan_alpha = slope of  -> default set at 0.333
  a1 = original upper slope angle -> default set at 1.107 (in rad.)
  a2 = terrace angle erosion profile -> default set at 0.099 (in rad.)
  a3 = cliff angle erosion profile -> default set at 0.321 (in rad.)

The motivation behind the default values, the explanation of the different variables and further documentation is placed in the scripts itself, and in the master thesis "Incorporating residual strength regarding clay resistance in dike safety assessment - a case study of the Ketelmeerdijk" of De Jong (2024)

  
 

