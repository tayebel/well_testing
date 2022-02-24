# well_testing
1. Program goal:
we set up a series of mock well tests, with a well opening and closing, and then simulate the evolution of the transient flow and pressure over time. This is useful when we do for the next step of the studies a well test analysis to identify the storage in the well, the finite-action and infinite-action flow time, and most importantly to derive the well permeability and skin factor from the well test data.

we will study the basic equations that describe the transient single-phase fluid flow in porous media and discuss the general idea of modeling this flow using python. 

We will include the various equations and conditions covered by this research in order to plot the flow and pressure profiles that are so beneficial in obtaining and identifying well and reservoir behavior and parameters.

2.	The insertion of functions used in the program design

2.1.	 The dimensionless parameters:
	
  We will begin by creating functions for the adimensional parameters (rD, tD, PD and QD) [Equation (4.3), (4.5), (4.8), (4.9), (4.22), (4.23)] as well as the condition where the flow begins to behave infinitely.

2.2.	 The solutions in variable flows and pressures: 

  The next step is to create a function defining the pressure and flow in a multi-flow and multi-pressure system respectively. [Equation (4.24), (4.26)]
  
2.3.Simulation of the multiple constant flow and multiple constant pressure test based on the superposition principle:  

	Finally, we will create the two functions ("simulate_multitirate_test" and "simulate_multipressure_test") that allow us to calculate q,t and Pwf and to plot their variations as a function of time.

2. Running  the following axample:
[Screenshot 2022-02-24 153424](https://user-images.githubusercontent.com/100303910/155544690-fddf1a1f-b72f-4d63-bf67-6759b2c7dc32.png)

![example _table](https://user-images.githubusercontent.com/100303910/155551762-b7662bbf-6df8-4c90-b9f7-b1dd5b5607d8.png)

 After running the example , we obtained promising results:
 
 ![result2](https://user-images.githubusercontent.com/100303910/155554125-9db511c4-dce8-43b7-ad38-642172bf8830.png)

![result](https://user-images.githubusercontent.com/100303910/155554104-ca2be24b-3780-4a7a-bcda-499ad5d1a8ec.png)


 
 
 
 
