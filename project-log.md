# WEEKLY LOG OF EP3051 PROJECT ICECUBE DATA ANALYSIS

Last edited 2506

Week 1: 0305 - 1005

Read few papers on neutrinos:

	https://arxiv.org/abs/2202.00694
	
	https://arxiv.org/abs/2112.13820
	
Data used for the project was downloaded from:
	https://www.atnf.csiro.au/research/pulsar/psrcat/ for ms pulsars
	
	https://icecube.wisc.edu/data-releases/2021/01/all-sky-point-source-icecube-data-years-2008-2018/ for icecube neutrinos
	
Read introductory articles on MCMCs and tried using emcee
	
Week 2: 1005 - 1705
Read a paper on neutrino physics https://arxiv.org/pdf/1910.11878.pdf

Installed and configured Linux on advise of Shantanu sir

Learnt about Monte Carlo simulations and MCMC's by writing few simple practice codes but yet to succsfully use emcee

Watched an Icecube colloqium https://www.youtube.com/watch?v=J5S5wN4Kfz4

Week 3: 1705 - 2405

	Assigned 1st task: a) plotting the sky distribution of IceCube neutrinos, 
			   b) energy spectrum of IceCube neutrinos as well as 
			   c) sky distribution of millisecond pulsars
	
Learnt astropy and it's required applications

Created github repo for the project

	Completed 1a:- Plotted the sky distribuutions in all biultin projections but focused on Hammer projection.

	1b incorrect attempt due to mis-interpretation of energy spectrum
	1c not done as psrcat files could not be compiled (later advised not needed)
	
Week 4: 2405-3105

	1c Plotted ms pulsars with Galactic latitude and logitude
	
	1b not complete as still did not understand energy spectrum meant histogram of energy so 
		advised to get an output similar to Fig 5 in https://arxiv.org/pdf/0711.0053.pdf
		Still could not get it.
		
Week 5: 3105 - 0706
	1b 
		Still plotting incorrect graphs by using formulas from https://arxiv.org/pdf/1910.11878.pdf to get N as a function of E
		Advised to plot hist(E) in matplotlib binned logarithmically in energy
		0606 - Finally plotted the graph of E binned in logE.
		       Used the logspace from log10(2.6) to log10(7) with 50 bin
	
Completed Task 1 on 06062022

Assigned Task 2/2a on 06062022:
	
	2a) Calculate the space angle between each neutrino detected by Icecube and all the millisecond pulsars.
	Then bin the data as a function of angles between 0 and 30 degrees (with bins of 5 degree each).
	
Cautioned against practical issues

	1. Need to consider the full IceCube data (not one year at a time as what you have done)
	
	2. Because of (1), calculating the angle for all neutrino-millisecond pulsar pairs could take a long time if done a brute force search.
	Advised to either vectorize in python or use a low-level language (such as C)

Week 6: 0706 - 1406

0906 - Completed Task 2/2a,  
0906 - Assigned task 2b:	
	Calculate all matches between IceCube neutrinos and millisecond pulsars within 3 degrees
		(This will constitute the total signal events.)
		
	To estimate the background, bin the data between 0 and 7.35 degrees into 6 equispaced cos (theta) bins.
		i.e, bin the  data into six equal cos(theta) bins between cos(0) and cos(7.35 degrees).
		(The bin closest to 0 corresponds to 3 degrees)
		The average of the last 5 bins will constitute the background.
	


1006 - 1106 - Completed task 2b
	      Rectified a mistake in Task 2a and updated it

1106 - Assigned Task 2c:
	
	Count the number of msp-neutrino matches within each neutrino error circle
	(Icecube provides the angular resolution / angerr for each neutrino which varies from neutrino to neutrino)
		I.e, if neutrino 1 has a resolution of x degrees, neutrino 2 y degrees, neutrino 3 z degrees,
		then calculate all MSP matches within x degrees of neutrino 1 + all MSP matches within y degrees of neutrino 2 
		    + all MSP matches within z degrees of neutrino 3 and so on
		
		
Week 7: 1406 - 2106

	Task 2c required high memory usage and took long to compile. So tried using GPU (PyCUDA, but was unsuccessful)

1906 - Completed Task2c (CPU execution) got 44972 matches
		Rectified a bug in Task 2a and 2b which prevented the code from iterating throughout the data. Waiting for next task
		
2006 - Assigned task 2d:
	Estimate the background by finding the number of coincidences within the same solid angle
	    corresponding to the error region but offset by 5 degrees.
	
	I.e for a given neutrino that has an error of X degrees, an event is considered background if it satisfies
	
	$|cos(\theta)-cos(5)|  < 0.5 x (1-cos X)$
				 (theta is the angle between pulsar and neutrino event)
				 
Week 8: 2106 - today 2506:

	Unsuccesfully tried GPU coding.

	BREAK 2106 - 2206:	Reinstalled OS as it crashed while reinstalling python

	2306 - First incorrect attempt: got 53101560 matches.
		Advised to re check
	2306 - Found the error: Did not take Cos(X) in the formula above. This resulted in abnormally large no.of matches as X is in degrees so 1-X frequently was < 0.
	
	2306 - Successfully completed task 2d/2c2 - got 31091 matches
	2506 - Updated Task 1a y plotting the ENTIRE sky distribution from 2008-18 in HAMMER projections using (RA, DEC) and (Gal_lat, Gal_lon)
	Waiting for next task
		
	
		



	

