In the ae4316P_20XX_data_group*.mat file you will find:

%% DATA STRUCTURES
Data for each subject is contained in a data structure with the following name format:

group<group name>_data_subj<subject number>  

The group name is either "NM" for the group who received pre-training in a fixed-base
simulator, or "M" for the control group. Note that you collected data for subject 1
in group M. 

In each of these structs you will find both the time histories and the estimated
pilot frequency responses for the 60 repeated measurements that were collected.
You can access the measured data in each Matlab data structure with the "." operator!

%% TIME HISTORIES %%
- the time vector t
- arrays that hold the 60 replications of the time histories of:
	- the error signal              e [deg]
	- the control signal            u [deg]	
	- the controlled pitch attitude x	[deg]

%% MEASURED PILOT FREQUENCY RESPONSES %%
- the frequency vector w in rad/s (consists of the frequencies of the twenty sinusoids in the disturbance forcing function signal)
- arrays that and for the 60 replications hold the measured frequency response:
	- the measured frequency response magnitude: magHp   [-]
	- the measured frequency response phase: 	   phaseHp [deg]

Good luck with the data analysis! 
