# TLE_fingerprint
Computation of TLE connectivity profiles and their distance to TLE functional and structural MRI fingerprint. Plots profile for the individual patient.

Morgan VL, Sainburg LE, Johnson GW, Janson A, Levine KK, Rogers BP, Chang C, Englot DJ. Presurgical temporal lobe epilepsy connectome fingerprint for seizure outcome prediction. Brain Communications 2022, in press.

__________________________________________________________________________________
TLE_fingerprint_distance.m

[MEFC_pat, MESC_pat, MEFCSC_pat, pat_profile] = TLE_fingerprint_distance(model,pat_data, side_1)

___________________________________________________________________________________
This function will plot patient profiles and the TLE fingerprint on a radar plot as shown in Morgan et al. Brain Communications 2022 Figure 5 subplots

INPUTS – 
model = matrix of the TLE fingerprint in 40 x 6 x 2
	40 are 40 permutations of seizure free patients
	6 are lobes = prefrontal, parietal, occipital, temporal,
	motor/somatosensory and subcortical
	2 are FC and SC

pat_data = matrix of the invidual patient data in n xn n x 2 = nodes x nodes x FC and SC

side_1 = side of focus - left = 0, right = 1

this program loads in the excel spreadsheet - node_info.csv which contains info about the nodes, their order, their side and their lobe. Interest column should be one for the nodes of interest.

pat_data(:,:,1)=all_subj_FC(36,:,:);
pat_data(:,:,2)=all_subj_SC(36,:,:);
model(:,:,1)=w_perm_FC_lobe(:,1:6,1);
model(:,:,2)=w_perm_FC_lobe(:,1:6,2);

OUTPUTS –
MEFC = functional connectivity distance to model profile
MESC = structural connectivity distance to model profile
total = total distance to model profile
pat_profile = 6 lobes x 2 (FC and SC) that can be plotted with plot_profile.m 
___________________________________________________________________________________________

plot_profile.m

plot_profile(model,pat_profile)
_______________________________________________________________________________________________

INPUTS - 
model = matrix of the TLE fingerprint in 40 x 6 x 2
	40 are 40 permutations of seizure free patients
	6 are lobes = prefrontal, parietal, occipital, temporal,
	motor/somatosensory and subcortical
	2 types of connectivy are FC and SC

pat_profile = matrix of the invidual patient profile output from TLE_fingerprint_distance.m matrix in 6 x 2
	6 are lobes = prefrontal, parietal, occipital, temporal, motor/somatosensory and subcortical
	2 types of connectivity are FC and SC

___________________________________________________________________________________________

Other files – 
Model_data.mat – matlab data with a sample patient (pat_data) and the model data (model) for use in functions above.

node_info.csv – information on the 117 regions of interest used for this research used to run the TLE_fingerprint_distance.m. File includes label names, side of brain, lobe of brain and atlas values.

For creating your own patient data:
1.	You need data in the regions in node_info.csv (or edit the csv file)
2.	nwT1W3D_0001.nii is a template image in MNI space
3.	MultiAtlas_subj.nii is an atlas of the regions in the node_info.csv that is computed from the nwT1W3D_0001.nii. This was created using MultiAtlas segmentation 
	https://github.com/VUIIS/Multi-Atlas-v3.0.0
4.	Connectomes must be in units of standard deviations from healthy controls. 

To do age correction of input patient data:
1. 	Start with patient data in regions ordered by node_info ID
2. 	FC should be in correlation r values
	SC should be in mrtrix3 stramline count values
3.	need age of patient in years and length of time series for FC
4.	run FCR_to_ageZ.m to convert FC r values to Z and then age correct resulting in std from age matched control
5.	run SC_to_age_log.m to convert SC to log values and then age correct resulting in std from age matched control 
6.	can test these by loading before_age_corr.mat and checking the first outputs with after_age_corr.mat	

