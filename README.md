# PhD-thesis
This repository contains the code files for the unpublished content of my PhD thesis.

FILE LIST
(Note on file extensions: py and pyx = python; R and rmd = R)

Chapter 2: plot and function files are located in the Chapter 4 folder; files from the publication are located in the dedicated repository
relationship_encountered_stones_to_depositions_sim_data.R
	numerical simulation for the discussion of Figure 3 in Franks and Deneubourg (1997)

run_large_sim
	simulations with all models (template-only, template + feedback, gradual model) across stone density values (runs sims, extracts results, and makes figures; funcs are in wall_funcs.py)

SA.pyx
	runs and analyses results from the sensitivity analysis (all models)


Chapter 3
analyseData.R
	combines brood data and building event data files and calculates per-min building rate

Activity analysis - overall variables - best models selection.rmd (output: Activity-analysis-by-ROI---best-models-selection.pdf):
	regression-model evaluation and comparison

best_starting_values_HMMs.R
	fits HMM regressions across a range of starting parameter values to ensure best fit is the global optimum

rates.csv
	processed data (used for model fit) 
	Column legend:
		D = per-min deposition rate calculated over a 15 min interval
		stone_dens = stone_dens value (averaged across interval)
		n_ants =  n of ants in the brood cluster (averaged across interval)
		distance = distance between building site and brood cluster centroid, averaged across buildig sites (averaged across interval)
		timepoints = time since start of building activity in the colony (min)
		ID = colony ID

for data-processing replication:
	Deposition_data-colonies_R54-R34-R29-R5.csv: raw building event data
		Colony = colony ID
		ROI = building site (numbered 1-3 in each colony)
		x,y = building site coordinates in video (in pixel)
		event = type of building event
		behaviour_type = subtype of building event
		time_in_video = time from video start at which the event takes place (hh:mm:ss format)
		cut_video_start = lenght of cut from original video
		time =  time from original video start at which the event takes place
		time_s, cut_from_video_start_s, cut_video_start_s = convertion to seconds of previous three columns
		n_stones = n of stones at building site when event starts (i.e., not including effect of current event or of other events occurring simultaneously)
		pixel/mm = n of pixels per mm in video
		radius_mm = building site radius converted to mm
		stone_density = n stones/mm at building site when event starts

	brood_measurements - colonies_R34-R54-R29-R5 - anonymous measurers.xlsx: brood measurement data taken at 15 min time intervals
		standardised_n_ants_in_cluster_1.5-ant = n of ants in cluster standardised within colony
		X_centroid_1.5-ant, Y_centroid_1.5-ant = coordinates of brood cluster centroid (in pixels)
		cluster_density_1.5-ant = density of brood cluster (ants/mm^2)
		measurer = measurer ID

	all_building_data - colonies_R54-R34-R29-R5.csv: combines above files
		normalised time =  time of event (s) from first building event in colony

	ROI_coordinates.xlsx
	pixel coordinates of the centre of each ROI


Chapter 4
ind_var_sim.pyx
	epigenetic individual variation file (runs sims and extracts results; funcs are in wall_funcs.py)

make_plots.R
	makes all plots from building simulations (Chapter 2 and Chapter 4) 

patrilines_sim.pyx 
	patriline-specific individual variation file (runs sims and extracts results; funcs are in wall_funcs.py)

Rq_sweep.pyx
	mutation-in-germline file (runs sims and extracts results; funcs are in wall_funcs.py)

wall_funcs.py
	contains all wall-building simulation functions
