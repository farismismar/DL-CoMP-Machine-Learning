For this code to work, please first obtain access to Vienna LTE-A Downlink System Level Simulator from https://www.nt.tuwien.ac.at/research/mobile-communications/vccs/vienna-lte-a-simulators/lte-a-downlink-system-level-simulator/

Also, the code uses `matlab2tikz` to plot figures, which can be separately downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz

Reproducibility is not guaranteed unfortunately if the data_files/channel_traces/\*.mat files are different or invalid.  These files are made by the simulator and not from my codes.

An important file LTE_sim_main.m is missing since it is part of the Vienna LTE-A simulator. I had to modify this file to implement the various CoMP algorithms, but of course cannot share the file due to copyrights.


## Version history
2/24/2019 Initial code release

7/30/2019 Version 2.  Fixed a few bugs and added code for both classifiers.
