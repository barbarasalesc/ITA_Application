# About the scripts on this directory

## Overview 
These scripts were developed as a hands-on educational activity for a Dark Matter School. The primary objective was to provide students with the necessary tools and resources to independently perform a Pulse Shape Discrimination (PSD) optimization analysis on events recorded by the PSci-arrays of the ReD experiment. I cannot provide the data used to this analysis, since it is proprietary to the collaboration.

Below is a description of each script and its role in the activity:

- create_file_traces.C: Extracts a custom "toy dataset" from the full ReD 252-Cf dataset. It filters and outputs only the relevant information required for the students to perform the PSD analysis.

- traceprofile.C: An analytical solution that computes the mean trace profile of the events. It allows students to observe electronic phenomena, such as signal reflection.

- trace-view.C: A visualization tool that plots the trace of an individual event, allowing students to inspect any specific event of their choice in detail.

- fprompt-analizer.C: Performs the core PSD optimization analysis. It uses a Figure of Merit (FoM) method to optimize the short and long charge integration time windows, providing the best discrimination between photon and neutron events.
