# SynTrackerVis: a Python-based web application for interactive visual analysis of SynTracker's results

### Version 1.0.6

## Overview

SynTrackerVis is a Python-based web application, designated for visual analyses and interactive exploration of the results obtained by the SynTracker pipeline.  

SynTrackerVis accepts a 'synteny_scores_per_region.csv' file, containing either one reference genome or multiple reference genomes.
It presents accordingly analyses for each genome separately and for multiple genomes together.  
A metadata file (matches to the samples that were previously compared by SynTracker) may also be provided
in order to enable a deeper analaysis of SynTracker's redults.  

SynTrackerVis allows the user to interactively select the plots that he would like to present and to change visual parameters in each plot.
It also enables to interactively select the metadata feature by which the samples should be grouped and coloured (for each plot separately).  

Each of the presented plots can be downloaded and saved as a high-resolution image, in several common file formats.

## Installation

Download SynTrakcerVis's latest release from: https://github.com/leylabmpi/SynTrackerVis/releases.

Extract the tar.gz file into the desired working-directory.

Create a new conda environment using the ‘SynTrackerVis.yml’ file (located in the root directory of SynTrakcerVis) 
using the following command:
      `conda env create -f SynTrackerVis.yml`

Activate the newly created environment: 
      `conda activate SynTracker_Vis`

## Open SynTrackerVis web-application

SynTrackerVis is built using [Panel - a powerful data exploration & web app framework for Python](https://panel.holoviz.org/index.html).
The web-application is served using the Bokeh server. 
To launch the server from the command-line and open the application in the browser, type the following command
(from within the activated conda SynTracker_Vis environment):

`panel serve syntracker_vis.py --port 5005 --websocket-max-message-size 524288000 --show &`

The appplication should be opened in the browser under the URL: http://localhost:5005/syntracker_vis.  
As long as the Bokeh server is running, SynTrackerVis application can be accessed using the above URL.  
Please note that several instances of SynTrackerVis can be opened simultaneously in different browser windows/tabs.  
It is also possible to launch more than several server processes simultaneously using different ports. 

**Stop the server**: In order to stop the Bokeh server, its running process should be killed.

## Run SynTrackerVis on a remote server and open it in the local browser

It is possible to run SynTrackerVis on a remote server or cloud service in cases of large datasets, in order to improve the performance.
This can be done according to the following steps:

1. **Installing SynTrackerVis:** SynTrackerVis shoul be installed on the remote server as explained above.

2. **Open SSH-tunnel:** Create an SSH tunnel from your local machine to the remote server:  
`ssh -L 5006:localhost:5006 user@remote-server`

3. **Start the Bokeh server on the remote server:**  
From within the activated conda SynTracker_Vis environment, run the following command:   
`panel serve syntracker_vis.py --port 5006 --websocket-max-message-size 524288000 &`  
Note that it is important not to use the --show option when SynTracker runs on a remote machine.

4. **Open the application in the local browser:** SynTrackerVis should be accessible under: http://localhost:5006/syntracker_vis .

## Input

#### Mandatory input:
SynTracker's output file 'synteny_scores_per_region.csv' for one or multiple genomes.  
Note that if the file is bigger than 300 Mb, it cannot be selected via the FileInput widget, but it's full path should
be typed into the TextInput field.

#### Optional input:
A metadata file in tab-delimited format. The first column must contain the sample IDs, that should match the sample IDs
in the uploaded SynTracker output file. The metadata may contain an unlimited number of columns (features).

## Visualization

When uploading a summary file which contains SynTracker's results for more than one reference genome, SynTrackerVis 
presents both Single Genome Visualization and Multiple Genomes Visualization in separate tabs. 

#### A. Single Genome Visualization
The analysis is performed for one reference genome at a time. The reference genome can be selected from a drop-down menu, 
containing all the genomes in the input file. The list of genomes may be sorted by the number of compared sample-pairs 
(in order to display the more abundant species first), or by the genome names.

#### B. Multiple Genomes Visualization
The analysis is performed for all the reference genomes together or for a selected set.

### Customizing the plots

Each plot allows the user to interactively set/change some visual parameters, like color, colormap, etc. 
In some of the plots it is possible to show/hide elements and to set other parameters which influence the data visualization. 

#### Including metadata in the plots

When a metadata file is uploaded, it is possible to incorporate it into most of the APSS-based plots. The user 
can interactively select a feature from the provided metadata features, by which the presented data will be grouped and colored.  

### Saving the plots

The plots can be saved either as images or as text files, containing the underlying data of the plots.  
The user may enter the name of the file (including full path), or use the default name and path provided by SynTrackerVis 
(under the 'SynTrackerVis/Downloads/' directory).

1. **Save as image:** Each one of the plots can be saved as a high-resolution image in one of the following formats: png, pdf, svg, eps.

2. **Save data table:** The data table that was used to create the plot can be saved as text in a delimited format.

## Help pages

A detailed documentation about all the types of provided analyses and the visualisation options is found here: https://github.com/leylabmpi/SynTrackerVis/blob/main/SynTrackerVis_app/manual.md

The manual can also be viewed from the SynTrackerVis web-application, under the 'Help' tab.
