# SynTrackerVis: a Python-based web application for interactive visual analysis of SynTracker's results

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
