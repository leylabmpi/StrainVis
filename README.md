# StrainVis: a Python-based web application for interactive visual analysis of strain-tracking methods

### Version 1.2.0

## Overview

StrainVis is a Python-based web application for visual analyses and interactive exploration of the results obtained by 
the SynTracker pipeline or by other strain tracking methods, based on ANI (Average Nucleotide Identity).  

StrainVis accepts either SynTracker's output file 'synteny_scores_per_region.csv' or an ANI file, obtained by another method, 
containing either one reference genome or multiple reference genomes (usually, one reference genome per species).
It presents accordingly analyses for each species separately and for multiple species together.  
A metadata file (matches to the samples that were previously compared by SynTracker) may also be provided
in order to enable a deeper analysis of SynTracker's results.  

StrainVis allows the user to interactively select the presented plots and to change visual parameters in each plot.
It also enables to interactively select the metadata feature by which the samples should be grouped and coloured (for each plot separately).  

Each of the presented plots can be downloaded and saved as a high-resolution image, in several common image formats.

## Installation

Download the latest release of StrainVis from: https://github.com/leylabmpi/StrainVis/releases.

Extract the tar.gz file into the desired working-directory.

Create a new conda environment using the ‘strainvis.yml’ file (located in the root directory of StrainVis) 
using the following command:
      `conda env create -f strainvis.yml`

Activate the newly created environment: 
      `conda activate StrainVis_1_2`

## Open StrainVis web-application

StrainVis is built using [Panel - a powerful data exploration & web app framework for Python](https://panel.holoviz.org/index.html).
The web-application is served using the Bokeh server. 
To launch the server from the command-line and open the application in the browser, type the following command
(from within the activated conda StrainVis environment):

`python run_strainvis.py --port 5005 &`

The application should be opened in the browser under the URL: http://localhost:5005/strain_vis .  
As long as the Bokeh server is running, StrainVis application can be accessed using the above URL.   

**Stop the Bokeh server and StrainVis application**: The two python running processes, created by the `python run_strainvis.py` command should be killed.  

### Open several StrainVis sessions to analyse different datasets simultaneously

- **Open several StrainVis browser-sessions using the same server process:**
Once StrainVis aplication was strated from the command-line, it is possible to open as many user sessions as needed 
in different browser windows/tabs, using the same URL (for example: http://localhost:5005/strain_vis). 
Each session works separately and can process a different dataset, but since all of them are executed by the same process, 
if one session is busy processing a heavy-duty task, it will affect all the other user sessions.
- **Run several server processes using different ports:** In order to process different datasets at the same time using separated computational resources, 
it is needed to start several web-server instances listening to different ports.
It simply means to run the `python run_strainvis.py --port [port]` command using a different port number for each StrainVis instance.
Each web-server instance will be accessible under: http://localhost:<port_number>/strain_vis .

## Run StrainVis on a remote server and open it in a local browser

It is possible to run StrainVis on a remote server or cloud service in cases of large datasets, in order to improve the performance.
This can be done according to the following steps:

1. **Installing StrainVis:** StrainVis should be installed on the remote server as explained above.

2. **Open SSH-tunnel:** Create an SSH tunnel from your local machine to the remote server:  
`ssh -L 5006:localhost:5006 user@remote-server`

3. **Start the Bokeh server on the remote server:**  
From within the activated conda strain_vis environment, run the following command:   
`python run_strainvis.py --port 5006 &`  

4. **Open the application in the local browser:** StrainVis should be accessible under: http://localhost:5006/strain_vis .

## Input

#### Mandatory input:

StrainVis accepts two types of input files:
- SynTracker's output file 'synteny_scores_per_region.csv' for one or multiple species.
- A tab-delimited ANI file for one or multiple species in the following format:  
  'Ref_genome', 'Sample1', 'Sample2', 'ANI'

The user should select one of the following three modes of execution, depending on the available input data:
- **SynTracker output file**: Analyse the results obtained by the SynTracker pipeline only.
- **ANI file**: Analyse the ANI results obtained by another strain-tracking method.
- **Both SynTracker and ANI files**: Analyse both types of input data (each one separately and combined).  
In this case, the names of the analysed species and the sample IDs must be identical between the two input files.
  
Note that if the input files are bigger than 300 Mb, they cannot be selected via the FileInput widget, but their full path should
be typed into the TextInput field.

#### Optional input:
A metadata file in tab-delimited format. The first column must contain the sample IDs, that should match the sample IDs
in the uploaded input file(s). The metadata file may contain an unlimited number of columns (features).

#### Sample input:
Sample input files of all three kinds are found under the 'Input_example/' directory.

## Visual analyses

When uploading an input file containing more than one species, StrainVis 
presents both single species analyses and multiple species analyses in separate tabs. 

- **Single species analyses**:
The analyses are performed on one species at a time. The species can be selected from a drop-down menu, 
containing all the species in the input file. The list of species may be sorted by the number of compared sample-pairs 
(in order to display the more abundant species first), or by the species names.

- **Multiple species analyses**:
The analyses are performed on all the species together or for a selected subset.

### Customizing the plots

Each plot allows the user to interactively set/change some visual parameters, like color, colormap, etc. 
In some of the plots it is possible to show/hide elements and to set other parameters which influence the data visualization. 

#### Including metadata in the plots

When a metadata file is uploaded, the metadata can be incorporated into the plots. The user 
can interactively select a feature, by which the presented data will be colored, grouped or filtered.  

### Saving the plots

Each plot can be saved either as an image or as a text file, containing the underlying data for the plot.  
The user may enter the name of the file (including full path), or use the default name and path provided by StrainVis 
(by default the files are saved under the 'StrainVis/Downloads/' directory).

1. **Save as image:** Each one of the plots can be saved as a high-resolution image in one of the following formats: png, pdf, svg, eps.

2. **Save data table:** The data table that was used to create the plot can be saved as text in a delimited format.

## Help pages

A detailed documentation about all types of provided analyses and customization options is found here: https://github.com/leylabmpi/StrainVis/blob/main/StrainVis_app/manual.md

The manual can also be viewed from the StrainVis web-application, under the 'Help' tab.

## Citation

If you use StrainVis please cite:

**StrainVis: interactive visual strain-level analysis of microbiome data**  
Paz I, Ley RE and Enav H.  
bioRxiv (2026). DOI: https://doi.org/10.64898/2026.03.11.711087

## Contact

If you encounter any issues or require assistance using StrainVis, please send an email to:   
inbal.paz@tuebingen.mpg.de or hagay.enav@tuebingen.mpg.de .
