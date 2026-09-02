# StrainVis: a Python-based web application for interactive visual analysis of strain-tracking methods

### Version 1.4.0

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

1. Download the latest release of StrainVis from: https://github.com/leylabmpi/StrainVis/releases.

2. Extract the tar.gz file into the desired working-directory.

3. Create a new conda environment for StrainVis by one of the following two methods:
   - **Option 1:** From the command-line, using the ‘StrainVis.yml’ file:  
     `conda env create -f strainvis.yml`  
   - **Option 2:** Double-click the setup execution file (`setup_conda_env.command` for MacOS/Linux or `setup_conda_env_windows.bat` for Windows).  
   Please note that creating the conda environment should be done once after StrainVis installation. 
   It can later be activated and used for StrainVis execution.

## Start the server and open StrainVis web-application

StrainVis web-application is built using [Panel - a powerful data exploration & web app framework for Python](https://panel.holoviz.org/index.html) 
and served in the browser using the [Bokeh server](https://docs.bokeh.org/en/latest/docs/user_guide/server.html).  
Once the server has been started, StrainVis application is accessible in the browser under the URL: http://localhost:5005/strain_vis (the port number can be set when starting the server).  
As long as the Bokeh server is running, the application can be accessed using the above URL.  

Please note that after a few hours without interacting with the application in the browser (or if your computer goes to sleep), 
it will automatically disconnect to free up memory. If a 'connection lost' message is displayed or the page stops responding, simply refresh your browser and re-upload your dataset to start a new session.

There are two available methods to launch the bokeh server and start StrainVis application:

### Option 1: Execute from the command-line (terminal)

1. Activate StrainVis environment: `conda activate StrainVis_1_4`
2. From the activated environment type: `python run_strainvis.py &`  
   **Optional arguments:**  
   `--port [PORT]`: Port to listen on (any port number that is not currently used, default is 5005)  
   `--show`: Open the application in the browser automatically (default is False)  
StrainVis web-application should be accessible in the browser under the URL: http://localhost:PORT/strain_vis


   **Stop the server:** when StrainVis is executed from the command-line, the two running python processes, created by the `python run_strainvis.py` command should be killed.

### Option 2: Use an execution file

Use one of the following wrapper scripts, which activate the conda environment, start the server and open StrainVis application in the browser automatically:

- In MacOS/Linux: double-click the `run_strainvis.command` execution file.
- In Windows: double-click the `run_strainvis_windows.bat` execution file.

**Stop the server:** closing the terminal window, that was opened by the execution file will simply kill the process and terminate StrainVis application.

### Open several StrainVis browser-sessions simultaneously

Once StrainVis application was started from the command-line, it is possible to open as many user sessions as needed 
in different browser windows/tabs, using the same URL (for example: http://localhost:5005/strain_vis). 
Each session works separately and can process a different dataset. However, since all sessions are executed by the same process,
    if one session is busy processing a heavy-duty task, it will affect the responsiveness of all the other sessions.

### Open several web-server instances using different ports

In order to process different datasets at the same time using separated computational resources, 
it is recommended to start several web-server instances listening to different ports.
It simply means to start the server for each required instance (by one of the two methods detailed above) using a different port number.
Each StrainVis instance will be accessible under: http://localhost:PORT/strain_vis .  

**Stop each web-server instance:** each StrainVis instance can be stopped by killing its two related running python processes.

### Run StrainVis on a remote server and open it in a local browser

It is possible to run StrainVis on a remote server or cloud service if your datasets are stored there.
The UI reactivity, however, is usually better in a local setup.  
Running StrainVis on a remote server can be done by the following steps:

1. **Install StrainVis:** StrainVis should be installed on the remote server as explained above.

2. **Open SSH-tunnel:** Create an SSH tunnel from your local machine to the remote server:  
`ssh -L PORT:localhost:PORT user@remote-server`

3. **Start the Bokeh server on the remote server:**  
Start the server from the command-line, as explained in the above 'Option 1'. 

4. **Open the application in the local browser:** StrainVis should be accessible under: http://localhost:PORT/strain_vis .

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
