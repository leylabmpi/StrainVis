# StrainVis: a Python-based web application for interactive visual analysis of strain-tracking methods

### Version 1.0.0

## Overview

StrainVis is a Python-based web application for visual analyses and interactive exploration of the results obtained by the SynTracker pipeline or other ANI-based strain tracking methods.

StrainVis accepts a 'synteny_scores_per_region.csv' file, containing either one reference genome or multiple reference genomes (usually, one reference genome per species).
It presents accordingly analyses for each species separately and for multiple species together.  
A metadata file (matches to the samples that were previously compared by SynTracker) may also be provided
in order to enable a deeper analysis of SynTracker's results.  

StrainVis allows the user to interactively select the plots that he would like to present and to change visual parameters in each plot.
It also enables to interactively select the metadata feature by which the samples should be grouped and coloured (for each plot separately).  

Each of the presented plots can be downloaded and saved as a high-resolution image, in several common file formats.

## Installation

Download StrainVis's latest release from: https://github.com/leylabmpi/StrainVis/releases.

Extract the tar.gz file into the desired working-directory.

Create a new conda environment using the ‘StrainVis.yml’ file (located in the root directory of SynTrakcerVis) 
using the following command:  
      `conda env create -f StrainVis.yml`

Activate the newly created environment: 
      `conda activate StrainVis`

## Open StrainVis web-application

StrainVis is built using [Panel - a powerful data exploration & web app framework for Python](https://panel.holoviz.org/index.html).
The web-application is served using the Bokeh server. 
To launch the server from the command-line and open the application in the browser, type the following command
(from within the activated conda strain_vis environment):

`panel serve strain_vis.py --port 5005 --websocket-max-message-size 524288000 --show &`

The application should be opened in the browser under the URL: http://localhost:5005/strain_vis.  
As long as the Bokeh server is running, StrainVis application can be accessed using the above URL.  
Please note that several instances of StrainVis can be opened simultaneously in different browser windows/tabs.  
It is also possible to launch more than several server processes simultaneously using different ports. 

**Stop the server**: In order to stop the Bokeh server, its running process should be killed.

### Open several StrainVis sessions

- **Several user-sessions using the same web-server instance:** It is possible to open as many user sessions as needed under the same server instance that was started using the 'panel serve...' command.
All sessions will be accessible under the same URL (for example: http://localhost:5005/strain_vis). 
Each session works separately and can process a different dataset, but since all of them are executed by the same process, 
if one session is busy processing a heavy-duty task, it will affect all the other user sessions.
- **Several web-server instances:** In order to process different datasets at the same time using separated computational resources, 
it is needed to start several web-server instances listening to different ports.
It simply means to run the 'panel serve...' command using a different port number for each StrainVis instance.
Each web-server instance will be accessible under: http://localhost:<port_number>/strain_vis.


## Run StrainVis on a remote server and open it in a local browser

It is possible to run StrainVis on a remote server or cloud service in cases of large datasets, in order to improve the performance.
This can be done according to the following steps:

1. **Installing StrainVis:** StrainVis shoul be installed on the remote server as explained above.

2. **Open SSH-tunnel:** Create an SSH tunnel from your local machine to the remote server:  
`ssh -L 5006:localhost:5006 user@remote-server`

3. **Start the Bokeh server on the remote server:**  
From within the activated conda strain_vis environment, run the following command:   
`panel serve strain_vis.py --port 5006 --websocket-max-message-size 524288000 &`  
Note that it is important not to use the --show option when SynTracker runs on a remote machine.

4. **Open the application in the local browser:** StrainVis should be accessible under: http://localhost:5006/strain_vis .

## Input

#### Mandatory input:
SynTracker's output file 'synteny_scores_per_region.csv' for one or multiple species.  
Note that if the file is bigger than 300 Mb, it cannot be selected via the FileInput widget, but it's full path should
be typed into the TextInput field.

#### Optional input:
A metadata file in tab-delimited format. The first column must contain the sample IDs, that should match the sample IDs
in the uploaded SynTracker output file. The metadata may contain an unlimited number of columns (features).

## Visualization

When uploading a summary file which contains SynTracker's results for more than one species, StrainVis 
presents both Single Species Visualization and Multiple Species Visualization in separate tabs. 

#### A. Single Species Visualization
The analysis is performed for one species at a time. The species can be selected from a drop-down menu, 
containing all the species in the input file. The list of species may be sorted by the number of compared sample-pairs 
(in order to display the more abundant species first), or by the species names.

#### B. Multiple Species Visualization
The analysis is performed for all the species together or for a selected set.

### Customizing the plots

Each plot allows the user to interactively set/change some visual parameters, like color, colormap, etc. 
In some of the plots it is possible to show/hide elements and to set other parameters which influence the data visualization. 

#### Including metadata in the plots

When a metadata file is uploaded, it is possible to incorporate it into most of the APSS-based plots. The user 
can interactively select a feature from the provided metadata features, by which the presented data will be grouped and colored.  

### Saving the plots

The plots can be saved either as images or as text files, containing the underlying data of the plots.  
The user may enter the name of the file (including full path), or use the default name and path provided by StrainVis 
(under the 'StrainVis/Downloads/' directory).

1. **Save as image:** Each one of the plots can be saved as a high-resolution image in one of the following formats: png, pdf, svg, eps.

2. **Save data table:** The data table that was used to create the plot can be saved as text in a delimited format.

## Single Species Visualization

The single species visualization is combined of two types of analyses:

1. **APSS-based analyses:** Analyses based on the APSS (Average Pairwise Synteny Scores) that are calculated 
using N sub-sampled regions (according to the user's choice). 

2. **Synteny per position:** Presenting different measures based on the synteny scores for each position in the reference genome.

### *APSS-based analyses*

### Initial bar-plots
The initial presented bar-plots allow the user to interactively change the number of sub-sampled regions
and see how it influences the number of sample-pairs comparisons and the total number of samples which can be taken into account 
in the following APSS-based analyses.  
The numbers of sub-sampled regions available for selection are: 40, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350, 400. 
In case the number of compared pairs, when selecting 40 regions, is lower than 100, the option 'All regions' is also available for selection.  
Setting the desired subsampling number is done using the slider widget.
That means that all the compared pairs from all the available regions are included in the downstream analyses.

By clicking the 'Display plots using the selected number of regions' button, the following downstream analyses are being calculated and displayed.

### APSS distribution plot

This plot shows the APSS distribution among all the compared sample-pairs.  
The plot type can be changed between scatter (jitter) plot and boxplot, so as the plot color(s).  

**Including metadata:** Upon feature selection, the comparisons can be divided into two categories: same / different feature.  
for example: if the selected feature is 'country', the two categories are 'same country' and 'different country'.  
The features are derived from the upoaded metadata file and can be interactively selected from a drop-down menu.

### Clustered heatmap plot

This plot presents the APSS (Average Pairwise Synteny Scores) of the included pairwise comparisons as a clustered heatmap.
The colormap for the scores can be interactively selected from a list of available colormaps.

**Including metadata:** An additional column, coloring the rows by a requested metadata feature, 
can be added when checking the 'Use metadata for coloring' option.  
- Color rows by: select a metadata feature, by which the rows in the additional column will be colored.  
- Select colormap: select a colormap from the drop-down menu to color the different groups of the selected feature.  
When selecting the 'Define custom colormap' option, the 'Custom colormap' text input widget becomes active. 
- Custom colormap: here the user can enter a list of colors, separated by commas. 
The colors can be provided as standard names (like: red, blue) or Hex-RGB values (like: #FF0000).  
A detailed guide for color notations: https://www.w3.org/TR/css-color-4/#named-colors .


Please note that in case the number of samples exceeds the limit of 150, the heatmap plot cannot be 
displayed properly. The scoring matrix is provided for download and can be used in another visualization program.

### Network plot

This plot presents the samples as a 2D network using the Fruchterman-Reingold force-directed graph layout. 
The layout algorithm clusters the nodes (in this case, the samples) iteratively, using the APSS as the weight attribute for the network.
The higher the APSS between two samples, the closer the two nodes will be positioned in space.   
The APSS threshold (by default, the mean APSS) determines whether two samples are connected in the network or not. 
This influences the clustering process of the network, as connected nodes are being clustered closer together than unconnected nodes.  
Please note that if the number of samples exceeds the limit of 300, the network plot cannot be well displayed including 
all its interactive features. It ia possible to download the network data in tsv format and visualize it using another program.

#### Interaction with the network plot:

The network plot is created using the Bokeh backend. The Bokeh interface allows the following interactions with the plot using the mouse:  
- Panning (dragging and moving the plot)
- Zoom-in / zoom-out 
- Select a specific area and display it
- Reset to initial display
- Hover tool: displays the relevant information about each sample when hovering the network nodes. 

#### Customization options of the network (unrelated to metadata):

- **Threshold for network connections:** Enables to set the APSS threshold to define whether two nodes (samples) are 
connected in the network or not (If the APSS of a comparison between two samples is equals to or greater than the threshold, 
the samples are connected in the network). The threshold affects the clustering of the network when increasing the number of iterations.
A higher threshold means less connections between the nodes, which results in a more fragmented network, composed of a larger number of smaller clusters.  
StrainVis provides three pre-defined threshold options:   
Mean APSS (among all pairwise comparisons)  
Mean APSS + 1 std  
Mean APSS + 2 std (when the value <= 0.99)  
By selecting the 'Define another threshold' option, it is possible to set a different threshold (between 0.5 and 1.0) 
using the 'Define threshold' widget.
- **Number of iterations:** Number of clustering iterations of the network, performed using the Fruchterman-Reingold algorithm.
With each iteration, the nodes that have higher APSS score, become closer to each other in the 2D space and the network becomes more clustered.
The slider widget allows to select the number of iterations between 50 and 500.
- **Initialize nodes positions:** Clicking this button assigns the nodes new initial positions in the 2D space (by random) 
and starts the clustering process from the beginning, performing the selected number of iterations.
- **Nodes / Edges color:** Enable to set a unified color for the nodes or for the edges of the network, using a color-picker widget.
- **Show sample names:** When checked, the sample names are displayed on top of the network nodes.

#### Customization options when metadata is provided:

- **Color nodes by:** Select a metadata feature by which the nodes will be grouped and colored. 
A legend for the colored groups is added to the network plot (in case the number of groups does not exceed 10).
- **Continuous feature:** If the selected feature has numeric continuous values (for example: age), this checkbox should be checked.
The variety of colormaps available for selection changes accordingly and a colorbar is added to the network plot.  
Note that only features which numerical values can be checked as continuous.
- **Select colormap for nodes:** Select a colormap from the drop-down menu to color the nodes by the different groups of the selected metadata feature.
A different set of colormaps is provided for categorical data and for continuous data.
- **Define custom colormap:** For categorical features, it is also possible to define a custom list of colors for the different groups.
This option becomes active when selecting the 'Define custom colormap' option from the colormaps drop-down menu.
A list of colors, separated by commas, can be entered to the text-input widget. 
The colors can be provided as standard names (like: red, blue) or Hex-RGB values (like: #FF0000).  
A detailed guide for color notations: https://www.w3.org/TR/css-color-4/#named-colors .
- **Color edges by feature (same/different):** Checking this option enables to select a feature, by which the edges (connections) 
 are divided into two categories and can be colored differently. 
One category is the connections between samples that belong to the same group and the other is the connections between samples that belong to a different group (of the selected feature).
For example, if the selected feature is 'country', the edges between nodes that belong to the same country can be colored differently 
than the edges between the nodes that belong to a different country.
- **Color edges by:** Select the feature, by which the edges are divided into two categories.
- **Same color / Different color:** Select the color that will be applied on each category of edges.

#### Saving the network

- **Save as image:** When using the 'Plot download options' and selecting one of the four available image formats, 
the network graph area is saved in its initial presentation (without considering interactions such as zoom-in/out, etc.).  
In order to save the exact current view (including graph interactive modifications), it is possible to use the Bokeh interface 'Save' button,
 which exports the graph in png format only (it is then saved in the browser's default Downloads location).
- **Save data table:** The network data is saved as a tab-delimited format file, containing the following columns: Sample1, Sample2, APSS, weight.  
The weight column reflects the score of each pair of samples after applying the APSS connections threshold. 

### *Synteny per position analyses*

These analyses are presented for each contig of the selected reference genome. If there is more than one contig, 
it can be selected from the list of contigs using a drop-down menu.
By default, the contigs are sorted by their length, but they can be sorted by their names as well.  

The synteny per position plot displays the following four analyses on top of each other, 
where each analysis can be marked as shown or hidden in the plot. 

1. **Average synteny scores per region:** This line-plot shows the average synteny score for each region in the reference genome.
It is calculated based on all the sample-pairs comparisons derived from a specific region.  

2. **All synteny scores per region:** This horizontal line plot shows all the synteny scores for each region on the reference genome.  

3. **Hypervariable regions:** Highlight regions that meet the following criteria:  
They are among the bottom 10% regions with the lowest average synteny scores.  
They appear in at least 10% of the compared sample-pairs.

4. **Hyperconserved regions:** Highlight regions that meet the following criteria:  
They are among the top 10% regions with the highest average synteny scores.  
They appear in at least 50% of the compared sample-pairs.

#### Customization options of (unrelated to metadata):

- **Set contig length range:** Set the start and end positions of the contig that will be presented in the plot.
- **Set new range button:** Clicking this button updates the plot to present the newly set range.
- **Reset range button:** Clicking this button resets the range to show the whole contig in the plot.
- **Colors:** The color of each one of the plots can be set separately, as well as the alpha transparency of the highlighted conserved / variable regions.

#### Customization options when metadata is provided:

- **Filter plot by metadata:** The metadata feature, by which the presented data will be filtered, can be selected from the drop-down menu.
- **Include the following groups in the plot:** It is possible to select one or more groups to be included in the plot.
- **Filter plot button:** Clicking this button updates the plot, so that only pairwise comparisons, originating from the selected groups of the selected feature, will be included in the plot.
- **Reset filteration button:** Clicking this button resets the filtering and updates the plot so that all data is shown.

## Multiple Species Visualization

The multiple Species visualisation tab is active when the input file contains more than one species. 
It presents APSS-based analysis for all the species or for a selected subset.

### Setting the species that will be included in the analysis
- **All species:** Include all available species.
- **Select a subset of species:** Using the multi-select widget, it is possible to select specific species to be included in the analysis.
By default, the species are sorted by their abundance (the number of compared pairs), but they can be sorted alphabetically by their names as well.
- **Update species selection button:** Clicking this button updates the set of species that are included in the analysis.

### Initial bar-plots
The initial presented bar-plots allow the user to interactively change the number of sub-sampled regions
and see how it influences the number of compared sample-pairs and the total number of species that can be taken into account 
in the following APSS-based analysis.  
The numbers of sub-sampled regions available for selection are: 40, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350, 400. 
In case the number of compared pairs, when selecting 40 regions, is lower than 100, the option 'All regions' is also available for selection.
That means that all the compared pairs from all the available regions are included in the downstream analysis.  
Setting the desired subsampling number is done using the slider widget.

By clicking the 'Display plots using the selected number of regions' button, the following downstream analyses are being calculated and displayed.

### APSS distribution among species plot

This plot shows the APSS distribution among all the compared sample-pairs of each included species as a boxplot.  
The plot color can be changed using the color-picker widget.  

**Including metadata:** When checking the 'Use metadata in plot' checkbox, it is possible to select a feature, 
by which the comparisons can be divided into two categories: same / different.  
for example: if the selected feature is 'country', the two categories are 'same country' and 'different country'.  
The features are derived from the uploaded metadata file and can be interactively selected from a drop-down menu.  
The colors of the same / different feature categories can be changed using the color-picker widgets.  
When using metadata, the P-values of the comparisons for the selected feature can be downloaded in addition to the APSS table.
