import gc
import panel as pn
import pandas as pd
import numpy as np
import io
import os
import re
from functools import partial
import holoviews as hv
import networkx as nx
import random
from bokeh.io import export_svgs, export_png
import SynTrackerVis_app.config as config
import SynTrackerVis_app.data_manipulation_single as ds
import SynTrackerVis_app.data_manipulation_multi as dm
import SynTrackerVis_app.plots_single_genome as ps
import SynTrackerVis_app.plots_multi_genomes as pm
import SynTrackerVis_app.widgets as widgets
pn.extension('floatpanel')
#pn.extension(disconnect_notification='Connection lost, try reloading the page!')


def change_disabled_state_inverse(chkbox_state):
    if chkbox_state:
        return False
    else:
        return True


def change_disabled_state_straight(chkbox_state):
    if chkbox_state:
        return True
    else:
        return False


def change_continuous_state(chkbox_state):
    if chkbox_state:
        return config.Bokeh_continuous_colormap_dict
    else:
        return config.Bokeh_categorical_colormap_dict


def change_continuous_state_for_value(chkbox_state):
    if chkbox_state:
        return config.Bokeh_continuous_colormap_dict['Turbo256']
    else:
        return config.Bokeh_categorical_colormap_dict['Category20']


def enable_submit(file_input, text_input):
    if file_input is not None or text_input != "":
        return False
    else:
        return True


def generate_rand_pos():
    rand = random.random() * 2 - 1
    return rand


class SynTrackerVisApp:

    def __init__(self):

        # Variables
        self.is_file_path = 0
        self.filename = ""
        self.input_file_loaded = 0
        self.is_metadata = 0
        self.metadata_dict = dict()
        self.metadata_features_list = []
        self.sample_ids_column_name = ""
        self.genomes_select = ""
        self.number_of_genomes = 0
        self.ref_genomes_list = []
        self.ref_genome = ""
        self.sampling_size = ""
        self.score_per_region_all_genomes_df = pd.DataFrame()
        self.score_per_region_selected_genome_df = pd.DataFrame()
        self.avg_score_per_region_selected_genome_df = pd.DataFrame()
        self.df_for_network = pd.DataFrame()
        self.avg_scores_all_genomes_all_sizes_dict = dict()
        self.calculated_avg_genome_size_dict = dict()
        self.working_directory = os.getcwd()
        self.downloads_dir_path = self.working_directory + config.downloads_dir
        self.contigs_dict = dict()
        self.contigs_list_by_name = []
        self.contigs_list_by_length = []
        self.avg_score_genome = 0
        self.std_score_genome = 0

        # Bootstrap template
        self.template = pn.template.VanillaTemplate(
            title='SynTrackerVis',
            site_url="syntracker_vis"
        )
        self.main_area = pn.Column(styles=config.main_area_style)
        self.template.main.append(self.main_area)

        # Widgets
        self.text_input = pn.widgets.TextInput(name='', placeholder='Enter the file path here...')
        self.input_file = pn.widgets.FileInput(accept='.csv, .tab, .txt')
        self.metadata_file = pn.widgets.FileInput(accept='.csv, .tsv, .tab, .txt')
        self.gene_annotation_file = pn.widgets.FileInput()

        self.submit_button = pn.widgets.Button(name='Submit', button_type='primary',
                                               disabled=pn.bind(enable_submit, file_input=self.input_file,
                                                                text_input=self.text_input, watch=True))
        self.submit_button.on_click(self.load_input_file)

        self.new_file_button = pn.widgets.Button(name='Process a new input file', button_type='primary')
        self.new_file_button.on_click(self.create_new_session)

        self.sample_sizes_slider = pn.widgets.DiscreteSlider(name='Subsampled regions', options=config.sampling_sizes,
                                                             bar_color='white')

        self.show_single_plots_button = pn.widgets.Button(name='Display analysis using the selected number of regions',
                                                          button_type='primary', margin=(25, 0))
        self.show_single_plots_button.on_click(self.create_single_genome_plots_by_avg)

        # Panel layouts
        self.single_multi_genome_tabs = pn.Tabs(dynamic=True, styles=config.single_multi_tabs_style)
        self.main_single_column = pn.Column(styles=config.main_column_style)
        self.single_tabs = pn.Tabs(dynamic=True, styles=config.single_tabs_style)
        self.activated_coverage_tab = 0
        self.plots_by_size_single_column = pn.Column()
        self.main_multi_column = pn.Column(styles=config.main_column_style)
        self.plots_by_ref_column = pn.Column(styles=config.main_column_style)
        self.selected_contig_column = pn.Column(styles={'padding': "10px 0 0 0"})

        # Plots cards
        self.clustermap_card = pn.Card(title='Clustermap', styles=config.plot_card_style, header_background="#2e86c1",
                                       header_color="#ffffff")
        self.jitter_card = pn.Card(title='Average synteny scores distribution', styles=config.plot_card_style,
                                   header_background="#2e86c1", header_color="#ffffff")
        self.network_card = pn.Card(title='Network', styles=config.plot_card_style, header_background="#2e86c1",
                                    header_color="#ffffff")

        download_text = 'Save file as: (if no full path, the file is saved under \'Downloads/\')'

        # Clustermap elements
        self.clustermap_plot = ""
        self.clustermap_cmap = pn.widgets.Select(value=config.clustermap_colormaps_list[0],
                                                 options=config.clustermap_colormaps_list,
                                                 name="Select colormap from the following list:")
        self.clustermap_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                         options=config.matplotlib_file_formats,
                                                         name="Select image format:")
        self.save_clustermap_file_path = pn.widgets.TextInput(name=download_text)
        self.download_clustermap_column = pn.Column()

        # Jitter plot elements
        self.jitter_plot = ""
        self.jitter_color = pn.widgets.ColorPicker(name='Select color', value='#3b89be')
        self.jitter_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                     options=config.matplotlib_file_formats, name="Select image format:")
        self.save_jitter_file_path = pn.widgets.TextInput(name=download_text)
        self.download_jiter_column = pn.Column()

        # Network plot elements
        self.network_plot = ""
        self.use_metadata = pn.widgets.Checkbox(name='Use metadata for coloring', value=False)
        self.color_edges_by_feature = pn.widgets.Checkbox(name='Color edges by feature (same/different)', value=False)
        self.metadata_colorby_card = pn.Card(title='Set the coloring by metadata', header_background="#ffffff",
                                             styles={'background': "#ffffff", 'margin': "10px", 'width': "315px"},
                                             hide_header=True, collapsed=pn.bind(change_disabled_state_inverse,
                                                                                 chkbox_state=self.use_metadata,
                                                                                 watch=True))
        self.network_node_color = pn.widgets.ColorPicker(name='Nodes color:', value='#459ED9',
                                                         disabled=pn.bind(change_disabled_state_straight,
                                                                          chkbox_state=self.use_metadata,
                                                                          watch=True))
        self.network_edge_color = pn.widgets.ColorPicker(name='Edges color:', value='#000000',
                                                         disabled=pn.bind(change_disabled_state_straight,
                                                                          chkbox_state=self.color_edges_by_feature,
                                                                          watch=True))
        self.nodes_color_by = pn.widgets.Select(options=['Select feature'], name="Color nodes by:",
                                                width=120)
        self.is_continuous = pn.widgets.Checkbox(name='Continuous feature', value=False)
        self.nodes_colormap = pn.widgets.ColorMap(name="Select colormap for nodes:",
                                                  options=pn.bind(change_continuous_state,
                                                                  chkbox_state=self.is_continuous,
                                                                  watch=True),
                                                  value=pn.bind(change_continuous_state_for_value,
                                                                chkbox_state=self.is_continuous,
                                                                watch=True))
        self.edges_color_by = pn.widgets.Select(options=['Select feature'],
                                                name="Color edges by:",
                                                width=120,
                                                disabled=pn.bind(change_disabled_state_inverse,
                                                                 chkbox_state=self.color_edges_by_feature,
                                                                 watch=True)
                                                )
        self.network_within_color = pn.widgets.ColorPicker(name='Same color:', value='#000000',
                                                           disabled=pn.bind(change_disabled_state_inverse,
                                                                            chkbox_state=self.color_edges_by_feature,
                                                                            watch=True))
        self.network_between_color = pn.widgets.ColorPicker(name='Different color:', value='#000000',
                                                            disabled=pn.bind(change_disabled_state_inverse,
                                                                             chkbox_state=self.color_edges_by_feature,
                                                                             watch=True)
                                                            )
        self.show_labels_chkbox = pn.widgets.Checkbox(name='Show sample names', value=False)
        self.network_iterations = pn.widgets.DiscreteSlider(name='Number of iterations',
                                                            options=config.network_iterations_options,
                                                            bar_color='white')
        self.network_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                      options=config.bokeh_file_formats,
                                                      name="Select image format:")
        self.save_network_file_path = pn.widgets.TextInput(name=download_text)
        self.download_network_column = pn.Column()
        self.network_pane = ""
        self.pos_dict = dict()
        self.nodes_list = []
        self.network = ""

        # Coverage plots elements
        self.coverage_plot = ""
        self.start_pos_input = pn.widgets.TextInput(name='Start position')
        self.end_pos_input = pn.widgets.TextInput(name='End position')
        self.avg_plot_chkbox = pn.widgets.Checkbox(name='Show avg. scores', value=True)
        self.avg_plot_color = pn.widgets.ColorPicker(name='Color:', value='#000080',
                                                     disabled=pn.bind(change_disabled_state_inverse,
                                                                      chkbox_state=self.avg_plot_chkbox,
                                                                      watch=True))
        self.coverage_plot_chkbox = pn.widgets.Checkbox(name='Show all synteny scores', value=True)
        self.coverage_plot_color = pn.widgets.ColorPicker(name='Color:', value='#ba55d3',
                                                          disabled=pn.bind(change_disabled_state_inverse,
                                                                           chkbox_state=self.coverage_plot_chkbox,
                                                                           watch=True))
        self.hypervar_chkbox = pn.widgets.Checkbox(name='Highlight hypervariable regions', value=True)
        self.hypervar_color = pn.widgets.ColorPicker(name='Color:', value='#00ffff',
                                                     disabled=pn.bind(change_disabled_state_inverse,
                                                                      chkbox_state=self.hypervar_chkbox,
                                                                      watch=True))
        self.hypervar_alpha_slider = pn.widgets.FloatSlider(name='Alpha transparency', start=0, end=1, step=0.1,
                                                            value=0.2,
                                                            disabled=pn.bind(change_disabled_state_inverse,
                                                                             chkbox_state=self.hypervar_chkbox,
                                                                             watch=True))
        self.hypercons_chkbox = pn.widgets.Checkbox(name='Highlight hyperconserved regions', value=True)
        self.hypercons_color = pn.widgets.ColorPicker(name='Color:', value='#fa8072',
                                                      disabled=pn.bind(change_disabled_state_inverse,
                                                                       chkbox_state=self.hypercons_chkbox,
                                                                       watch=True))
        self.hypercons_alpha_slider = pn.widgets.FloatSlider(name='Alpha transparency', start=0, end=1, step=0.1,
                                                             value=0.2,
                                                             disabled=pn.bind(change_disabled_state_inverse,
                                                                              chkbox_state=self.hypercons_chkbox,
                                                                              watch=True))
        self.coverage_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                       options=config.matplotlib_file_formats,
                                                       name="Select image format:")
        self.save_coverage_file_path = pn.widgets.TextInput(name=download_text)
        self.download_coverage_column = pn.Column()

        # Build the initial layout
        mandatory_input_title = "Mandatory input"
        self.main_area.append(pn.pane.Markdown(mandatory_input_title,
                                               styles={'font-size': "20px", 'color': config.title_purple_color,
                                                       'margin-bottom': "0"}))

        file_upload_title = "Upload SynTracker's output table 'synteny_scores_per_region' for one or multiple genomes:"
        text_input_title = "Or, if the file size is bigger than 300 Mb, enter it's full path here and press the " \
                           "Enter key:"
        self.main_area.append(pn.pane.Markdown(file_upload_title, styles={'font-size': "17px", 'margin-bottom': "0",
                                                                          'margin-top': "0"}))
        self.main_area.append(self.input_file)
        self.main_area.append(pn.pane.Markdown(text_input_title, styles={'font-size': "17px", 'margin-bottom': "0",
                                                                         'margin-top': "0"}))
        self.main_area.append(self.text_input)
        self.main_area.append(pn.Spacer(height=20))

        optional_input_title = "Optional input"
        self.main_area.append(pn.pane.Markdown(optional_input_title,
                                               styles={'font-size': "20px", 'color': config.title_purple_color,
                                                       'margin-bottom': "0"}))

        metadata_upload_title = "Upload a metadata file for the compared genomes (delimited format: csv or tsv):"
        self.main_area.append(pn.pane.Markdown(metadata_upload_title, styles={'font-size': "17px", 'margin-bottom': "0",
                                                                              'margin-top': "0"}))
        self.main_area.append(self.metadata_file)
        metadata_note = "Please note: the metadata file may contain an unlimited number of columns (metadata levels)." \
                        "\nThe first column must contain the sample IDs (identical to the sample IDs that appear in " \
                        "SynTracker's output file)."
        self.main_area.append(pn.pane.Markdown(metadata_note, styles={'font-size': "15px", 'margin-bottom': "0",
                                                                      'margin-top': "0",
                                                                              'color': config.title_red_color}))

        self.main_area.append(pn.Spacer(height=30))

        self.main_area.append(self.submit_button)

        #pn.state.onload(self.init_query_params)

        self.template.servable()

    def init_parameters(self):
        del self.score_per_region_all_genomes_df
        del self.score_per_region_selected_genome_df
        del self.avg_score_per_region_selected_genome_df

        gc.collect()

    def create_new_session(self, event):
        if self.input_file_loaded:
            pn.state.location.unsync(self.input_file)
        if self.genomes_select is not None and self.genomes_select != "":
            pn.state.location.unsync(self.genomes_select)
        self.init_parameters()
        pn.state.location.reload = True

    def init_query_params(self):
        query = {'': ''}
        #pn.state.location.update_query(input_file=None)
        pn.state.location.update_query(query)

    def submit_new_file_button(self):
        button_column = pn.Column(pn.Spacer(height=30), self.new_file_button)
        return button_column

    def get_number_of_genomes(self):

        # If the input file contains no 'Ref_genome column', add it with artificial genome name
        if 'Ref_genome' not in self.score_per_region_all_genomes_df.columns:
            self.score_per_region_all_genomes_df['Ref_genome'] = 'Reference Genome'

        self.number_of_genomes = self.score_per_region_all_genomes_df.groupby(['Ref_genome']).ngroups
        self.ref_genomes_list = list(self.score_per_region_all_genomes_df.groupby(['Ref_genome']).groups)

    def select_ref_genome(self, ref_genome):
        self.ref_genome = ref_genome
        print("Selected ref genome = " + self.ref_genome)

    def select_sampling_size(self, size):
        self.sampling_size = size
        print("Selected subsampling size = " + self.sampling_size)

    def load_input_file(self, event):
        print("\nIn load_input_file")

        # File path was given - read directly
        if self.text_input.value != "":
            self.is_file_path = 1

            # Verify that the input is an existing file path
            if os.path.exists(self.text_input.value):
                print("File path: " + self.text_input.value)
                self.filename = os.path.basename(self.text_input.value)
                print("File name: " + self.filename)
                self.display_results_page()

                pn.state.location.update_query(input_file=self.filename)
                pn.state.location.sync(self.input_file, {'filename': 'input_file'})
                print("Location: " + str(pn.state.location.pathname))
                print("Query params: " + str(pn.state.location.query_params))
                self.input_file_loaded = 1

                #self.start_process()

            else:
                title = "The requested input file does not exist, please enter again a valid file path"
                self.display_error_page(title)

        # File was given via FileInput widget
        else:
            # Filename is None - usually because of server problems
            if self.input_file.filename is None:
                print("File name: None")
                title = "Cannot upload the requested file (probably server problems) - please try again by entering " \
                        "the file's full path"
                self.display_error_page(title)

            else:
                self.filename = self.input_file.filename
                content_length = len(self.input_file.value) / 1000
                print("Content length in Kb: " + str(content_length))

                # There is filename but no content - usually happens when the file is too big
                if content_length == 0:
                    title = "Cannot upload the requested file (probably too big) - please try again by entering " \
                            "the file's full path"
                    self.display_error_page(title)

                # File has content
                else:
                    self.display_results_page()

                    pn.state.location.update_query(input_file=self.filename)
                    pn.state.location.sync(self.input_file, {'filename': 'input_file'})
                    print("Location: " + str(pn.state.location.pathname))
                    print("Query params: " + str(pn.state.location.query_params))
                    self.input_file_loaded = 1

                    #self.start_process()

        # Check if the user provided a metadata file
        if self.metadata_file.filename is not None:
            content_length = len(self.metadata_file.value) / 1000
            print("\nUploaded metadata file")
            print("Content length of metadata file in Kb: " + str(content_length))

            # There is filename but no content - usually happens when the file is too big
            if content_length == 0:
                title = "Cannot upload the requested metadata file (probably too big)"
                self.display_error_page(title)
            else:
                self.is_metadata = 1

        if self.input_file_loaded:
            self.start_process()

    def display_error_page(self, message):
        self.main_area.clear()

        self.main_area.append(pn.pane.Markdown(message, styles={'font-size': "20px"}))
        self.main_area.append(self.submit_new_file_button())

    def display_results_page(self):
        self.main_area.clear()

        title = "Loading input file: " + self.filename
        self.main_area.append(pn.pane.Markdown(title, styles={'font-size': "20px"}))

    def start_process(self):
        self.ref_genomes_list = []
        self.single_multi_genome_tabs.clear()
        self.main_single_column.clear()
        self.single_tabs.clear()
        self.plots_by_size_single_column.clear()

        # Read the file directly from the path
        col_set = ['Ref_genome', 'Sample1', 'Sample2', 'Region', 'Synteny_score']
        if self.is_file_path == 1:
            self.score_per_region_all_genomes_df = pd.read_csv(self.text_input.value,
                                                               usecols=lambda c: c in set(col_set))

        # Read the content of the uploaded file
        else:
            self.score_per_region_all_genomes_df = pd.read_csv(io.BytesIO(self.input_file.value),
                                                               usecols=lambda c: c in set(col_set))

        # Extract the number of genomes from the input file
        self.get_number_of_genomes()
        print("\nNumber of genomes: " + str(self.number_of_genomes))
        print("\nGenomes list: " + str(self.ref_genomes_list))
        
        # If a metadata file was uploaded - read the file into a DF
        if self.is_metadata:
            metadata_df = pd.read_table(io.BytesIO(self.metadata_file.value))
            print("\nMetadata before validation:")
            print(metadata_df)

            # If some samples are missing from the metadata - fill them with np.nan values
            self.metadata_dict, self.metadata_features_list, self.sample_ids_column_name = \
                dm.complete_metadata(self.score_per_region_all_genomes_df, metadata_df)

        # Initialize the dict that saves the dataframes for all the combinations of ref_genomes and sizes
        for genome in self.ref_genomes_list:
            self.avg_scores_all_genomes_all_sizes_dict[genome] = dict()
            for size in config.sampling_sizes:
                self.avg_scores_all_genomes_all_sizes_dict[genome][size] = pd.DataFrame()

        # Initialize the dict that saves the combinations of ref_genomes and sizes that have already been calculated
        for genome in self.ref_genomes_list:
            self.calculated_avg_genome_size_dict[genome] = dict()
            for size in config.sampling_sizes:
                self.calculated_avg_genome_size_dict[genome][size] = 0

        # Input file contains only one ref-genome -> present only single genome visualization
        if self.number_of_genomes == 1:
            self.ref_genome = self.ref_genomes_list[0]

            # Create the single-genome visualization layout
            self.main_single_column = self.create_single_genome_column(self.ref_genome)

            self.main_area.clear()
            self.main_area.append(self.main_single_column)

        # Input file contains more than one ref-genome -> display two tabs, for single- and multi-genome visualizations
        else:
            # Initially, the ref-genome is the first on the list
            self.ref_genome = self.ref_genomes_list[0]

            # Display a select widget to select the reference genome to visualize
            self.genomes_select = pn.widgets.Select(name='Select a reference genome to process:',
                                                    options=self.ref_genomes_list, styles={'margin': "0"})
            self.main_single_column.append(self.genomes_select)

            # Create the single-genome visualization layout for the selected ref-genome
            single_genome_column = pn.bind(self.create_single_genome_column, self.genomes_select)
            self.main_single_column.append(single_genome_column)

            pn.state.location.sync(self.genomes_select, {'value': 'ref_genome'})

            # Create the multiple-genome visualization layout
            multi_genome_layout = self.create_multi_genome_layout()

            self.single_multi_genome_tabs.append(('Single genome visualization', self.main_single_column))
            self.single_multi_genome_tabs.append(('Multiple genomes visualization', multi_genome_layout))

            self.main_area.clear()

        self.main_area.append(self.single_multi_genome_tabs)

        self.main_area.append(self.submit_new_file_button())

    def create_single_genome_column(self, ref_genome):
        # Initialize the options and clear the plots by avg area
        self.plots_by_size_single_column.clear()
        self.sample_sizes_slider.value = config.sampling_sizes[0]
        self.contigs_dict = {}
        self.contigs_list_by_length = []
        self.contigs_list_by_name = []


        # Get the selected reference genome
        self.select_ref_genome(ref_genome)

        ref_genome_title = "Reference Genome: " + self.ref_genome

        # Get the score-per-region table for the selected genome only
        self.score_per_region_selected_genome_df = ds.return_selected_genome_table(self.score_per_region_all_genomes_df,
                                                                                   self.ref_genome)

        # Get the dictionary of score-per-region tables, for each contig and the sorted contigs list
        self.contigs_dict, self.contigs_list_by_name, self.contigs_list_by_length = \
            ds.create_score_per_region_sorted_contigs_table(self.score_per_region_selected_genome_df)

        # Calculate the average score and std for the whole genome (all contigs)
        self.avg_score_genome = self.score_per_region_selected_genome_df['Synteny_score'].mean()
        self.std_score_genome = self.score_per_region_selected_genome_df['Synteny_score'].std()
        print("\nAverage score for the genome = " + str(self.avg_score_genome))
        print("Std of score for the genome = " + str(self.std_score_genome))

        # Create the df for plot presenting the number of pairs vs. subsampled regions
        pairs_num_per_sampling_size_selected_genome_df = \
            ds.create_pairs_num_per_sampling_size(self.score_per_region_selected_genome_df)

        # If the number of pairs with 40 sampled regions is smaller than 100, present also the 'All regions' bar
        # If not, do not present this bar (and remove this option from the slider)
        pairs_at_40 = pairs_num_per_sampling_size_selected_genome_df['Number of pairs'].iloc[1]
        if pairs_at_40 >= config.min_pairs_for_all_regions:
            self.sample_sizes_slider.options = config.sampling_sizes_wo_all
            self.sample_sizes_slider.value = config.sampling_sizes_wo_all[0]
            is_all_regions = 0
        else:
            is_all_regions = 1

        # Create the number of pairs vs. subsampled regions bar plot
        pairs_vs_sampling_size_bar_plot = pn.bind(ps.plot_pairs_vs_sampling_size_bar,
                                                  df=pairs_num_per_sampling_size_selected_genome_df,
                                                  sampling_size=self.sample_sizes_slider, is_all_regions=is_all_regions)

        # Create a markdown for the pairs lost percent (binded to the slider)
        binded_text = pn.bind(widgets.create_pairs_lost_text, pairs_num_per_sampling_size_selected_genome_df,
                              self.sample_sizes_slider)

        pairs_plot_column = pn.Column(pn.pane.Markdown(refs=binded_text), pairs_vs_sampling_size_bar_plot)

        initial_plots_column = pn.Column(
            pairs_plot_column,
            pn.Spacer(height=20),
            self.sample_sizes_slider,
            self.show_single_plots_button,
            self.plots_by_size_single_column,
            styles={'padding': "20px"}
            )

        self.single_tabs.clear()
        self.single_tabs.append(('Sample-pairs comparisons', initial_plots_column))

        plots_by_ref_column = self.create_coverage_plots_tab()
        self.single_tabs.append(('Coverage plots', plots_by_ref_column))

        # Create the layout of the initial plots and widgets
        single_genome_column = pn.Column(
            pn.pane.Markdown(ref_genome_title, styles={'font-size': "20px", 'color': config.title_purple_color,
                                                       'margin': "0"}),
            self.single_tabs
        )

        return single_genome_column

    def changed_active_single_tab(self, event):
        # Create the coverage plots column only the first time that this tab becomes active
        if self.single_tabs.active == 1 and self.activated_coverage_tab == 0:
            print("\nchanged_active_single_tab: calling create_coverage_plots_tab")
            self.activated_coverage_tab = 1
            self.create_coverage_plots_tab()

    def create_single_genome_plots_by_avg(self, event):

        self.sampling_size = self.sample_sizes_slider.value
        print("Selected subsampling size = " + self.sampling_size)

        self.plots_by_size_single_column.clear()
        self.jitter_card.clear()
        self.clustermap_card.clear()
        self.network_card.clear()
        self.metadata_colorby_card.clear()
        self.network_iterations.value = config.network_iterations_options[0]

        # Check if the requested genome and size have already been calculated. If so, fetch the specific dataframe
        if self.calculated_avg_genome_size_dict[self.ref_genome][self.sampling_size]:
            print("\nThe selected size (" + self.sampling_size + ") has already been calculated - retrieve it.")
            selected_genome_and_size_avg_df = \
                self.avg_scores_all_genomes_all_sizes_dict[self.ref_genome][self.sampling_size]

        else:
            # Calculate and return the dataframe with average scores for the selected genome and sampling size
            print("\nThe selected size (" + self.sampling_size + ") has not been calculated yet - calculate it.")
            selected_genome_and_size_avg_df = \
                ds.calculate_avg_scores_selected_genome_size(self.score_per_region_selected_genome_df, self.ref_genome,
                                                             self.sampling_size)
            # Save the dataframe in the main dictionary
            self.avg_scores_all_genomes_all_sizes_dict[self.ref_genome][self.sampling_size] = \
                selected_genome_and_size_avg_df
            self.calculated_avg_genome_size_dict[self.ref_genome][self.sampling_size] = 1

        # No data at the selected sampling size
        if selected_genome_and_size_avg_df.empty:
            text = "The data obtained using " + self.sampling_size + " subsampled regions is not sufficient for " \
                                                                     "further processing"
            self.plots_by_size_single_column.append(pn.pane.Markdown(text, styles={'font-size': "18px",
                                                                                   'color': config.title_red_color,
                                                                                   'margin': "0"}))

        # Enough data -> creating plots
        else:
            if self.sampling_size == 'All_regions':
                size_title = "Presenting plots using average synteny scores from all available regions:"
            else:
                size_title = "Presenting plots using average synteny scores from " + self.sampling_size + \
                             " subsampled regions:"

            self.plots_by_size_single_column.append(pn.pane.Markdown(size_title, styles={'font-size': "18px",
                                                                                         'margin': "0 0 10px 0",
                                                                                         'padding': "0"}))

            # Add the plots to the layout
            self.create_jitter_pane(selected_genome_and_size_avg_df)
            self.create_clustermap_pane(selected_genome_and_size_avg_df)
            self.create_network_pane(selected_genome_and_size_avg_df)

            plots_column = pn.Column(self.jitter_card, pn.Spacer(height=30), self.clustermap_card, pn.Spacer(height=20),
                                     self.network_card)
            self.plots_by_size_single_column.append(plots_column)

    def create_jitter_pane(self, selected_genome_and_size_avg_df):
        styling_title = "Plot styling options:"
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.jitter_color)

        save_file_title = "Download image options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_jitter)

        jitter_file = "Jitterplot_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_jitter_file_path.placeholder = jitter_file

        self.download_jiter_column = pn.Column(pn.pane.Markdown(save_file_title, styles={'font-size': "15px",
                                                                                         'font-weight': "bold",
                                                                                         'color': config.title_blue_color,
                                                                                         'margin': "0"}),
                                               self.jitter_image_format, self.save_jitter_file_path,
                                               download_button, pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_jiter_column)

        self.jitter_plot = pn.bind(ps.create_jitter_plot, avg_df=selected_genome_and_size_avg_df,
                                   color=self.jitter_color)

        #jitter_pane = pn.pane.Bokeh(self.jitter_plot)
        jitter_pane = pn.pane.Matplotlib(self.jitter_plot, height=500, dpi=300, tight=True, format='png')

        jitter_row = pn.Row(controls_col, pn.Spacer(width=250), jitter_pane, styles={'padding': "15px"})
        self.jitter_card.append(jitter_row)

    def download_jitter(self, event):
        fformat = self.jitter_image_format.value

        # Set the directory for saving
        if self.save_jitter_file_path.value == "":
            jitter_file_path = self.downloads_dir_path + self.save_jitter_file_path.placeholder + "." + fformat
        else:
            jitter_file_path = self.save_jitter_file_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, jitter_file_path, re.IGNORECASE):
                jitter_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(jitter_file_path):
                jitter_file_path = self.downloads_dir_path + jitter_file_path

        self.jitter_plot().savefig(jitter_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image was downloaded successfully under:\n" + jitter_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=20)
        self.download_jiter_column.pop(4)
        self.download_jiter_column.append(download_floatpanel)

    def create_clustermap_pane(self, selected_genome_and_size_avg_df):

        styling_title = "Clustermap styling options:"
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.clustermap_cmap)

        save_file_title = "Download image options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_clustermap)

        clustermap_file = "Clustermap_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_clustermap_file_path.placeholder = clustermap_file

        self.download_clustermap_column = pn.Column(pn.pane.Markdown(save_file_title, styles={'font-size': "15px",
                                                                                              'font-weight': "bold",
                                                                                              'color': config.title_blue_color,
                                                                                              'margin': "0"}),
                                                    self.clustermap_image_format, self.save_clustermap_file_path,
                                                    download_button, pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_clustermap_column)

        # Transform the data into a scoring matrix
        pivoted_df = selected_genome_and_size_avg_df.pivot(columns='Sample1', index='Sample2', values='Avg_synteny_score')
        scores_matrix = pivoted_df.combine_first(pivoted_df.T)

        # Check the number of columns in the matrix
        col_num = len(scores_matrix.columns)
        print("\ncreate_clustermap_pane: number of columns = " + str(col_num))

        # If the num of columns exceeds the defined maximum, do not create the clustermap and display an error message
        if col_num > config.max_clustermap_cols:
            error = "The number of compared pairs is too high to be well presented in a clustermap plot."
            clustermap_row = pn.Row(pn.pane.Markdown(error, styles={'font-size': "16px",
                                                                    'color': config.title_red_color}))

        else:
            self.clustermap_plot = pn.bind(ps.create_clustermap, matrix=scores_matrix,
                                           cmap=self.clustermap_cmap)

            clustermap_pane = pn.pane.Matplotlib(self.clustermap_plot, height=600, dpi=300, tight=True, format='png')
            clustermap_row = pn.Row(controls_col, pn.Spacer(width=110), clustermap_pane, styles={'padding': "15px"})

        self.clustermap_card.append(clustermap_row)

    def download_clustermap(self, event):
        fformat = self.clustermap_image_format.value

        # Set the directory for saving
        if self.save_clustermap_file_path.value == "":
            clustermap_file_path = self.downloads_dir_path + self.save_clustermap_file_path.placeholder + "." + fformat
        else:
            clustermap_file_path = self.save_clustermap_file_path.value

            # Add a .png suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, clustermap_file_path, re.IGNORECASE):
                clustermap_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(clustermap_file_path):
                clustermap_file_path = self.downloads_dir_path + clustermap_file_path

        self.clustermap_plot().savefig(clustermap_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image was downloaded successfully under:\n" + clustermap_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=20)
        self.download_clustermap_column.pop(4)
        self.download_clustermap_column.append(download_floatpanel)

    def create_network_pane(self, selected_genome_and_size_avg_df):
        init_button = pn.widgets.Button(name='Initialize nodes positions', button_type='primary',
                                        button_style='outline')

        styling_title = "Network customization options:"
        no_metadata_colors_row = pn.Row(self.network_node_color, pn.Spacer(width=10), self.network_edge_color)
        continuous_col = pn.Column(pn.Spacer(height=20), self.is_continuous)
        nodes_color_by_row = pn.Row(self.nodes_color_by, continuous_col)
        edges_color_by_row = pn.Row(self.edges_color_by, self.network_within_color, pn.Spacer(width=3),
                                    self.network_between_color)
        metadata_coloring_col = pn.Column(nodes_color_by_row,
                                          self.nodes_colormap,
                                          pn.Spacer(height=10),
                                          self.color_edges_by_feature,
                                          edges_color_by_row,
                                          styles={'padding': "10x"})
        self.metadata_colorby_card.append(metadata_coloring_col)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                no_metadata_colors_row,
                                pn.Spacer(height=5),
                                self.use_metadata,
                                self.metadata_colorby_card,
                                self.show_labels_chkbox,
                                self.network_iterations,
                                init_button)

        save_file_title = "Download image options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_network)

        network_file = "Network_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                            self.network_iterations.value + "_iterations"
        self.save_network_file_path.placeholder = network_file

        self.download_network_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                  styles={'font-size': "15px",
                                                                          'font-weight': "bold",
                                                                          'color': config.title_blue_color,
                                                                          'margin': "0"}),
                                                 self.network_image_format, self.save_network_file_path,
                                                 download_button, pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_network_column)

        ########################################################
        # Create a table for the network
        df_for_network = selected_genome_and_size_avg_df.loc[:, ['Sample1', 'Sample2', 'Avg_synteny_score']].copy()

        df_for_network.loc[(df_for_network['Avg_synteny_score'] < 0), 'Avg_synteny_score'] = 0
        df_for_network.loc[(df_for_network['Avg_synteny_score'] == 1), 'Avg_synteny_score'] = 0.999999

        # Set a score threshold for the connections (below it zero the weight).
        # Currently the threshold is the mean APSS
        APSS_connections_threshold = df_for_network.loc[:, 'Avg_synteny_score'].mean()
        # df_for_network['weight'] = np.negative(np.log(1 - df_for_network['Avg_synteny_score']))
        df_for_network['weight'] = np.where(df_for_network['Avg_synteny_score'] >= APSS_connections_threshold,
                                            np.negative(np.log(1 - df_for_network['Avg_synteny_score'])), 0)
        print(df_for_network)
        print("\nMean avg score: " + str(APSS_connections_threshold) + "\n")

        # Create a network using networkx
        network = nx.from_pandas_edgelist(df_for_network, source='Sample1', target='Sample2', edge_attr='weight')
        self.nodes_list = list(network.nodes)
        self.generate_rand_positions()

        # There is metadata
        if self.is_metadata:
            # Update the color nodes by- drop-down menu with the available metadata features
            self.nodes_color_by.options = self.metadata_features_list
            self.nodes_color_by.value = self.metadata_features_list[0]
            self.edges_color_by.options = self.metadata_features_list
            self.edges_color_by.value = self.metadata_features_list[0]

            # Insert the features information as nodes attributes
            for node in self.nodes_list:
                for feature in self.metadata_features_list:
                    network.nodes[node][feature] = self.metadata_dict[feature][node]

                # Add node attribute 'SampleID' for the hover tooltip
                network.nodes[node]['SampleID'] = node

        # No metadata
        else:
            self.use_metadata.disabled = True

        # Create the network plot using the selected parameters
        self.network_plot = pn.bind(ps.cretae_network_plot, network=network, is_metadata=self.use_metadata,
                                    nodes_feature=self.nodes_color_by, is_continuous=self.is_continuous,
                                    cmap=self.nodes_colormap,
                                    node_color=self.network_node_color, edge_color=self.network_edge_color,
                                    is_edge_colorby=self.color_edges_by_feature, edges_feature=self.edges_color_by,
                                    within_edge_color=self.network_within_color,
                                    between_edge_color=self.network_between_color,
                                    iterations=self.network_iterations, pos_dict=self.pos_dict,
                                    show_labels=self.show_labels_chkbox, metadata_dict=self.metadata_dict)
        self.network_pane = pn.pane.HoloViews(self.network_plot, height=600, width=700, sizing_mode="fixed")

        init_button.on_click(partial(self.init_positions, network))

        network_row = pn.Row(controls_col, pn.Spacer(width=10), self.network_pane, styles={'padding': "15px"})
        self.network_card.append(network_row)

    def download_network(self, event):
        fformat = self.network_image_format.value

        # Set the directory for saving
        if self.save_network_file_path.value == "":
            network_file_path = self.downloads_dir_path + self.save_network_file_path.placeholder + "." + fformat
        else:
            network_file_path = self.save_network_file_path.value

            # Add a .png suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, network_file_path, re.IGNORECASE):
                network_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(network_file_path):
                network_file_path = self.downloads_dir_path + network_file_path

        # Save the network plot in the requested format
        hv.extension('bokeh')
        hv.plotting.bokeh.ElementPlot.width = 1200
        hv.plotting.bokeh.ElementPlot.height = 1000
        fig = hv.render(self.network_plot())

        # Save the figure in png format
        if fformat == 'png':
            export_png(fig, filename=network_file_path)

        # Save the figure in svg format
        else:
            fig.output_backend = "svg"
            export_svgs(fig, filename=network_file_path)

        download_message = "The image was downloaded successfully under:\n" + network_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=20)
        self.download_network_column.pop(4)
        self.download_network_column.append(download_floatpanel)

    def init_positions(self, network, event):
        self.generate_rand_positions()
        self.update_network_plot(network)

    def generate_rand_positions(self):
        self.pos_dict = {}

        for node in self.nodes_list:
            pos_x = generate_rand_pos()
            pos_y = generate_rand_pos()
            pos_tuple = (pos_x, pos_y)
            self.pos_dict[node] = pos_tuple

    # Update the network plot using the selected parameters and the new positions dict
    def update_network_plot(self, network):
        self.network_plot = pn.bind(ps.cretae_network_plot, network=network, is_metadata=self.use_metadata,
                                    nodes_feature=self.nodes_color_by, is_continuous=self.is_continuous,
                                    cmap=self.nodes_colormap,
                                    node_color=self.network_node_color, edge_color=self.network_edge_color,
                                    is_edge_colorby=self.color_edges_by_feature, edges_feature=self.edges_color_by,
                                    within_edge_color=self.network_within_color,
                                    between_edge_color=self.network_between_color,
                                    iterations=self.network_iterations, pos_dict=self.pos_dict,
                                    show_labels=self.show_labels_chkbox, metadata_dict=self.metadata_dict)
        self.network_pane.object = self.network_plot

    def create_coverage_plots_tab(self):

        plots_by_ref_column = pn.Column(styles=config.main_column_style)

        contigs_num = len(self.contigs_list_by_name)

        # If there is more than one contig - display a select widget to select the contig
        if contigs_num > 1:

            # Create a drop-down menu of the contigs
            contig_select = pn.widgets.Select(name='Select a contig:', options=self.contigs_list_by_name,
                                              value=self.contigs_list_by_name[0], styles={'margin': "0"})
            # When the selection of contig is changed, a new column is created for the selected contig
            contig_select_watcher = contig_select.param.watch(partial(self.changed_contig, contig_select),
                                                              'value', onlychanged=False)

            sorting_options = ['Contig name', 'Contig length']
            sorting_select = pn.widgets.Select(name='Sort by:', options=sorting_options,
                                               styles={'margin': "0"})
            # When the selection of sorting method is changed, sort the contigs in the contigs-select drop-down menu
            # accordingly
            sorting_select_watcher = sorting_select.param.watch(partial(self.changed_sorting, sorting_select,
                                                                        contig_select),
                                                                'value', onlychanged=False)

            contig_select_row = pn.Row(contig_select, pn.Spacer(width=20), sorting_select)

            plots_by_ref_column.append(contig_select_row)

            contig_select.param.trigger('value')

        # Create the coverage plot for the Ref-genome
        else:
            self.create_selected_contig_column(self.contigs_list_by_name[0])

        plots_by_ref_column.append(self.selected_contig_column)

        return plots_by_ref_column

    def changed_sorting(self, sorting_select, contig_select, event):
        sorting_method = sorting_select.value
        print("\nChanged contigs sorting method. Sort by: " + sorting_method)

        # Change the order of the contigs selection and present the first one in the new order
        if sorting_method == "Contig name":
            contig_select.options = self.contigs_list_by_name
            contig_select.value = self.contigs_list_by_name[0]
        else:
            contig_select.options = self.contigs_list_by_length
            contig_select.value = self.contigs_list_by_length[0]

    def changed_contig(self, contig_name, event):
        print("\nChanged_contig, contig name:")
        print(contig_name.value)
        self.create_selected_contig_column(contig_name.value)

    def create_selected_contig_column(self, contig_name):

        score_per_pos_contig = self.contigs_dict[contig_name]

        self.selected_contig_column.clear()

        contig_name_title = "Contig name: " + contig_name
        self.selected_contig_column.append(pn.pane.Markdown(contig_name_title,
                                                            styles={'font-size': "17px",
                                                                    'color': config.title_purple_color,
                                                                    'margin': "5px 5px 0px 5px",
                                                                    'padding-bottom': "0px"}))

        # Find contig length by the last position
        score_per_pos_contig['Position'] = score_per_pos_contig['Position'].astype(int)
        score_per_pos_contig = score_per_pos_contig.sort_values('Position')
        print("\nLast position: " + str(score_per_pos_contig.iloc[-1]['Position']))
        contig_length = score_per_pos_contig.iloc[-1]['Position'] + config.region_length
        contig_length_title = "Contig length: " + str(contig_length) + " bp"
        self.selected_contig_column.append(pn.pane.Markdown(contig_length_title,
                                                            styles={'font-size': "16px", 'margin': "0px 3px 3px 5px",
                                                                    'padding-top': "0px"}))

        # Create the customization column

        # Set the length range
        start_pos = '0'
        end_pos = str(contig_length)
        self.start_pos_input.placeholder = start_pos
        self.start_pos_input.value = start_pos
        self.end_pos_input.placeholder = end_pos
        self.end_pos_input.value = end_pos

        reset_range_button = pn.widgets.Button(name='Reset range', button_type='primary',
                                               styles={'margin-top': "22px"})
        reset_range_button.on_click(partial(self.reset_range, start_pos, end_pos))
        pos_range_cust_row = pn.Row(pn.pane.Markdown("Set genome length range:",
                                                     styles={'font-size': "14px", 'margin': "10px 5px 5px 5px"}),
                                    pn.Spacer(width=10), self.start_pos_input, pn.Spacer(width=5), self.end_pos_input,
                                    pn.Spacer(width=5), reset_range_button,
                                    styles={'margin-left': "5px"})

        avg_plot_chkbox_col = pn.Column(pn.Spacer(height=20), self.avg_plot_chkbox)
        avg_plot_cust_row = pn.Row(avg_plot_chkbox_col, pn.Spacer(width=10), self.avg_plot_color,
                                   styles={'margin-left': "5px"})

        coverage_plot_chkbox_col = pn.Column(pn.Spacer(height=20), self.coverage_plot_chkbox)
        coverage_plot_cust_row = pn.Row(coverage_plot_chkbox_col, pn.Spacer(width=10), self.coverage_plot_color,
                                        styles={'margin-left': "5px"})

        hypervar_chkbox_col = pn.Column(pn.Spacer(height=20), self.hypervar_chkbox)
        hypervar_cust_row = pn.Row(hypervar_chkbox_col, pn.Spacer(width=10), self.hypervar_color, pn.Spacer(width=5),
                                   self.hypervar_alpha_slider, styles={'margin-left': "5px"})

        hypercons_chkbox_col = pn.Column(pn.Spacer(height=20), self.hypercons_chkbox)
        hypercons_cust_row = pn.Row(hypercons_chkbox_col, pn.Spacer(width=10), self.hypercons_color, pn.Spacer(width=5),
                                    self.hypercons_alpha_slider, styles={'margin-left': "5px"})

        styling_title = "Customization options:"
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "10px 5px 5px 5px"}),
                                pos_range_cust_row,
                                avg_plot_cust_row,
                                coverage_plot_cust_row,
                                hypervar_cust_row,
                                hypercons_cust_row)

        save_file_title = "Download image options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(partial(self.download_coverage_plot, contig_name))

        coverage_file = "Coverage_plot_" + self.ref_genome + "_" + contig_name
        self.save_coverage_file_path.placeholder = coverage_file

        self.download_coverage_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                   styles={'font-size': "15px", 'font-weight': "bold",
                                                                           'color': config.title_blue_color,
                                                                           'margin': "10px 5px 5px 5px"}),
                                                  self.coverage_image_format, self.save_coverage_file_path,
                                                  pn.Spacer(height=10),
                                                  download_button, pn.pane.Markdown())

        # Prepare the data structures necessary for the coverage plots of a specific contig
        # Calculate the average synteny scores for each position
        avg_score_per_pos_contig = \
            score_per_pos_contig[['Contig_name', 'Position', 'Synteny_score']]. \
                sort_values(['Position']).groupby(['Position']).mean(numeric_only=True).reset_index(). \
                rename(columns={"Synteny_score": "Avg_synteny_score"})
        #print("\nAverage df:")
        #print(avg_score_per_pos_contig)

        # Fill the missing positions with score=0 (default jump=5000)
        avg_score_per_pos_contig = avg_score_per_pos_contig. \
            merge(how='right', on='Position',
                  right=pd.DataFrame(
                      {'Position': np.arange(avg_score_per_pos_contig.iloc[0]['Position'],
                                             avg_score_per_pos_contig.iloc[-1]['Position'] + 2 * config.region_length,
                                             config.region_length)
                       })).sort_values(by='Position').reset_index(). \
            drop(['index'], axis=1)
        avg_score_per_pos_contig['Position'] = avg_score_per_pos_contig['Position'].astype(int)

        #print("\nAfter filling missing positions:")
        #print(avg_score_per_pos_contig)

        # Create the data for highlighting hypervariable and hyperconserved regions
        # For the calculation, consider the whole data (full range)
        height = 1
        bottom = 0
        min_score = score_per_pos_contig['Synteny_score'].min()
        if min_score < 0:
            height += abs(min_score) + 0.05
            bottom = min_score - 0.05
        print("Min score = " + str(min_score))
        print("Height = " + str(height))
        avg_score_per_pos_contig['Hypervariable'] = np.where(
            avg_score_per_pos_contig['Avg_synteny_score'] < (self.avg_score_genome - 0.75 * self.std_score_genome),
            height, 0)
        avg_score_per_pos_contig['Hyperconserved'] = np.where(
            avg_score_per_pos_contig['Avg_synteny_score'] > (self.avg_score_genome + 0.75 * self.std_score_genome),
            height, 0)
        #print(avg_score_per_pos_contig)

        # Create the coverage plot for the selected contig of the Ref-genome with the customized parameters
        self.coverage_plot = pn.bind(
            ps.create_coverage_plot, contig_name=contig_name, score_per_pos_contig=score_per_pos_contig,
            avg_score_per_pos_contig=avg_score_per_pos_contig,
            start_pos=self.start_pos_input, end_pos=self.end_pos_input,
            show_avg=self.avg_plot_chkbox, avg_color=self.avg_plot_color,
            show_scores=self.coverage_plot_chkbox, scores_color=self.coverage_plot_color,
            show_hyper_var=self.hypervar_chkbox, hyper_var_color=self.hypervar_color,
            hyper_var_alpha=self.hypervar_alpha_slider,
            show_hyper_cons=self.hypercons_chkbox, hyper_cons_color=self.hypercons_color,
            hyper_cons_alpha=self.hypercons_alpha_slider, bottom_val=bottom)

        coverage_pane = pn.pane.Matplotlib(self.coverage_plot, height=600, dpi=300, tight=True, format='png')

        coverage_plot_row = pn.Row(styles={'padding': "5px 0 0 0"})
        coverage_plot_row.append(coverage_pane)

        self.selected_contig_column.append(coverage_plot_row)
        self.selected_contig_column.append(pn.Spacer(width=20))
        self.selected_contig_column.append(styling_col)
        self.selected_contig_column.append(self.download_coverage_column)

    def reset_range(self, start, end, event):
        self.start_pos_input.value = start
        self.end_pos_input.value = end

    def download_coverage_plot(self, contig_name, event):
        fformat = self.coverage_image_format.value

        # Set the directory for saving
        if self.save_coverage_file_path.value == "":
            coverage_file_path = self.downloads_dir_path + self.save_coverage_file_path.placeholder + "." + fformat
        else:
            coverage_file_path = self.save_coverage_file_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, coverage_file_path, re.IGNORECASE):
                coverage_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(coverage_file_path):
                coverage_file_path = self.downloads_dir_path + coverage_file_path

        self.coverage_plot().savefig(coverage_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image was downloaded successfully under:\n" + coverage_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_coverage_column.pop(5)
        self.download_coverage_column.append(download_floatpanel)

    def create_multi_genome_layout(self):

        multi_genome_layout = pn.Column()

        return multi_genome_layout



