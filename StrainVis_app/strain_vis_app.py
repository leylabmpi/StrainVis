import gc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from dna_features_viewer import GraphicFeature, GraphicRecord
import panel as pn
import pandas as pd
import numpy as np
import io
import os
import re
import time
import threading
from functools import partial
import networkx as nx
import random
from scipy.stats import mannwhitneyu, ranksums, spearmanr
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import StrainVis_app.config as config
import StrainVis_app.data_manipulation_single as ds
import StrainVis_app.data_manipulation_multi as dm
import StrainVis_app.plots_single_genome as ps
import StrainVis_app.plots_multi_genomes as pm
import StrainVis_app.widgets as widgets
import matplotlib
import textwrap
matplotlib.use('agg')
pn.extension(disconnect_notification='Connection lost, try reloading the page!')
pn.extension('floatpanel')


def destroyed(session_context):
    print("\n\n\nThe session is closed...")


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


def change_collapse_state(selection_val):
    if selection_val == 'All species':
        return True
    else:
        return False


def change_disabled_state_threshold(value):
    if value == 'Define another threshold':
        return False
    else:
        return True


def change_disabled_state_custom_colormap(value):
    if value == config.categorical_colormap_dict['Define custom colormap']:
        return False
    else:
        return True


def enable_submit(syn_file_input, syn_text_input, ani_file_input, ani_text_input):
    if syn_file_input is not None or syn_text_input != "" or ani_file_input is not None or ani_text_input != "":
        return False
    else:
        return True


def generate_rand_pos():
    rand = random.random() * 2 - 1
    return rand


def return_p_value_ranksums(data1, data2):
    stat, p_val = ranksums(data1, data2)
    return stat, p_val


def return_p_value_mannwhitneyu(data1, data2):
    u, p_val = mannwhitneyu(data1, data2)
    return u, p_val


def return_significance(row):
    if str(row['P_value']) != "NaN":
        # e-10
        if row['P_value'] < 0.0000000005:
            return "***"
        # e-5
        elif row['P_value'] < 0.00005:
            return "**"
        # e-2
        elif row['P_value'] < 0.05:
            return "*"
        else:
            return "NS"
    else:
        return ""


def is_missing_value(val):
    str_val = str(val)
    if str_val == 'nan' or str_val == 'NA' or str_val == 'NAN' or re.search(r'^unknown', str_val, flags=re.IGNORECASE):
        return True
    else:
        return False


def category_by_feature(row, feature, metadata_dict):
    if is_missing_value(metadata_dict[feature][row['Sample1']]) \
            or is_missing_value(metadata_dict[feature][row['Sample2']]):
        return 'Unknown'
    else:
        if metadata_dict[feature][row['Sample1']] == metadata_dict[feature][row['Sample2']]:
            return 'Same ' + feature
        else:
            return 'Different ' + feature


def wrap_label(s, width=12):
    if not s:
        return ""
    return "\n".join(textwrap.wrap(s, width=width, break_long_words=True))


def set_features_in_range(features, x_min, x_max):
    filtered_features = []

    for feature in features:
        # Append only features within the selected sequence range
        if feature.end >= x_min and feature.start <= x_max:
            filtered_features.append(feature)

    return filtered_features


def load_genes_from_gff(gff_text):
    #print("\nload_genes_from_gff:")
    contigs_list = []
    contigs_dict = dict()
    for line in gff_text.splitlines():
        if line.startswith("#"):
            # Found new contig
            if re.search("^##sequence-region", line):
                m = re.search(r'^##sequence-region\s+(\S+)\s+\d+\s+\d+', line)
                if m:
                    contig = m.group(1)
                    contigs_list.append(contig)
                    contigs_dict[contig] = []  # Initialize the new contig's features array
                    #print("Found new contig: " + contig)

            # If found a FASTA section - stop reading the file
            elif re.search("^##FASTA", line):
                break

            else:
                continue

        # Real informative line
        else:
            seqid, source, feature_type, start, end, score, strand, phase, attrs = line.strip().split("\t")

            # Consider only CDC records
            #if feature_type != "gene":
            if feature_type != "CDS":
                continue

            # Verify that the current record belongs to a contig from the list
            if seqid in contigs_list:
                is_gene = 0
                contig = seqid

                strand = 1 if strand == "+" else -1

                for item in attrs.split(";"):
                    if item.startswith("ID="):
                        cds_id = item.replace("ID=", "")
                    elif item.startswith("gene="):
                        gene = item.replace("gene=", "")
                        is_gene = 1

                # If there is a gene name, use it as label. If not, use the CDS ID
                if is_gene:
                    label = gene
                else:
                    label = cds_id
                wrapped_label = wrap_label(label, width=14)  # Wrap labels to make them narrow and stacked

                contigs_dict[contig].append(
                    GraphicFeature(
                        start=int(start),
                        end=int(end),
                        strand=strand,
                        color="#66c2a5",
                        label=wrapped_label,
                        fontdict={'fontsize': 7}
                    )
                )

    return contigs_dict


def ensure_bottom_margin_for_texts(fig, ax_list, pad_pixels=6, min_bottom=0.05, max_bottom=0.30):
    """
    Increase figure bottom margin so that texts living below the lowest axis in ax_list
    are not clipped.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    ax_list : sequence of matplotlib.axes.Axes
        Axes to inspect for text artists (we check .texts, xticklabels, xlabel)
    pad_pixels : int
        Extra padding in pixels between lowest text bbox and the axes area.
    min_bottom : float
        Minimum fraction of figure height to use as bottom margin (safety).
    max_bottom : float
        Maximum allowed bottom fraction.
    """
    # Force a draw so text positions are updated
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_h = fig.bbox.height  # pixels

    # Collect bottom-most pixel coordinate among relevant texts
    y_bottoms = []

    for ax in ax_list:
        # axis texts (DFV labels)
        for txt in ax.texts:
            if not txt.get_visible():
                continue
            try:
                bbox = txt.get_window_extent(renderer)
            except Exception:
                continue
            # bbox.y0 is the lower pixel coordinate of the text (origin at bottom-left)
            y_bottoms.append(bbox.y0)

        # xtick labels on this axis
        for lbl in ax.get_xticklabels():
            if not lbl.get_visible():
                continue
            bbox = lbl.get_window_extent(renderer)
            y_bottoms.append(bbox.y0)

        # x-axis label (xlabel)
        lbl = ax.xaxis.get_label()
        if lbl.get_visible():
            bbox = lbl.get_window_extent(renderer)
            y_bottoms.append(bbox.y0)

    if not y_bottoms:
        return  # nothing to do

    lowest_px = min(y_bottoms)  # smallest y (closest to bottom)
    # convert to fraction of figure height
    lowest_frac = lowest_px / fig_h

    # We need the subplot bottom to be *above* the lowest_frac by pad_pixels
    pad_frac = pad_pixels / fig_h
    required_bottom = lowest_frac + pad_frac

    # clamp to reasonable range
    required_bottom = max(required_bottom, min_bottom)
    required_bottom = min(required_bottom, max_bottom)

    # current bottom
    left, bottom, right, top = fig.subplotpars.left, fig.subplotpars.bottom, fig.subplotpars.right, fig.subplotpars.top

    # If required_bottom is larger than current bottom, increase it
    if required_bottom > bottom:
        plt.subplots_adjust(bottom=required_bottom)


class StrainVisApp:

    def __init__(self):

        # Variables
        self.is_syntracker_file_path = 0
        self.is_ani_file_path = 0
        self.input_mode = "SynTracker"  # Options: SynTracker / ANI / both
        self.active_input_mode = "SynTracker"  # In case of both, can be: SynTracker / ANI
        self.score_type = "APSS"  # Options: APSS / ANI
        self.syntracker_filename = ""
        self.ani_filename = ""
        self.input_file_loaded = 0
        self.syntracker_loaded = 0
        self.ani_loaded = 0
        self.is_metadata = 0
        self.valid_metadata = 1
        self.metadata_dict = dict()
        self.groups_per_feature_dict = dict()
        self.metadata_features_list = []
        self.number_of_genomes = 0
        self.number_of_genomes_syntracker = 0
        self.number_of_genomes_ani = 0
        self.ref_genomes_list = []
        self.ref_genomes_list_by_pairs_num = []
        self.ref_genomes_list_syntracker = []
        self.ref_genomes_list_by_pairs_num_syntracker = []
        self.ref_genomes_list_ani = []
        self.ref_genomes_list_by_pairs_num_ani = []
        self.annotation_per_ref_genome_dict = dict()
        self.selected_genomes_subset = []
        self.sorted_selected_genomes_subset = []
        self.sorted_selected_genomes_subset_ani = []
        self.selected_subset_species_num = 0
        self.species_num_at_sampling_size = 0
        self.species_num_in_subset_syntracker = 0
        self.species_num_in_subset_ani = 0
        self.ref_genome = ""
        self.sampling_size = ""
        self.sampling_size_multi = ""
        self.score_per_region_all_genomes_df = pd.DataFrame()
        self.score_per_region_genomes_subset_df = pd.DataFrame()
        self.score_per_region_selected_genome_df = pd.DataFrame()
        self.ani_scores_all_genomes_df = pd.DataFrame()
        self.ani_scores_genomes_subset_df = pd.DataFrame()
        self.ani_scores_selected_genome_df = pd.DataFrame()
        self.score_per_pos_contig = pd.DataFrame()
        self.score_per_pos_contig_filtered = pd.DataFrame()
        self.avg_score_per_pos_contig = pd.DataFrame()
        self.avg_score_per_pos_contig_filtered = pd.DataFrame()
        self.genomes_subset_selected_size_APSS_df = pd.DataFrame()
        self.pairs_num_per_sampling_size_multi_genomes_df = pd.DataFrame()
        self.boxplot_p_values_df = pd.DataFrame()
        self.boxplot_p_values_df_ani = pd.DataFrame()
        self.APSS_by_genome_all_sizes_dict = dict()
        self.APSS_all_genomes_all_sizes_dict = dict()
        self.calculated_APSS_genome_size_dict = dict()
        self.calculated_APSS_all_genomes_size_dict = dict()
        self.working_directory = os.getcwd()
        self.downloads_dir_path = self.working_directory + config.downloads_dir
        self.contigs_list_by_name = []
        self.contigs_list_by_length = []
        self.total_pairs_genome = 0
        self.avg_score_genome = 0
        self.std_score_genome = 0
        self.top_percentile = 0
        self.bottom_percentile = 0
        self.median_counts = 0
        self.median_counts_filtered = 0
        self.bottom_percentile_counts = 0
        self.bottom_percentile_counts_filtered = 0
        self.genomes_select_watcher = ""
        self.genomes_sort_select_watcher = ""
        self.genomes_sort_select_multi_watcher = ""
        self.threshold_select_watcher = ""
        self.threshold_input_watcher = ""
        self.highlight_sample_watcher = ""
        self.threshold_select_ani_watcher = ""
        self.threshold_input_ani_watcher = ""
        self.highlight_sample_ani_watcher = ""
        self.feature_select_watcher = ""
        self.feature_select_ani_watcher = ""
        self.continuous_network_watcher = ""
        self.colormap_watcher = ""
        self.custom_colormap_watcher = ""
        self.nodes_colorby_watcher = ""
        self.nodes_highlight_by_watcher = ""
        self.continuous_network_ani_watcher = ""
        self.colormap_ani_watcher = ""
        self.custom_colormap_ani_watcher = ""
        self.nodes_colorby_ani_watcher = ""
        self.nodes_highlight_by_ani_watcher = ""
        self.feature_colormap_watcher = ""
        self.custom_colormap_clustermap_watcher = ""
        self.feature_colormap_ani_watcher = ""
        self.custom_colormap_clustermap_ani_watcher = ""
        self.continuous_clustermap_watcher = ""
        self.continuous_clustermap_ani_watcher = ""
        self.color_by_feature_clustermap_watcher = ""
        self.color_by_feature_clustermap_ani_watcher = ""
        self.jitter_feature_watcher = ""
        self.jitter_feature_ani_watcher = ""
        self.synteny_per_pos_feature_select_watcher = ""
        self.visited_multi_genome_tab = 0
        self.activated_synteny_per_pos_tab = 0
        self.clicked_button_display_APSS = 0
        self.visited_ANI_tab = 0
        self.clicked_button_display_APSS_multi = 0
        self.visited_ANI_tab_multi = 0

        # Bootstrap template
        self.template = pn.template.VanillaTemplate(
            title='StrainVis',
            site_url="strain_vis",
        )

        # Apply custom CSS to adjust button font size in the header
        button_css = '''
            /* Target buttons in the header of the Vanilla template */
            .bk-btn { 
                font-size: 20px !important;
                color: white; 
                margin-right: 20px
            }
            .bk-btn:hover { 
                color: #87cefa !important; 
            }
            .bk-btn.bk-active {
                color: #87cefa !important; 
            }
        '''

        self.header_buttons = pn.widgets.ToggleGroup(name='header_buttons', options=['Home', 'Help'], behavior="radio",
                                                     button_type='primary', stylesheets=[button_css])
        self.header_buttons.param.watch(self.load_correct_page, 'value')
        self.menu_row = pn.Row(self.header_buttons, styles=config.menu_row_style)
        self.header_container = pn.Column(
            self.menu_row,
            styles=config.header_container_style)
        self.template.header.append(self.header_container)

        self.main_container = pn.Column(sizing_mode='stretch_width')
        self.main_area = pn.Column(styles=config.main_area_style)
        
        # Reading the manual.md into a variable and display it in a markdown pane
        manual_file_path = self.working_directory + config.manual_file
        with open(manual_file_path, 'r') as manual:
            manual_content = manual.read()
        self.help_area = pn.Column(styles=config.main_area_style)
        #self.help_area.append(pn.pane.Markdown(manual_content, styles={'font-size': "16px"}))
        self.help_area.append(pn.pane.HTML(manual_content, styles={'font-size': "16px"}))

        self.main_container.append(self.main_area)
        self.template.main.append(self.main_container)

        # Homepage layouts
        self.mandatory_input_card = pn.Card(hide_header=True, styles={'margin': "10px", 'padding': "10px",
                                                                      'width': "1000px"})
        self.optional_input_card = pn.Card(hide_header=True, styles={'margin': "10px", 'padding': "10px",
                                                                     'width': "1000px"})
        self.ani_upload_row = pn.Row(styles={'font-size': "16px", 'margin-bottom': "0", 'margin-top': "0"})

        # Widgets
        radio_group_css = '''
            .bk-input-group.bk-inline > * {
                margin-right: 20px;  /* Adjust spacing between buttons */
                font-size: 16px;
            }
        '''
        self.input_type_radio_group = pn.widgets.ToggleGroup(name='RadioBoxGroup',
                                                             options=['SynTracker output file',
                                                                      'ANI file',
                                                                      'Both SynTracker and ANI files (for the same species)'],
                                                             behavior="radio", widget_type="box", inline=True,
                                                             stylesheets=[radio_group_css])
        self.input_type_radio_group.param.watch(self.update_input_card, 'value', onlychanged=True)
        self.SynTracker_text_input = pn.widgets.TextInput(name='', placeholder='Enter SynTracker file path here...')
        self.SynTracker_input_file = pn.widgets.FileInput(accept='.csv, .tab, .txt')
        self.ANI_text_input = pn.widgets.TextInput(name='', placeholder='Enter ANI file path here...')
        self.ANI_input_file = pn.widgets.FileInput(accept='.tsv, .tab, .txt')
        self.metadata_file = pn.widgets.FileInput(accept='.csv, .tsv, .tab, .txt')

        self.submit_button = pn.widgets.Button(name='Submit', button_type='primary',
                                               disabled=pn.bind(enable_submit,
                                                                syn_file_input=self.SynTracker_input_file,
                                                                syn_text_input=self.SynTracker_text_input,
                                                                ani_file_input=self.ANI_input_file,
                                                                ani_text_input=self.ANI_text_input,
                                                                watch=True))
        self.submit_button.on_click(self.load_input_file)

        self.new_file_button = pn.widgets.Button(name='Process a new input file', button_type='primary')
        self.new_file_button.on_click(self.create_new_session)

        self.genomes_select = pn.widgets.Select(name='Select a species to process:', value=None,
                                                options=[], styles={'margin': "0"})
        self.genomes_sort_select = pn.widgets.Select(name='Sort by:', value=config.genomes_sorting_options[0],
                                                     options=config.genomes_sorting_options, styles={'margin': "0"})
        self.sample_sizes_slider = pn.widgets.DiscreteSlider(name='Subsampled regions', options=config.sampling_sizes,
                                                             bar_color='white')
        self.show_single_plots_button = pn.widgets.Button(name='Display plots using the selected number of regions',
                                                          button_type='primary', margin=(25, 0))
        self.show_single_plots_button.on_click(self.create_single_genome_plots_by_APSS)

        self.all_or_subset_radio = pn.widgets.RadioBoxGroup(name='all_or_subset',
                                                            options=['All species', 'Select a subset of species'],
                                                            inline=False, styles={'font-size': "16px"})
        self.genomes_select_card = pn.Card(title='Species subset selection',
                                           collapsed=pn.bind(change_collapse_state,
                                                             selection_val=self.all_or_subset_radio, watch=True),
                                           styles={'margin': "5px 0 5px 10px", 'width': "420px"}
                                           )
        self.genomes_subset_select = pn.widgets.MultiSelect(options=[], height=300,
                                                            disabled=pn.bind(change_collapse_state,
                                                                             selection_val=self.all_or_subset_radio,
                                                                             watch=True))
        self.genomes_sort_select_multi = pn.widgets.Select(name='Sort species by:',
                                                           value=config.genomes_sorting_options[0],
                                                           options=config.genomes_sorting_options, width=200,
                                                           styles={'margin': "0"})
        self.update_genomes_selection_button = pn.widgets.Button(name='Update species selection/sorting',
                                                                 button_type='primary',
                                                                 styles={'margin': "12px 0 12px 10px"})
        self.update_genomes_selection_button.on_click(self.update_genomes_selection)
        self.sample_sizes_slider_multi = pn.widgets.DiscreteSlider(name='Subsampled regions',
                                                                   options=config.sampling_sizes, bar_color='white',
                                                                   styles={'font-size': "16px", 'width': "450px",
                                                                           'text-align': "center"})
        self.show_box_plot_multi_button = pn.widgets.Button(name='Display plot using the selected number of regions',
                                                            button_type='primary', margin=(25, 0))
        self.show_box_plot_multi_button.on_click(self.create_multi_genomes_plots_by_APSS)

        # Panel layouts
        single_multi_tabs_css = '''
                    .bk-tab {
                        font-size: 24px;
                        color: #808080;
                        background-color: #f9f9f9;
                        padding: 5px 10px 5px 10px;
                        margin-right: 5px;
                    }
                    .bk-tab.bk-active {
                        background-color: #0072b5;
                        color: #ffffff;
                    }
        '''
        single_tabs_css = '''
                            .bk-tab {
                                font-size: 20px;
                                color: #808080;
                                background-color: #f9f9f9;
                                padding: 5px 10px 5px 10px;
                                margin-right: 5px;
                            }
                            .bk-tab.bk-active {
                                background-color: #800080;
                                color: #ffffff;
                            }
        '''
        syntracker_ani_tabs = '''
                            .bk-tab {
                                font-size: 20px;
                                color: #808080;
                                background-color: #fdfdfd;
                                padding: 5px 10px 5px 10px;
                                margin-right: 5px;
                            }
                            .bk-tab.bk-active {
                                background-color: #0072b5;
                                color: #ffffff;
                            }
        '''

        self.single_multi_genome_tabs = pn.Tabs(dynamic=True, stylesheets=[single_multi_tabs_css])
        self.synteny_ani_single_tabs = pn.Tabs(dynamic=True, stylesheets=[syntracker_ani_tabs])
        self.synteny_ani_multi_tabs = pn.Tabs(dynamic=True, stylesheets=[syntracker_ani_tabs])
        self.main_single_column = pn.Column(styles=config.main_column_style)
        self.main_multi_column = pn.Column(styles=config.main_column_style)
        self.ref_genome_column = pn.Column()
        self.synteny_single_initial_plots_column = pn.Column(styles=config.synteny_only_single_style)
        self.ani_single_plots_column = pn.Column(styles=config.ani_or_multi_style)
        self.combined_single_plots_column = pn.Column(styles=config.both_mode_other_style)
        self.synteny_single_tabs = pn.Tabs(dynamic=True, stylesheets=[single_tabs_css])
        self.APSS_analyses_single_column = pn.Column(styles=config.single_tabs_style)
        self.synteny_per_pos_plot_column = pn.Column(styles=config.single_tabs_style)
        self.plots_by_size_single_column = pn.Column()
        self.synteny_multi_initial_plots_column = pn.Column(styles=config.ani_or_multi_style)
        self.ani_multi_plots_column = pn.Column(styles=config.ani_or_multi_style)
        self.plots_by_size_multi_column = pn.Column()
        self.selected_contig_column = pn.Column(styles={'padding': "10px 0 0 0"})

        # Plots cards
        self.clustermap_card = pn.Card(title='Clustered heatmap plot', styles=config.plot_card_style,
                                       header_background="#2e86c1", header_color="#ffffff")
        self.jitter_card = pn.Card(title='APSS distribution plot', styles=config.plot_card_style,
                                   header_background="#2e86c1", header_color="#ffffff")
        self.clustermap_card_ani = pn.Card(title='Clustered heatmap plot', styles=config.plot_card_style,
                                       header_background="#2e86c1", header_color="#ffffff")
        self.jitter_card_ani = pn.Card(title='ANI distribution plot', styles=config.plot_card_style,
                                   header_background="#2e86c1", header_color="#ffffff")
        self.network_card = pn.Card(title='Network plot', styles=config.plot_card_style, header_background="#2e86c1",
                                    header_color="#ffffff")
        self.network_card_ani = pn.Card(title='Network plot', styles=config.plot_card_style, header_background="#2e86c1",
                                        header_color="#ffffff")
        self.combined_scatter_card = pn.Card(title='ANI/APSS scatter plot', styles=config.plot_card_style,
                                             header_background="#2e86c1", header_color="#ffffff")
        self.box_plot_card = pn.Card(title='APSS distribution among species', styles=config.plot_card_style,
                                     header_background="#2e86c1", header_color="#ffffff")
        self.box_plot_card_ani = pn.Card(title='ANI distribution among species', styles=config.plot_card_style,
                                         header_background="#2e86c1", header_color="#ffffff")

        download_image_text = 'Save image as: (if no full path, the file is saved under \'Downloads/\')'
        download_table_text = 'Save data table as: (if no full path, the file is saved under \'Downloads/\')'

        # Clustermap elements
        self.clustermap_plot = ""
        self.clustermap_pane = ""
        self.scores_matrix = pd.DataFrame()
        self.clustermap_cmap = pn.widgets.Select(value=config.clustermap_colormaps_list[0],
                                                 options=config.clustermap_colormaps_list,
                                                 name="Select colormap from the following list:")
        self.clustermap_method = pn.widgets.Select(value=config.clustering_methods[0],
                                                   options=config.clustering_methods,
                                                   name="Select a distance metric for the clustering:")
        self.clustermap_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                         options=config.matplotlib_file_formats,
                                                         name="Select image format:")
        self.save_clustermap_file_path = pn.widgets.TextInput(name=download_image_text)
        self.download_clustermap_column = pn.Column()
        download_matrix_text = 'Save matrix as: (if no full path, the file is saved under \'Downloads/\')'
        self.save_matrix_file_path = pn.widgets.TextInput(name=download_matrix_text)
        self.download_matrix_column = pn.Column()
        self.use_metadata_clustermap = pn.widgets.Checkbox(name='Use metadata for coloring', value=False,
                                                           styles={'font-size': "14px"})
        self.metadata_clustermap_card = pn.Card(title='', header_background="#ffffff",
                                                styles={'background': "#ffffff", 'margin': "10px", 'width': "335px"},
                                                hide_header=True,
                                                collapsed=pn.bind(change_disabled_state_inverse,
                                                                  chkbox_state=self.use_metadata_clustermap,
                                                                  watch=True))
        self.color_by_feature = pn.widgets.Select(options=['Select feature'], name="Color rows by:", width=150)
        self.is_continuous_clustermap = pn.widgets.Checkbox(name='Continuous feature', value=False)
        self.feature_colormap = pn.widgets.ColorMap(name="Select colormap:",
                                                    options=config.categorical_colormap_dict,
                                                    value=config.categorical_colormap_dict['cet_glasbey'])
        self.custom_colormap_input_clustermap = pn.widgets.TextInput(
            name='Custom colormap: enter a list of colors separated by comma:',
            placeholder='color1, color2, color3, etc...',
            disabled=pn.bind(change_disabled_state_custom_colormap, value=self.feature_colormap, watch=True))

        # Clustermap ANI elements
        self.clustermap_plot_ani = ""
        self.clustermap_pane_ani = ""
        self.scores_matrix_ani = pd.DataFrame()
        self.clustermap_cmap_ani = pn.widgets.Select(value=config.clustermap_colormaps_list[0],
                                                     options=config.clustermap_colormaps_list,
                                                     name="Select colormap from the following list:")
        self.clustermap_method_ani = pn.widgets.Select(value=config.clustering_methods[0],
                                                       options=config.clustering_methods,
                                                       name="Select a distance metric for the clustering:")
        self.clustermap_image_format_ani = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                             options=config.matplotlib_file_formats,
                                                             name="Select image format:")
        self.save_clustermap_file_path_ani = pn.widgets.TextInput(name=download_image_text)
        self.download_clustermap_column_ani = pn.Column()
        self.save_matrix_file_path_ani = pn.widgets.TextInput(name=download_matrix_text)
        self.download_matrix_column_ani = pn.Column()
        self.use_metadata_clustermap_ani = pn.widgets.Checkbox(name='Use metadata for coloring', value=False,
                                                               styles={'font-size': "14px"})
        self.metadata_clustermap_card_ani = pn.Card(title='', header_background="#ffffff",
                                                    styles={'background': "#ffffff", 'margin': "10px", 'width': "335px"},
                                                    hide_header=True,
                                                    collapsed=pn.bind(change_disabled_state_inverse,
                                                                      chkbox_state=self.use_metadata_clustermap_ani,
                                                                      watch=True))
        self.color_by_feature_ani = pn.widgets.Select(options=['Select feature'], name="Color rows by:", width=150)
        self.is_continuous_clustermap_ani = pn.widgets.Checkbox(name='Continuous feature', value=False)
        self.feature_colormap_ani = pn.widgets.ColorMap(name="Select colormap:",
                                                        options=config.categorical_colormap_dict,
                                                        value=config.categorical_colormap_dict['cet_glasbey'])
        self.custom_colormap_input_clustermap_ani = pn.widgets.TextInput(
            name='Custom colormap: enter a list of colors separated by comma:',
            placeholder='color1, color2, color3, etc...',
            disabled=pn.bind(change_disabled_state_custom_colormap, value=self.feature_colormap_ani, watch=True))

        # Jitter plot elements
        self.jitter_plot = ""
        self.jitter_pane = ""
        self.df_for_jitter = pd.DataFrame()
        self.use_metadata_jitter = pn.widgets.Checkbox(name='Use metadata in plot', value=False,
                                                       styles={'font-size': "14px"})
        self.jitter_color = pn.widgets.ColorPicker(name='Select color', value='#3b89be',
                                                   disabled=pn.bind(change_disabled_state_straight,
                                                                    chkbox_state=self.use_metadata_jitter,
                                                                    watch=True))
        self.jitter_type_select = pn.widgets.Select(options=config.catplot_types, width=150,
                                                    name="Select plot type:")
        self.metadata_jitter_card = pn.Card(styles={'background': "#ffffff", 'margin': "10px", 'width': "300px"},
                                            hide_header=True, collapsed=pn.bind(change_disabled_state_inverse,
                                                                                chkbox_state=self.use_metadata_jitter,
                                                                                watch=True))
        self.jitter_feature_select = pn.widgets.Select(options=['Select feature'], width=150,
                                                       name="Separate plot by following feature:",
                                                       disabled=pn.bind(change_disabled_state_inverse,
                                                                        chkbox_state=self.use_metadata_jitter,
                                                                        watch=True))
        self.jitter_same_color = pn.widgets.ColorPicker(name='Same color:', value=config.same_color,
                                                        disabled=pn.bind(change_disabled_state_inverse,
                                                                         chkbox_state=self.use_metadata_jitter,
                                                                         watch=True))
        self.jitter_different_color = pn.widgets.ColorPicker(name='Different color:', value=config.diff_color,
                                                             disabled=pn.bind(change_disabled_state_inverse,
                                                                              chkbox_state=self.use_metadata_jitter,
                                                                              watch=True))
        self.jitter_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                     options=config.matplotlib_file_formats, name="Select image format:")
        self.save_jitter_file_path = pn.widgets.TextInput(name=download_image_text)
        self.download_jitter_column = pn.Column()
        self.save_jitter_table_path = pn.widgets.TextInput(name=download_table_text)
        self.download_jitter_table_column = pn.Column()

        # Jitter plot ANI elements
        self.jitter_plot_ani = ""
        self.jitter_pane_ani = ""
        self.df_for_jitter_ani = pd.DataFrame()
        self.use_metadata_jitter_ani = pn.widgets.Checkbox(name='Use metadata in plot', value=False,
                                                           styles={'font-size': "14px"})
        self.jitter_color_ani = pn.widgets.ColorPicker(name='Select color', value='#3b89be',
                                                       disabled=pn.bind(change_disabled_state_straight,
                                                                        chkbox_state=self.use_metadata_jitter_ani,
                                                                        watch=True))
        self.jitter_type_select_ani = pn.widgets.Select(options=config.catplot_types, width=150,
                                                        name="Select plot type:")
        self.metadata_jitter_card_ani = pn.Card(styles={'background': "#ffffff", 'margin': "10px", 'width': "300px"},
                                                hide_header=True,
                                                collapsed=pn.bind(change_disabled_state_inverse,
                                                                  chkbox_state=self.use_metadata_jitter_ani,
                                                                  watch=True))
        self.jitter_feature_select_ani = pn.widgets.Select(options=['Select feature'], width=150,
                                                           name="Separate plot by following feature:",
                                                           disabled=pn.bind(change_disabled_state_inverse,
                                                                            chkbox_state=self.use_metadata_jitter_ani,
                                                                            watch=True))
        self.jitter_same_color_ani = pn.widgets.ColorPicker(name='Same color:', value=config.same_color,
                                                            disabled=pn.bind(change_disabled_state_inverse,
                                                                            chkbox_state=self.use_metadata_jitter_ani,
                                                                            watch=True))
        self.jitter_different_color_ani = pn.widgets.ColorPicker(name='Different color:', value=config.diff_color,
                                                                 disabled=pn.bind(
                                                                     change_disabled_state_inverse,
                                                                     chkbox_state=self.use_metadata_jitter_ani,
                                                                     watch=True))
        self.jitter_image_format_ani = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                         options=config.matplotlib_file_formats,
                                                         name="Select image format:")
        self.save_jitter_file_path_ani = pn.widgets.TextInput(name=download_image_text)
        self.download_jitter_column_ani = pn.Column()
        self.save_jitter_table_path_ani = pn.widgets.TextInput(name=download_table_text)
        self.download_jitter_table_column_ani = pn.Column()

        # Combined scatter plot elements
        self.sample_sizes_combined_scatter_slider = pn.widgets.DiscreteSlider(name='Subsampled regions',
                                                                              options=config.sampling_sizes,
                                                                              bar_color='white')
        self.combined_scatter_plot = ""
        self.df_for_combined_scatter = pd.DataFrame()
        self.use_metadata_combined_scatter = pn.widgets.Checkbox(name='Use metadata in plot', value=False,
                                                                 styles={'font-size': "14px"})
        self.combined_scatter_color = pn.widgets.ColorPicker(name='Select color', value='#3b89be',
                                                             disabled=pn.bind(change_disabled_state_straight,
                                                                              chkbox_state=self.use_metadata_combined_scatter,
                                                                              watch=True))
        self.metadata_combined_scatter_card = pn.Card(styles={'background': "#ffffff", 'margin': "10px", 'width': "300px"},
                                                      hide_header=True,
                                                      collapsed=pn.bind(change_disabled_state_inverse,
                                                                        chkbox_state=self.use_metadata_combined_scatter,
                                                                        watch=True))
        self.combined_scatter_feature_select = pn.widgets.Select(options=['Select feature'], width=150,
                                                                 name="Separate plot by following feature:",
                                                                 disabled=pn.bind(change_disabled_state_inverse,
                                                                                  chkbox_state=self.use_metadata_combined_scatter,
                                                                                  watch=True))
        self.combined_scatter_same_color = pn.widgets.ColorPicker(name='Same color:', value=config.same_color,
                                                                  disabled=pn.bind(change_disabled_state_inverse,
                                                                                   chkbox_state=self.use_metadata_combined_scatter,
                                                                                   watch=True))
        self.combined_scatter_different_color = pn.widgets.ColorPicker(name='Different color:', value=config.diff_color,
                                                                       disabled=pn.bind(change_disabled_state_inverse,
                                                                                        chkbox_state=self.use_metadata_combined_scatter,
                                                                                        watch=True))
        self.combined_scatter_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                               options=config.matplotlib_file_formats,
                                                               name="Select image format:")
        self.save_combined_scatter_file_path = pn.widgets.TextInput(name=download_image_text)
        self.download_combined_scatter_column = pn.Column()
        self.save_combined_scatter_table_path = pn.widgets.TextInput(name=download_table_text)
        self.download_combined_scatter_table_column = pn.Column()

        # Network plot elements
        self.network = ""
        self.network_plot_hv = ""
        self.network_plot_matplotlib = ""
        self.df_for_network = pd.DataFrame()
        self.use_metadata_network = pn.widgets.Checkbox(name='Or, use metadata for coloring', value=False,
                                                        styles={'font-size': "14px"})
        self.layout_parameters_card = pn.Card(title='Layout parameters', header_background="#ffffff",
                                              styles={'background': "#ffffff", 'margin': "5px 0 5px 10px",
                                                      'width': "350px", 'padding': "7px"})
        self.metadata_colorby_card = pn.Card(title='Set the coloring by metadata', header_background="#ffffff",
                                             styles={'background': "#ffffff", 'margin': "5px 0 5px 10px",
                                                     'width': "350px"},
                                             hide_header=True, collapsed=pn.bind(change_disabled_state_inverse,
                                                                                 chkbox_state=self.use_metadata_network,
                                                                                 watch=True))
        self.network_node_color = pn.widgets.ColorPicker(name='Nodes color:', value='#459ED9',
                                                         disabled=pn.bind(change_disabled_state_straight,
                                                                          chkbox_state=self.use_metadata_network,
                                                                          watch=True))
        self.network_edge_color = pn.widgets.ColorPicker(name='Edges color:', value='#000000',
                                                         disabled=pn.bind(change_disabled_state_straight,
                                                                          chkbox_state=self.use_metadata_network,
                                                                          watch=True))
        self.nodes_color_by = pn.widgets.Select(options=['Select feature'], name="Color nodes by:", width=130)
        self.is_continuous_network = pn.widgets.Checkbox(name='Continuous feature', value=False)
        self.nodes_colormap = pn.widgets.ColorMap(name="Select colormap for nodes:",
                                                  options=config.categorical_colormap_dict,
                                                  value=config.categorical_colormap_dict['cet_glasbey'],
                                                  styles={'margin': "5px 0 0 15px"})
        self.custom_colormap_input = pn.widgets.TextInput(name='Custom colormap: enter a list of colors separated by '
                                                               'comma',
                                                          placeholder='color1, color2, color3, etc...',
                                                          disabled=pn.bind(change_disabled_state_custom_colormap,
                                                                           value=self.nodes_colormap,
                                                                           watch=True),
                                                          styles={'margin': "5px 0 0 15px"}
                                                          )
        self.highlight_nodes_by_feature = pn.widgets.Checkbox(name='Highlight nodes by feature', value=False,
                                                              styles={'margin': "0 0 0 15px"})
        self.nodes_highlight_by = pn.widgets.Select(options=['Select feature'], name="Highlight nodes by:", width=140,
                                                    disabled=pn.bind(change_disabled_state_inverse,
                                                                     chkbox_state=self.highlight_nodes_by_feature,
                                                                     watch=True)
                                                    )
        self.nodes_highlight_group_select = pn.widgets.Select(options=['Select group'], name="Select group:", width=140,
                                                              disabled=pn.bind(change_disabled_state_inverse,
                                                              chkbox_state=self.highlight_nodes_by_feature,
                                                              watch=True)
                                                             )
        self.color_edges_by_feature = pn.widgets.Checkbox(name='Color edges by feature (same/different)', value=False,
                                                          styles={'margin': "0 0 0 15px"})
        self.edges_color_by = pn.widgets.Select(options=['Select feature'],
                                                name="Color edges by:", width=100,
                                                disabled=pn.bind(change_disabled_state_inverse,
                                                                 chkbox_state=self.color_edges_by_feature,
                                                                 watch=True))
        self.network_within_color = pn.widgets.ColorPicker(name='Same color:', value='#000000',
                                                           disabled=pn.bind(change_disabled_state_inverse,
                                                                            chkbox_state=self.color_edges_by_feature,
                                                                            watch=True))
        self.network_between_color = pn.widgets.ColorPicker(name='Different color:', value='#000000',
                                                            disabled=pn.bind(change_disabled_state_inverse,
                                                                             chkbox_state=self.color_edges_by_feature,
                                                                             watch=True)
                                                            )
        self.show_labels_chkbox = pn.widgets.Checkbox(name='Show sample names', value=False,
                                                      styles={'font-size': "14px"})
        self.all_or_highlighted_select = pn.widgets.Select(options=['All', 'Highlighted only'], width=140,
                                                           disabled=pn.bind(change_disabled_state_inverse,
                                                                            chkbox_state=self.show_labels_chkbox,
                                                                            watch=True))
        self.highlight_sample_chkbox = pn.widgets.Checkbox(name='Highlight node(s):', value=False,
                                                           styles={'font-size': "14px", 'padding-top': "6px"})
        self.highlight_sample_input = pn.widgets.TextInput(
                                                           placeholder='Enter sampleID(s)', width=205,
                                                           disabled=pn.bind(change_disabled_state_inverse,
                                                           chkbox_state=self.highlight_sample_chkbox, watch=True),
                                                          )
        self.network_threshold_select = pn.widgets.Select(name="Threshold for network connections:", width=200,
                                                          options=[])
        self.network_threshold_input = pn.widgets.FloatInput(name='Define threshold:',
                                                             value=config.APSS_connections_threshold_default, step=0.01,
                                                             start=0.5, end=1.0, width=100,
                                                             disabled=True
                                                             )
        self.network_iterations = pn.widgets.DiscreteSlider(name='Number of iterations',
                                                            options=config.network_iterations_options,
                                                            bar_color='white')
        self.network_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                      options=config.matplotlib_file_formats,
                                                      name="Select image format:")
        self.save_network_plot_path = pn.widgets.TextInput(name=download_image_text)
        self.download_network_column = pn.Column()
        self.save_network_table_path = pn.widgets.TextInput(name=download_table_text)
        self.download_network_table_column = pn.Column()
        self.network_pane = ""
        self.pos_dict = dict()
        self.nodes_list = []
        self.network = ""
        self.APSS_connections_threshold = config.APSS_connections_threshold_default

        # Network plot ANI elements
        self.network_ani = ""
        self.network_plot_hv_ani = ""
        self.network_plot_matplotlib_ani = ""
        self.df_for_network_ani = pd.DataFrame()
        self.use_metadata_network_ani = pn.widgets.Checkbox(name='Or, use metadata for coloring', value=False,
                                                            styles={'font-size': "14px"})
        self.layout_parameters_card_ani = pn.Card(title='Layout parameters', header_background="#ffffff",
                                                  styles={'background': "#ffffff", 'margin': "5px 0 5px 10px",
                                                          'width': "350px", 'padding': "7px"})
        self.metadata_colorby_card_ani = pn.Card(title='Set the coloring by metadata', header_background="#ffffff",
                                                 styles={'background': "#ffffff", 'margin': "5px 0 5px 10px",
                                                         'width': "350px"}, hide_header=True,
                                                 collapsed=pn.bind(change_disabled_state_inverse,
                                                                   chkbox_state=self.use_metadata_network_ani,
                                                                   watch=True))
        self.network_node_color_ani = pn.widgets.ColorPicker(name='Nodes color:', value='#459ED9',
                                                             disabled=pn.bind(change_disabled_state_straight,
                                                                              chkbox_state=self.use_metadata_network_ani,
                                                                              watch=True))
        self.network_edge_color_ani = pn.widgets.ColorPicker(name='Edges color:', value='#000000',
                                                             disabled=pn.bind(change_disabled_state_straight,
                                                                              chkbox_state=self.use_metadata_network_ani,
                                                                              watch=True))
        self.nodes_color_by_ani = pn.widgets.Select(options=['Select feature'], name="Color nodes by:", width=130)
        self.is_continuous_network_ani = pn.widgets.Checkbox(name='Continuous feature', value=False)
        self.nodes_colormap_ani = pn.widgets.ColorMap(name="Select colormap for nodes:",
                                                      options=config.categorical_colormap_dict,
                                                      value=config.categorical_colormap_dict['cet_glasbey'],
                                                      styles={'margin': "5px 0 0 15px"})
        self.custom_colormap_input_ani = pn.widgets.TextInput(name='Custom colormap: enter a list of colors separated '
                                                                   'by comma',
                                                              placeholder='color1, color2, color3, etc...',
                                                              disabled=pn.bind(change_disabled_state_custom_colormap,
                                                                               value=self.nodes_colormap_ani,
                                                                               watch=True),
                                                              styles={'margin': "5px 0 0 15px"}
                                                             )
        self.highlight_nodes_by_feature_ani = pn.widgets.Checkbox(name='Highlight nodes by feature', value=False,
                                                                  styles={'margin': "0 0 0 15px"})
        self.nodes_highlight_by_ani = pn.widgets.Select(options=['Select feature'], name="Highlight nodes by:", width=140,
                                                        disabled=pn.bind(change_disabled_state_inverse,
                                                                         chkbox_state=self.highlight_nodes_by_feature_ani,
                                                                         watch=True)
                                                        )
        self.nodes_highlight_group_select_ani = pn.widgets.Select(options=['Select group'], name="Select group:",
                                                                  width=140,
                                                                  disabled=pn.bind(change_disabled_state_inverse,
                                                                                   chkbox_state=self.highlight_nodes_by_feature_ani,
                                                                                   watch=True)
                                                                 )
        self.color_edges_by_feature_ani = pn.widgets.Checkbox(name='Color edges by feature (same/different)',
                                                              value=False, styles={'margin': "0 0 0 15px"})
        self.edges_color_by_ani = pn.widgets.Select(options=['Select feature'],
                                                    name="Color edges by:", width=100,
                                                    disabled=pn.bind(change_disabled_state_inverse,
                                                                     chkbox_state=self.color_edges_by_feature_ani,
                                                                     watch=True))
        self.network_within_color_ani = pn.widgets.ColorPicker(name='Same color:', value='#000000',
                                                               disabled=pn.bind(change_disabled_state_inverse,
                                                                                chkbox_state=self.color_edges_by_feature_ani,
                                                                                watch=True))
        self.network_between_color_ani = pn.widgets.ColorPicker(name='Different color:', value='#000000',
                                                                disabled=pn.bind(change_disabled_state_inverse,
                                                                                 chkbox_state=self.color_edges_by_feature_ani,
                                                                                 watch=True)
                                                            )
        self.show_labels_chkbox_ani = pn.widgets.Checkbox(name='Show sample names', value=False,
                                                          styles={'font-size': "14px"})
        self.all_or_highlighted_select_ani = pn.widgets.Select(options=['All', 'Highlighted only'], width=140,
                                                               disabled=pn.bind(change_disabled_state_inverse,
                                                                                chkbox_state=self.show_labels_chkbox_ani,
                                                                                watch=True))
        self.highlight_sample_chkbox_ani = pn.widgets.Checkbox(name='Highlight node(s):', value=False,
                                                               styles={'font-size': "14px", 'padding-top': "6px"})
        self.highlight_sample_input_ani = pn.widgets.TextInput(placeholder='Enter sampleID(s)', width=205,
                                                               disabled=pn.bind(change_disabled_state_inverse,
                                                                                chkbox_state=self.highlight_sample_chkbox_ani,
                                                                                watch=True),
                                                               )
        self.network_threshold_select_ani = pn.widgets.Select(name="Threshold for network connections:", width=200,
                                                              options=[])
        self.network_threshold_input_ani = pn.widgets.FloatInput(name='Define threshold:',
                                                                 value=config.ANI_connections_threshold_default,
                                                                 step=0.001, start=0.800, end=0.999, width=100,
                                                                 disabled=True
                                                                 )
        self.network_iterations_ani = pn.widgets.DiscreteSlider(name='Number of iterations',
                                                                options=config.network_iterations_options,
                                                                bar_color='white')
        self.network_image_format_ani = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                          options=config.matplotlib_file_formats,
                                                          name="Select image format:")
        self.save_network_plot_path_ani = pn.widgets.TextInput(name=download_image_text)
        self.download_network_column_ani = pn.Column()
        self.save_network_table_path_ani = pn.widgets.TextInput(name=download_table_text)
        self.download_network_table_column_ani = pn.Column()
        self.network_pane_ani = ""
        self.pos_dict_ani = dict()
        self.nodes_list_ani = []
        self.network_ani = ""
        self.ani_connections_threshold = config.ANI_connections_threshold_default

        # Box-plot elements
        self.box_plot = ""
        self.box_plot_pane = ""
        self.use_metadata_box_plot = pn.widgets.Checkbox(name='Use metadata in plot', value=False,
                                                         styles={'font-size': "14px"})
        self.box_plot_color = pn.widgets.ColorPicker(name='Select color', value=config.diff_color,
                                                     disabled=pn.bind(change_disabled_state_straight,
                                                                      chkbox_state=self.use_metadata_box_plot,
                                                                      watch=True))
        self.metadata_box_plot_card = pn.Card(styles={'background': "#ffffff", 'margin': "10px", 'width': "300px"},
                                              hide_header=True, collapsed=pn.bind(change_disabled_state_inverse,
                                                                                  chkbox_state=self.use_metadata_box_plot,
                                                                                  watch=True))
        self.box_plot_feature_select = pn.widgets.Select(options=['Select feature'], width=150,
                                                         name="Separate plot by following feature:",
                                                         disabled=pn.bind(change_disabled_state_inverse,
                                                                          chkbox_state=self.use_metadata_box_plot,
                                                                          watch=True))
        self.box_plot_same_color = pn.widgets.ColorPicker(name='Same color:', value=config.same_color,
                                                          disabled=pn.bind(change_disabled_state_inverse,
                                                                           chkbox_state=self.use_metadata_box_plot,
                                                                           watch=True))
        self.box_plot_different_color = pn.widgets.ColorPicker(name='Different color:', value=config.diff_color,
                                                               disabled=pn.bind(change_disabled_state_inverse,
                                                                                chkbox_state=self.use_metadata_box_plot,
                                                                                watch=True))
        self.box_plot_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                       options=config.matplotlib_file_formats,
                                                       name="Select image format:")
        self.save_box_plot_file_path = pn.widgets.TextInput(name=download_image_text)
        self.download_box_plot_column = pn.Column()
        self.save_boxplot_table_path = pn.widgets.TextInput(name=download_table_text)
        self.download_boxplot_table_column = pn.Column()
        self.save_pvalues_table_path = pn.widgets.TextInput(
            name='Save P-values table as: (if no full path, the file is saved under \'Downloads/\')')
        self.download_pvalues_table_column = pn.Column()
        self.download_pvalues_button = pn.widgets.Button(name='Download P-values table in tsv format',
                                                         button_type='primary',
                                                         disabled=pn.bind(change_disabled_state_inverse,
                                                                          chkbox_state=self.use_metadata_box_plot,
                                                                          watch=True))
        self.download_pvalues_button.on_click(self.download_pvalues_table)
        self.download_multi_col = pn.Column

        # Box-plot ANI elements
        self.box_plot_ani = ""
        self.box_plot_pane_ani = ""
        self.use_metadata_box_plot_ani = pn.widgets.Checkbox(name='Use metadata in plot', value=False,
                                                             styles={'font-size': "14px"})
        self.box_plot_color_ani = pn.widgets.ColorPicker(name='Select color', value=config.diff_color,
                                                         disabled=pn.bind(change_disabled_state_straight,
                                                                          chkbox_state=self.use_metadata_box_plot_ani,
                                                                          watch=True))
        self.metadata_box_plot_card_ani = pn.Card(styles={'background': "#ffffff", 'margin': "10px", 'width': "300px"},
                                                  hide_header=True, collapsed=pn.bind(change_disabled_state_inverse,
                                                                                      chkbox_state=self.use_metadata_box_plot_ani,
                                                                                      watch=True))
        self.box_plot_feature_select_ani = pn.widgets.Select(options=['Select feature'], width=150,
                                                             name="Separate plot by following feature:",
                                                             disabled=pn.bind(change_disabled_state_inverse,
                                                                              chkbox_state=self.use_metadata_box_plot_ani,
                                                                              watch=True))
        self.box_plot_same_color_ani = pn.widgets.ColorPicker(name='Same color:', value=config.same_color,
                                                              disabled=pn.bind(change_disabled_state_inverse,
                                                                               chkbox_state=self.use_metadata_box_plot_ani,
                                                                               watch=True))
        self.box_plot_different_color_ani = pn.widgets.ColorPicker(name='Different color:', value=config.diff_color,
                                                                   disabled=pn.bind(change_disabled_state_inverse,
                                                                                    chkbox_state=self.use_metadata_box_plot_ani,
                                                                                    watch=True))
        self.box_plot_image_format_ani = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                           options=config.matplotlib_file_formats,
                                                           name="Select image format:")
        self.save_box_plot_file_path_ani = pn.widgets.TextInput(name=download_image_text)
        self.download_box_plot_column_ani = pn.Column()
        self.save_boxplot_table_path_ani = pn.widgets.TextInput(name=download_table_text)
        self.download_boxplot_table_column_ani = pn.Column()
        self.save_pvalues_table_path_ani = pn.widgets.TextInput(
            name='Save P-values table as: (if no full path, the file is saved under \'Downloads/\')')
        self.download_pvalues_table_column_ani = pn.Column()
        self.download_multi_col_ani = pn.Column()
        self.download_pvalues_button_ani = pn.widgets.Button(name='Download P-values table in tsv format',
                                                             button_type='primary',
                                                             disabled=pn.bind(change_disabled_state_inverse,
                                                                              chkbox_state=self.use_metadata_box_plot_ani,
                                                                              watch=True))
        self.download_pvalues_button_ani.on_click(self.download_pvalues_table_ani)

        # synteny_per_pos plots elements
        self.visited_synteny_per_pos_tab = 0
        self.finished_initial_synteny_per_pos_plot = 0
        self.ax_for_synteny_per_pos_plot = ""
        self.ax_annotations = ""
        self.synteny_per_pos_plot = ""
        self.coverage_plot = ""
        self.line_avg_plot = ""
        self.hypervar_bars = ""
        self.hypercons_bars = ""
        self.synteny_per_pos_pane = ""
        self.contig_select = pn.widgets.Select(name='Select a contig:', options=[], styles={'margin': "0"})
        self.contig_select_watcher = ""
        self.contig_name = ""
        self.contig_length = ""
        self.filter_plot_by_metadata = 0
        self.sorting_select = pn.widgets.Select(name='Sort by:', options=config.contig_sorting_options,
                                                styles={'margin': "0"})
        self.sorting_select_watcher = ""
        self.start_pos_input = pn.widgets.TextInput(name='Start position', styles={'font-size': "14px"})
        self.end_pos_input = pn.widgets.TextInput(name='End position', styles={'font-size': "14px"})
        self.avg_plot_chkbox = pn.widgets.Checkbox(name='Show average synteny scores', value=True,
                                                   styles={'font-size': "14px"})
        self.avg_plot_chkbox_watcher = ""
        self.avg_plot_color = pn.widgets.ColorPicker(name='Color:', value='#000080',
                                                     disabled=pn.bind(change_disabled_state_inverse,
                                                                      chkbox_state=self.avg_plot_chkbox,
                                                                      watch=True))
        self.avg_plot_color_watcher = ""
        self.coverage_plot_chkbox = pn.widgets.Checkbox(name='Show all synteny scores', value=True,
                                                        styles={'font-size': "14px"})
        self.coverage_plot_chkbox_watcher = ""
        self.coverage_plot_color = pn.widgets.ColorPicker(name='Color:', value='#ba55d3',
                                                          disabled=pn.bind(change_disabled_state_inverse,
                                                                           chkbox_state=self.coverage_plot_chkbox,
                                                                           watch=True))
        self.coverage_plot_color_watcher = ""
        self.hypervar_chkbox = pn.widgets.Checkbox(name='Highlight hypervariable regions', value=True,
                                                   styles={'font-size': "14px"})
        self.hypervar_chkbox_watcher = ""
        self.hypervar_color = pn.widgets.ColorPicker(name='Color:', value=config.variable_color,
                                                     disabled=pn.bind(change_disabled_state_inverse,
                                                                      chkbox_state=self.hypervar_chkbox,
                                                                      watch=True))
        self.hypervar_color_watcher = ""
        self.hypercons_chkbox = pn.widgets.Checkbox(name='Highlight hyperconserved regions', value=True,
                                                    styles={'font-size': "14px"})
        self.hypercons_chkbox_watcher = ""
        self.hypercons_color = pn.widgets.ColorPicker(name='Color:', value=config.conserved_color,
                                                      disabled=pn.bind(change_disabled_state_inverse,
                                                                       chkbox_state=self.hypercons_chkbox,
                                                                       watch=True))
        self.hypercons_color_watcher = ""
        self.alpha_slider = pn.widgets.FloatSlider(name='Alpha transparency', start=0, end=1, step=0.1, value=0.3,
                                                   styles={'font-size': "14px"})
        self.alpha_slider_watcher = ""
        self.synteny_per_pos_image_format = pn.widgets.Select(value=config.matplotlib_file_formats[0],
                                                       options=config.matplotlib_file_formats,
                                                       name="Select image format:")
        self.save_synteny_per_pos_plot_path = pn.widgets.TextInput(name=download_image_text)
        self.download_synteny_per_pos_plot_column = pn.Column()
        self.save_synteny_per_pos_table_path = pn.widgets.TextInput(name=download_table_text)
        self.download_synteny_per_pos_table_column = pn.Column()

        self.filter_by_metadata_card = pn.Card(title='Filter plot by metadata', collapsed=True,
                                               styles={'margin': "5px 0 5px 10px", 'width': "1050px"})
        self.synteny_per_pos_feature_select = pn.widgets.Select(options=['Select feature'], width=200,
                                                                name="Filter plot by the following feature:")
        self.synteny_per_pos_groups_select = pn.widgets.MultiSelect(options=[], width=350, height=200, size=5,
                                                                    name='Include the following groups in the plot:')
        self.filter_synteny_per_pos_button = pn.widgets.Button(name='Filter plot', button_type='primary')
        self.filter_synteny_per_pos_button.on_click(self.filter_synteny_per_pos_plot)
        self.reset_filter_synteny_per_pos_button = pn.widgets.Button(name='Reset filteration', button_type='primary')
        self.reset_filter_synteny_per_pos_button.on_click(self.reset_filter_synteny_per_pos_plot)

        self.add_annotation_card = pn.Card(title='Add annotation data', collapsed=True,
                                           styles={'margin': "5px 0 5px 10px", 'width': "1050px"})
        self.upload_annotation_file = pn.widgets.FileInput(accept='.gff')
        self.gff_filename = ""
        self.annotation_upload_col = pn.Column()
        self.show_annotations_col = pn.Column()
        self.show_annotations_row = pn.Row()
        self.submit_annotation_file_button = pn.widgets.Button(name='Upload file', button_type='primary')
        self.submit_annotation_file_button.on_click(self.read_annotation_file)
        self.show_annotations_chkbox = pn.widgets.Checkbox(name='Show annotated genes for current contig', value=False,
                                                           styles={'font-size': "14px"})
        self.annotation_help_button = pn.widgets.ButtonIcon(icon="help", size="1.5em",
                                                            description="Genes plot help", margin=(3, 0, 0, 0))
        self.annotation_help_button.on_click(self.show_annotations_help_float_panel)
        self.show_annotations_watcher = ""
        self.show_annotations = False

        # Build the initial layout
        mandatory_input_title = "Mandatory input"
        self.main_area.append(pn.pane.Markdown(mandatory_input_title,
                                               styles={'font-size': "24px", 'color': config.title_purple_color,
                                                       'margin-bottom': "0", 'padding-bottom': "0"}))

        input_type_title = "Select the type of input to analyse using StrainVis:"
        self.main_area.append(pn.pane.Markdown(input_type_title, styles={'font-size': "18px", 'margin-bottom': "0",
                                                                         'margin-top': "0", 'padding-top': "0"}))
        self.main_area.append(self.input_type_radio_group)

        file_input_title = "Upload SynTracker's output table 'synteny_scores_per_region.csv' for one or multiple " \
                           "species:"
        text_input_title = "Or, if the file size is bigger than 300 Mb, enter it's full path here:"
        self.mandatory_input_card.append(pn.pane.Markdown(file_input_title, styles={'font-size': "16px",
                                                                                    'margin-bottom': "0",
                                                                                    'margin-top': "0"}))
        self.mandatory_input_card.append(self.SynTracker_input_file)
        self.mandatory_input_card.append(pn.pane.Markdown(text_input_title, styles={'font-size': "16px",
                                                                                    'margin-bottom': "0",
                                                                                    'margin-top': "0"}))
        self.mandatory_input_card.append(self.SynTracker_text_input)
        self.main_area.append(self.mandatory_input_card)

        self.main_area.append(pn.Spacer(height=20))

        optional_input_title = "Optional input"
        self.main_area.append(pn.pane.Markdown(optional_input_title,
                                               styles={'font-size': "24px", 'color': config.title_purple_color,
                                                       'margin-bottom': "0"}))

        # Metadata input
        metadata_upload_title = "Upload a metadata file for the compared samples (in tab delimited format):"
        metadata_help_button = pn.widgets.ButtonIcon(icon="help", size="1.5em", description="Metadata help",
                                                     margin=(17, 5, 0, 0))
        metadata_help_button.on_click(self.show_metadata_help_float_panel)
        self.metadata_upload_row = pn.Row(pn.pane.Markdown(metadata_upload_title, styles={'font-size': "16px",
                                                                                          'margin-bottom': "0",
                                                                                          'margin-top': "0"}),
                                          metadata_help_button, pn.Spacer(width=1))
        self.optional_input_card.append(self.metadata_upload_row)
        self.optional_input_card.append(self.metadata_file)
        self.main_area.append(self.optional_input_card)

        self.main_area.append(pn.Spacer(height=30))

        self.main_area.append(self.submit_button)

        pn.state.on_session_destroyed(destroyed)

        self.template.servable()

    def show_metadata_help_float_panel(self, event):
        metadata_note = "The metadata file may contain an unlimited number of columns (features).  " \
                        "\nThe first column must contain the sample IDs (identical to the sample IDs that appear in " \
                        "the input file(s))."
        markdown = pn.pane.Markdown(metadata_note, styles={'font-size': "14px", 'color': config.title_red_color})
        floatpanel = pn.layout.FloatPanel(markdown, name='metadata file help', margin=10, width=400)
        self.metadata_upload_row[2] = floatpanel

    def load_correct_page(self, event):
        if self.header_buttons.value == 'Home':
            self.load_home_page()
        else:
            self.load_help_page()

    def load_home_page(self):
        self.main_container.clear()
        self.main_container.append(self.main_area)

    def load_help_page(self):
        self.main_container.clear()
        self.main_container.append(self.help_area)

    def changed_main_tab(self, event):
        if self.menu_tabs.active == 1:
            self.main_container.clear()
            self.main_container.append(self.help_area)

        else:
            self.main_container.clear()
            self.main_container.append(self.main_area)

    def update_input_card(self, event):
        syn_file_input_title = "Upload SynTracker's output table 'synteny_scores_per_region.csv' for one or multiple " \
                               "species:"
        ani_file_input_title = "Upload tab-delimited ANI file for one or multiple species:"
        ani_help_button = pn.widgets.ButtonIcon(icon="help", size="1.5em", description="ANI help",
                                                margin=(17, 5, 0, 0))
        ani_help_button.on_click(self.show_ani_help_float_panel)
        text_input_title = "Or, if the file size is bigger than 300 Mb, enter it's full path here:"

        self.mandatory_input_card.clear()
        self.ani_upload_row.clear()
        self.SynTracker_input_file.value = None
        self.SynTracker_text_input.value = ""
        self.ANI_input_file.value = None
        self.ANI_text_input.value = ""

        # Only SynTracker input
        if self.input_type_radio_group.value == 'SynTracker output file':
            self.mandatory_input_card.append(pn.pane.Markdown(syn_file_input_title, styles={'font-size': "16px",
                                                                                            'margin-bottom': "0",
                                                                                            'margin-top': "0"}))
            self.mandatory_input_card.append(self.SynTracker_input_file)
            self.mandatory_input_card.append(pn.pane.Markdown(text_input_title, styles={'font-size': "16px",
                                                                                        'margin-bottom': "0",
                                                                                        'margin-top': "0"}))
            self.mandatory_input_card.append(self.SynTracker_text_input)

            self.input_mode = "SynTracker"
            self.active_input_mode = "SynTracker"
            self.score_type = "APSS"

        # Only ANI input
        elif self.input_type_radio_group.value == 'ANI file':
            self.ani_upload_row.append(pn.pane.Markdown(ani_file_input_title, styles={'font-size': "16px",
                                                                                      'margin-bottom': "0",
                                                                                      'margin-top': "0"}))
            self.ani_upload_row.append(ani_help_button)
            self.mandatory_input_card.append(self.ani_upload_row)
            self.mandatory_input_card.append(self.ANI_input_file)
            self.mandatory_input_card.append(pn.pane.Markdown(text_input_title, styles={'font-size': "16px",
                                                                                        'margin-bottom': "0",
                                                                                        'margin-top': "0"}))
            self.mandatory_input_card.append(self.ANI_text_input)

            self.input_mode = "ANI"
            self.active_input_mode = "ANI"
            self.score_type = "ANI"

        # Combined input
        else:
            self.mandatory_input_card.append(pn.pane.Markdown(syn_file_input_title, styles={'font-size': "16px",
                                                                                            'margin-bottom': "0",
                                                                                            'margin-top': "0"}))
            self.mandatory_input_card.append(self.SynTracker_input_file)
            self.mandatory_input_card.append(pn.pane.Markdown(text_input_title, styles={'font-size': "16px",
                                                                                        'margin-bottom': "0",
                                                                                        'margin-top': "0"}))
            self.mandatory_input_card.append(self.SynTracker_text_input)

            self.mandatory_input_card.append(pn.Spacer(height=10))
            self.mandatory_input_card.append(pn.layout.Divider())

            self.ani_upload_row.append(pn.pane.Markdown(ani_file_input_title, styles={'font-size': "16px",
                                                                                        'margin-bottom': "0",
                                                                                        'margin-top': "0"}))
            self.ani_upload_row.append(ani_help_button)
            self.mandatory_input_card.append(self.ani_upload_row)
            self.mandatory_input_card.append(self.ANI_input_file)
            self.mandatory_input_card.append(pn.pane.Markdown(text_input_title, styles={'font-size': "16px",
                                                                                        'margin-bottom': "0",
                                                                                        'margin-top': "0"}))
            self.mandatory_input_card.append(self.ANI_text_input)

            matching_note = "Please note: the names of species and the sample IDs must be identical between both " \
                            "input files."
            self.mandatory_input_card.append(pn.pane.Markdown(matching_note, styles={'font-size': "16px",
                                                                                     'margin-bottom': "0",
                                                                                     'margin-top': "0",
                                                                                     'color': config.title_red_color}))

            self.input_mode = "both"
            self.active_input_mode = "SynTracker"
            self.score_type = "APSS"

    def show_ani_help_float_panel(self, event):
        ani_text = "The tab-delimited ANI file should contain the following columns: " \
                   "'Ref_genome', 'Sample1', 'Sample2', 'ANI'."
        markdown = pn.pane.Markdown(ani_text, styles={'font-size': "14px", 'color': config.title_red_color})
        floatpanel = pn.layout.FloatPanel(markdown, name='ANI file help', margin=10, width=450)
        self.ani_upload_row.append(floatpanel)

    def init_parameters(self):
        del self.score_per_region_all_genomes_df
        del self.score_per_region_genomes_subset_df
        del self.score_per_region_selected_genome_df
        del self.ani_scores_all_genomes_df
        del self.ani_scores_genomes_subset_df
        del self.ani_scores_selected_genome_df
        del self.score_per_pos_contig
        del self.score_per_pos_contig_filtered
        del self.genomes_subset_selected_size_APSS_df
        del self.pairs_num_per_sampling_size_multi_genomes_df
        del self.boxplot_p_values_df
        del self.boxplot_p_values_df_ani
        del self.df_for_jitter
        del self.df_for_jitter_ani
        del self.scores_matrix
        del self.scores_matrix_ani
        del self.df_for_network
        del self.avg_score_per_pos_contig
        del self.avg_score_per_pos_contig_filtered
        del self.metadata_dict
        del self.groups_per_feature_dict
        del self.APSS_by_genome_all_sizes_dict
        del self.APSS_all_genomes_all_sizes_dict

        gc.collect()

    def create_new_session(self, event):
        self.init_parameters()
        pn.state.location.reload = True

    def submit_new_file_button(self):
        button_column = pn.Column(pn.Spacer(height=30), self.new_file_button)
        return button_column

    def load_input_file(self, event):
        print("\nIn load_input_file. Input mode = " + self.input_mode)

        # Verify that a SynTracker file was loaded
        if self.input_mode == "SynTracker" or self.input_mode == "both":

            # File path was given - read directly
            if self.SynTracker_text_input.value != "":
                self.is_syntracker_file_path = 1

                # Verify that the input is an existing file path
                if os.path.exists(self.SynTracker_text_input.value):
                    print("\n\nSynTracker input file path: " + self.SynTracker_text_input.value)
                    self.syntracker_filename = os.path.basename(self.SynTracker_text_input.value)
                    print("SynTracker input file name: " + self.syntracker_filename)
                    #self.display_results_page()
                    self.syntracker_loaded = 1

                else:
                    title = "The requested input file does not exist, please enter again a valid file path"
                    self.display_error_page(title)

            # File was given via FileInput widget
            else:
                # Filename is None - usually because of server problems
                if self.SynTracker_input_file.filename is None:
                    print("SynTracker input file name: None")
                    title = "Cannot upload the requested file (probably server problems) - please try again by entering " \
                            "the file's full path"
                    self.display_error_page(title)

                else:
                    self.syntracker_filename = self.SynTracker_input_file.filename
                    content_length = len(self.SynTracker_input_file.value) / 1000
                    print("Content length of SynTracker file in Kb: " + str(content_length))

                    # There is filename but no content - usually happens when the file is too big
                    if content_length == 0:
                        title = "Cannot upload the requested file (probably too big) - please try again by entering " \
                                "the file's full path"
                        self.display_error_page(title)

                    # File has content
                    else:
                        self.syntracker_loaded = 1

        # Verify that an ANI file was loaded
        if self.input_mode == "ANI" or self.input_mode == "both":

            # File path was given - read directly
            if self.ANI_text_input.value != "":
                self.is_ani_file_path = 1

                # Verify that the input is an existing file path
                if os.path.exists(self.ANI_text_input.value):
                    print("\n\nANI input file path: " + self.ANI_text_input.value)
                    self.ani_filename = os.path.basename(self.ANI_text_input.value)
                    print("ANI input file name: " + self.ani_filename)
                    #self.display_results_page()
                    self.ani_loaded = 1

                else:
                    title = "The requested input file does not exist, please enter again a valid file path"
                    self.display_error_page(title)

            # File was given via FileInput widget
            else:
                # Filename is None - usually because of server problems
                if self.ANI_input_file.filename is None:
                    print("ANI input file name: None")
                    title = "Cannot upload the requested file (probably server problems) - " \
                            "please try again by entering the file's full path"
                    self.display_error_page(title)

                else:
                    self.ani_filename = self.ANI_input_file.filename
                    content_length = len(self.ANI_input_file.value) / 1000
                    print("Content length of ANI file in Kb: " + str(content_length))

                    # There is filename but no content - usually happens when the file is too big
                    if content_length == 0:
                        title = "Cannot upload the requested file (probably too big) - please try again by entering " \
                                "the file's full path"
                        self.display_error_page(title)

                    # File has content
                    else:
                        #self.display_results_page()
                        self.ani_loaded = 1

        ## Verify that all necessary files were loaded according to the input mode
        # Input type: SynTracker
        if self.input_mode == "SynTracker":
            if self.syntracker_loaded:
                self.display_results_page()
                self.input_file_loaded = 1

        # Input type: ANI
        elif self.input_mode == "ANI":
            if self.ani_loaded:
                self.display_results_page()
                self.input_file_loaded = 1

        # Input type: both
        else:
            if self.syntracker_loaded and self.ani_loaded:
                self.display_results_page()
                self.input_file_loaded = 1
            else:
                title = "There is a problem uploading both input files - please try again"
                self.display_error_page(title)

        # Check if the user provided a metadata file
        if self.metadata_file.filename is not None:
            content_length = len(self.metadata_file.value) / 1000
            print("\nUploaded metadata file: " + self.metadata_file.filename)
            #print("Content length of metadata file in Kb: " + str(content_length))

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

        self.main_area.append(pn.pane.Markdown(message, styles={'font-size': "20px", 'color': config.title_red_color}))
        self.main_area.append(self.submit_new_file_button())

    def display_results_page(self):
        self.main_area.clear()

        if self.input_mode == "SynTracker":
            title = "Loading input file: " + self.syntracker_filename
        elif self.input_mode == "ANI":
            title = "Loading input file: " + self.ani_filename
        else:
            title = "Loading input files:\n" + self.syntracker_filename + "\n" + self.ani_filename

        self.main_area.append(pn.pane.Markdown(title, styles={'font-size': "20px"}, hard_line_break=True))

    def start_process(self):
        self.ref_genomes_list = []
        self.ref_genomes_list_by_pairs_num = []
        self.ref_genomes_list_syntracker = []
        self.ref_genomes_list_by_pairs_num_syntracker = []
        self.ref_genomes_list_ani = []
        self.ref_genomes_list_by_pairs_num_ani = []
        self.annotation_per_ref_genome_dict = {}
        self.selected_genomes_subset = []
        self.sorted_selected_genomes_subset = []
        self.sorted_selected_genomes_subset_ani = []
        self.selected_subset_species_num = 0
        self.species_num_at_sampling_size = 0
        self.species_num_in_subset_syntracker = 0
        self.species_num_in_subset_ani = 0
        self.single_multi_genome_tabs.clear()
        self.main_single_column.clear()
        self.ref_genome_column.clear()
        self.synteny_single_initial_plots_column.clear()
        self.main_multi_column.clear()
        self.synteny_single_tabs.clear()
        self.synteny_ani_single_tabs.clear()
        self.synteny_ani_multi_tabs.clear()
        self.APSS_analyses_single_column.clear()
        self.synteny_per_pos_plot_column.clear()
        self.plots_by_size_single_column.clear()
        self.ani_single_plots_column.clear()
        self.combined_single_plots_column.clear()
        self.synteny_multi_initial_plots_column.clear()
        self.ani_multi_plots_column.clear()
        self.plots_by_size_multi_column.clear()

        # Read the SynTracker input
        if self.input_mode == "SynTracker" or self.input_mode == "both":

            # Read the file directly from the path
            before = time.time()
            if self.is_syntracker_file_path == 1:
                file = self.SynTracker_text_input.value

            # Read the content of the uploaded file
            else:
                file = io.BytesIO(self.SynTracker_input_file.value)

            self.score_per_region_all_genomes_df = pd.read_csv(file, usecols=lambda c: c in set(config.col_set))
            after = time.time()
            duration = after - before
            print("\nReading SynTracker input file took " + str(duration) + " seconds.\n")

            # Verify that the df is valid and contain the necessary columns
            if 'Sample1' not in self.score_per_region_all_genomes_df.columns or \
                'Sample2' not in self.score_per_region_all_genomes_df.columns or \
                'Region' not in self.score_per_region_all_genomes_df.columns or \
                'Synteny_score' not in self.score_per_region_all_genomes_df.columns:
                error = "The format of the uploaded SynTracker file is wrong - please verify that you upload " \
                        "SynTracker's output table 'synteny_scores_per_region.csv'."
                self.display_error_page(error)

            # If the input file contains no 'Ref_genome column', add it with artificial genome name
            if 'Ref_genome' not in self.score_per_region_all_genomes_df.columns:
                self.score_per_region_all_genomes_df['Ref_genome'] = 'Reference Genome'

            # Extract the number of genomes from the input file
            self.number_of_genomes_syntracker = self.score_per_region_all_genomes_df.groupby(['Ref_genome']).ngroups
            print("\nNumber of species in SynTracker file: " + str(self.number_of_genomes_syntracker))

            # Get the genomes list, sorted by name
            self.ref_genomes_list_syntracker = list(self.score_per_region_all_genomes_df.groupby(['Ref_genome']).groups)
            #print("\Species list in SynTracker file: " + str(self.ref_genomes_list_syntracker))

            # Calculate the number of pairs at 40 regions for all the genomes and create a sorted list of
            # genomes
            self.ref_genomes_list_by_pairs_num_syntracker = \
                dm.create_sorted_by_pairs_genomes_list_syntracker(self.score_per_region_all_genomes_df)

        # Read the ANI input
        if self.input_mode == "ANI" or self.input_mode == "both":

            # Read the file directly from the path
            before = time.time()
            if self.is_ani_file_path == 1:
                ani_scores_all_genomes_df = pd.read_table(self.ANI_text_input.value)

            # Read the content of the uploaded file
            else:
                ani_scores_all_genomes_df = pd.read_table(io.BytesIO(self.ANI_input_file.value))

            after = time.time()
            duration = after - before
            print("\nReading ANI input file took " + str(duration) + " seconds.\n")

            # Verify that the file contains 4 columns
            if ani_scores_all_genomes_df.shape[1] != 4:
                error = "The format of the uploaded ANI file is wrong - " \
                        "please verify that your file is tab-delimited and contains the following columns: " \
                        "'Ref_genome', 'Sample1', 'Sample2', 'ANI'."
                self.display_error_page(error)

            # To make sure that the column names are unified, replace the original names
            new_column_names = config.ANI_col_names
            ani_scores_all_genomes_df.columns = new_column_names

            # Filter out species with less than 10 pairwise comparisons and return the sorted-by-pairs-num genomes list
            self.ani_scores_all_genomes_df, self.ref_genomes_list_by_pairs_num_ani = \
                dm.filter_genomes_ani(ani_scores_all_genomes_df)

            # Sort the genomes list alphabetically
            self.ref_genomes_list_ani = sorted(self.ref_genomes_list_by_pairs_num_ani)

            # Get the number of genomes
            self.number_of_genomes_ani = len(self.ref_genomes_list_ani)
            print("\nNumber of species in ANI file: " + str(self.number_of_genomes_ani))
            #print("\Species list in ANI file: " + str(self.ref_genomes_list_ani))

        # Set the genomes list according to the input mode
        if self.input_mode == "SynTracker":
            self.ref_genomes_list = self.ref_genomes_list_syntracker
            self.ref_genomes_list_by_pairs_num = self.ref_genomes_list_by_pairs_num_syntracker
            self.number_of_genomes = self.number_of_genomes_syntracker
        elif self.input_mode == "ANI":
            self.ref_genomes_list = self.ref_genomes_list_ani
            self.ref_genomes_list_by_pairs_num = self.ref_genomes_list_by_pairs_num_ani
            self.number_of_genomes = self.number_of_genomes_ani
        else:
            # Unite the two genome lists (they are not necessarily identical)
            self.ref_genomes_list = sorted(set(self.ref_genomes_list_syntracker) | set(self.ref_genomes_list_ani))
            self.ref_genomes_list_by_pairs_num = list(dict.fromkeys(self.ref_genomes_list_by_pairs_num_syntracker +
                                                      self.ref_genomes_list_by_pairs_num_ani))
            self.number_of_genomes = len(self.ref_genomes_list)
            print("\nNumber of total species from both files: " + str(self.number_of_genomes))
            #print("Total species list from both files: " + str(self.ref_genomes_list))

        # Initialize the annotation dict
        for ref in self.ref_genomes_list:
            self.annotation_per_ref_genome_dict[ref] = None

        # If a metadata file was uploaded - read the file into a DF
        if self.is_metadata:
            before = time.time()
            metadata_df = pd.read_table(io.BytesIO(self.metadata_file.value))
            after = time.time()
            duration = after - before
            print("\nReading metadata file took " + str(duration) + " seconds.\n")
            #print("\nMetadata before validation:")
            #print(metadata_df)

            # Check if the provided metadata is valid and match the sample-IDs.
            # If some samples are missing from the metadata - fill them with np.nan values
            before = time.time()
            if self.input_mode == "SynTracker" or self.input_mode == "both":
                scores_df = self.score_per_region_all_genomes_df
            elif self.input_mode == "ANI":
                scores_df = self.ani_scores_all_genomes_df
            self.metadata_dict, self.groups_per_feature_dict, self.metadata_features_list, error = \
                dm.complete_metadata(scores_df, metadata_df)
            after = time.time()
            duration = after - before

            # There is some problem with the metadata file - print error
            if error != "":
                self.display_error_page(error)
                self.valid_metadata = 0

            else:
                print("\nFilling missing metadata took " + str(duration) + " seconds.\n")

        if self.valid_metadata:

            # Initialize the dictionaries that save the dataframes for all the combinations of ref_genomes and sizes
            # and whether they have already been calculated
            if self.input_mode != "ANI":
                for genome in self.ref_genomes_list:
                    self.calculated_APSS_genome_size_dict[genome] = dict()
                    for size in config.sampling_sizes:
                        self.APSS_all_genomes_all_sizes_dict[size] = pd.DataFrame()
                        self.calculated_APSS_all_genomes_size_dict[size] = 0

            # Input file contains only one ref-genome -> present only single genome visualization
            if self.number_of_genomes == 1:
                self.ref_genome = self.ref_genomes_list[0]
                print("\nReference genome = " + self.ref_genome)

                # Create the single-genome visualization layout
                self.create_single_genome_column()
                self.main_single_column.append(self.ref_genome_column)

                self.main_area.clear()
                self.main_area.append(self.main_single_column)

            # Input file contains more than one ref-genome -> display two tabs, for single- and multi-genome views
            else:
                self.ref_genome = self.ref_genomes_list_by_pairs_num[0]
                print("\nReference genome = " + self.ref_genome)
                #pn.state.location.sync(self.genomes_select, {'value': 'ref_genome'})
                self.genomes_select.options = self.ref_genomes_list_by_pairs_num
                self.genomes_select.value = self.ref_genome

                genome_select_row = pn.Row(self.genomes_select, pn.Spacer(width=20), self.genomes_sort_select)
                self.main_single_column.append(genome_select_row)

                # Create the single-genome visualization layout for the selected ref-genome
                self.create_single_genome_column()
                self.main_single_column.append(self.ref_genome_column)

                # Define watchers for the ref_genome_select and the genome_sorting_select widgets
                # This will create the single genome column with the selected reference genome
                self.genomes_select_watcher = self.genomes_select.param.watch(partial(self.select_ref_genome,
                                                                              self.genomes_select), 'value',
                                                                              onlychanged=True)
                self.genomes_sort_select_watcher = self.genomes_sort_select.param.watch(
                    self.changed_genomes_sorting, 'value', onlychanged=True)
                self.genomes_sort_select_multi_watcher = self.genomes_sort_select_multi.param.watch(
                    self.changed_multi_genomes_sorting, 'value', onlychanged=True)

                self.single_multi_genome_tabs.append(('Single species analyses', self.main_single_column))

                multi_message = "Preparing the multiple species visualization - please wait..."
                self.main_multi_column.append(pn.pane.Markdown(multi_message, styles={'font-size': "20px",
                                                                                      'margin': "0"}))
                self.single_multi_genome_tabs.append(('Multiple species analyses', self.main_multi_column))

                self.single_multi_genome_tabs.param.watch(self.changed_single_multi_tabs, 'active')

                self.main_area.clear()

            self.main_area.append(self.single_multi_genome_tabs)

            self.main_area.append(self.submit_new_file_button())

    def changed_single_multi_tabs(self, event):

        # Create the multiple-genome visualization layout when the user selects the multi-genome tab for the first time
        if self.single_multi_genome_tabs.active == 1 and self.visited_multi_genome_tab == 0:
            before = time.time()
            print("\nCalling create_multi_genome_column to create the multiple-species visualization")
            self.create_multi_genome_column()
            after = time.time()
            duration = after - before
            print("\ncreate_multi_genome_column took " + str(duration) + " seconds.\n")

            self.visited_multi_genome_tab = 1

    def changed_genomes_sorting(self, event):
        sorting_method = self.genomes_sort_select.value
        print("\nChanged species sorting method. Sort by: " + sorting_method)

        # Change the order of the species selection and present the first one of the new order
        # Input mode is 'both', Active tab is 'ANI'
        if self.input_mode == "both" and self.active_input_mode == "ANI":
            if sorting_method == "Species name":
                self.genomes_select.options = self.ref_genomes_list_ani
                self.genomes_select.value = self.ref_genomes_list_ani[0]
            else:
                self.genomes_select.options = self.ref_genomes_list_by_pairs_num_ani
                self.genomes_select.value = self.ref_genomes_list_by_pairs_num_ani[0]

        # All other cases
        else:
            if sorting_method == "Species name":
                self.genomes_select.options = self.ref_genomes_list
                self.genomes_select.value = self.ref_genomes_list[0]
            else:
                self.genomes_select.options = self.ref_genomes_list_by_pairs_num
                self.genomes_select.value = self.ref_genomes_list_by_pairs_num[0]

    def changed_multi_genomes_sorting(self, event):
        sorting_method = self.genomes_sort_select_multi.value
        print("\nChanged sorting method for species subset selection. Sort by: " + sorting_method)

        if sorting_method == "Species name":
            self.genomes_subset_select.options = self.ref_genomes_list
        else:
            self.genomes_subset_select.options = self.ref_genomes_list_by_pairs_num

        # Sort the selected values (if any) according to the order in options
        self.genomes_subset_select.value = [item for item in self.genomes_subset_select.options if item in
                                            self.genomes_subset_select.value]

    def select_ref_genome(self, ref_genome, event):
        self.ref_genome = ref_genome.value
        print("\n\nSelected species = " + self.ref_genome)

        self.init_ref_genome()

        self.create_single_genome_column()

    def init_ref_genome(self):

        # Clear layouts
        self.ref_genome_column.clear()
        self.synteny_single_initial_plots_column.clear()
        self.plots_by_size_single_column.clear()
        self.ani_single_plots_column.clear()
        self.combined_single_plots_column.clear()
        self.sample_sizes_slider.value = config.sampling_sizes[0]
        self.selected_contig_column.clear()
        self.APSS_analyses_single_column.clear()
        self.synteny_per_pos_plot_column.clear()
        self.sampling_size = ""
        self.gff_filename = ""

        if self.input_mode == "SynTracker" or self.input_mode == "both":
            # Stop watching the contig-related widgets
            if self.visited_synteny_per_pos_tab and self.finished_initial_synteny_per_pos_plot and \
                    len(self.contigs_list_by_name) > 1:
                # print("\nUnwatching contig_select and sorting_select widgets")
                self.contig_select.param.unwatch(self.contig_select_watcher)
                self.sorting_select.param.unwatch(self.sorting_select_watcher)
                self.avg_plot_chkbox.param.unwatch(self.avg_plot_chkbox_watcher)
                self.avg_plot_color.param.unwatch(self.avg_plot_color_watcher)
                self.coverage_plot_chkbox.param.unwatch(self.coverage_plot_chkbox_watcher)
                self.coverage_plot_color.param.unwatch(self.coverage_plot_color_watcher)
                self.hypervar_chkbox.param.unwatch(self.hypervar_chkbox_watcher)
                self.hypervar_color.param.unwatch(self.hypervar_color_watcher)
                self.hypercons_chkbox.param.unwatch(self.hypercons_chkbox_watcher)
                self.hypercons_color.param.unwatch(self.hypercons_color_watcher)
                self.alpha_slider.param.unwatch(self.alpha_slider_watcher)
                self.synteny_per_pos_feature_select.param.unwatch(self.synteny_per_pos_feature_select_watcher)
                self.show_annotations_chkbox.param.unwatch(self.show_annotations_watcher)

            # Initialize variables
            self.contigs_list_by_length = []
            self.contigs_list_by_name = []
            self.visited_synteny_per_pos_tab = 0
            self.finished_initial_synteny_per_pos_plot = 0
            self.synteny_single_tabs.active = 0
            self.highlight_sample_chkbox.value = False
            self.highlight_nodes_by_feature.value = False
            self.color_edges_by_feature.value = False
            if self.input_mode == "both":
                self.synteny_ani_single_tabs.active = 0
                #self.synteny_ani_multi_tabs.active = 0

        del self.score_per_pos_contig
        del self.score_per_pos_contig_filtered
        del self.df_for_jitter
        del self.df_for_jitter_ani
        del self.scores_matrix
        del self.scores_matrix_ani
        del self.df_for_network
        del self.avg_score_per_pos_contig
        del self.avg_score_per_pos_contig_filtered
        del self.APSS_by_genome_all_sizes_dict

        self.df_for_jitter = pd.DataFrame()
        self.df_for_jitter_ani = pd.DataFrame()
        self.scores_matrix = pd.DataFrame()
        self.scores_matrix_ani = pd.DataFrame()
        self.df_for_network = pd.DataFrame()
        self.score_per_pos_contig = pd.DataFrame()
        self.score_per_pos_contig_filtered = pd.DataFrame()
        self.avg_score_per_pos_contig = pd.DataFrame()
        self.avg_score_per_pos_contig_filtered = pd.DataFrame()
        self.APSS_by_genome_all_sizes_dict = dict()

        self.avg_plot_chkbox.value = True
        self.coverage_plot_chkbox.value = True
        self.hypervar_chkbox.value = True
        self.hypercons_chkbox.value = True
        self.filter_plot_by_metadata = 0

        gc.collect()

    def create_single_genome_column(self):
        ref_genome_title = "Species: " + self.ref_genome

        self.ref_genome_column.append(pn.pane.Markdown(ref_genome_title,
                                                       styles={'font-size': "22px", 'color': config.title_purple_color,
                                                               'margin': "0"}))

        if self.input_mode == "SynTracker":
            self.create_single_genome_column_syntracker_mode()
            self.ref_genome_column.append(self.synteny_single_initial_plots_column)

        elif self.input_mode == "ANI":
            self.create_single_genome_column_ANI_mode()
            self.ref_genome_column.append(self.ani_single_plots_column)

        else:
            self.create_single_genome_column_both_mode()
            self.ref_genome_column.append(self.synteny_ani_single_tabs)

    def create_single_genome_column_both_mode(self):
        before = time.time()
        print("\nStart create_single_genome_column_both_mode")

        self.synteny_single_initial_plots_column.styles = config.both_mode_SynTracker_single_style
        self.ani_single_plots_column.styles = config.both_mode_other_style

        # Run the task of creating the ANI single view in another thread.
        thread_ani = threading.Thread(target=self.create_single_genome_column_ANI_mode)
        thread_ani.start()  # Start the thread

        self.create_single_genome_column_syntracker_mode()
        self.create_single_genome_column_combined_mode()

        self.synteny_ani_single_tabs.clear()
        self.synteny_ani_single_tabs.append(('SynTracker', self.synteny_single_initial_plots_column))
        self.synteny_ani_single_tabs.append(('ANI', self.ani_single_plots_column))
        self.synteny_ani_single_tabs.append(('Combined', self.combined_single_plots_column))

        after = time.time()
        duration = after - before
        print("\ncreate_single_genome_column_both_mode took " + str(duration) + " seconds.\n")

    def create_single_genome_column_syntracker_mode(self):

        # Verify that the current ref-genome is found in the SynTracker input
        if self.ref_genome in self.ref_genomes_list_syntracker:
            before = time.time()

            synteny_per_pos_message = "Preparing the plot - please wait..."
            self.synteny_per_pos_plot_column.append(pn.pane.Markdown(synteny_per_pos_message,
                                                                     styles={'font-size': "20px",
                                                                             'margin': "0"}))

            # Get the score-per-region table for the selected genome only
            self.score_per_region_selected_genome_df = ds.return_selected_genome_table(self.score_per_region_all_genomes_df,
                                                                                       self.ref_genome)

            # Run the task of creating the initial synteny_per_pos plots tab (without the plot itself) in another thread.
            thread_synteny_per_pos = threading.Thread(target=self.create_initial_synteny_per_pos_plot_tab)
            thread_synteny_per_pos.start()  # Start the thread

            # Initialize the dictionary that holds the calculated sampleing sizes
            for size in config.sampling_sizes:
                self.APSS_by_genome_all_sizes_dict[size] = pd.DataFrame()
                self.calculated_APSS_genome_size_dict[size] = 0

            # Create the df for plot presenting the number of pairs vs. subsampled regions
            pairs_num_per_sampling_size_selected_genome_df = \
                ds.create_pairs_num_per_sampling_size(self.score_per_region_selected_genome_df)

            # If the number of pairs with 40 sampled regions is smaller than 100, present also the 'All regions' bar
            # If not (or if it is equal to the total number of pairs), do not present this bar
            # (and remove this option from the slider)
            pairs_at_40 = pairs_num_per_sampling_size_selected_genome_df['Number_of_pairs'].iloc[1]
            self.total_pairs_genome = pairs_num_per_sampling_size_selected_genome_df['Number_of_pairs'].iloc[0]
            if pairs_at_40 >= config.min_pairs_for_all_regions or pairs_at_40 == self.total_pairs_genome:
                self.sample_sizes_slider.options = config.sampling_sizes_wo_all
                self.sample_sizes_slider.value = config.sampling_sizes_wo_all[0]
                is_all_regions = 0
            else:
                is_all_regions = 1
                self.sample_sizes_slider.options = config.sampling_sizes
                self.sample_sizes_slider.value = config.sampling_sizes[0]

            # Create the number of pairs vs. subsampled regions bar plot
            pairs_vs_sampling_size_bar_plot = pn.bind(ps.plot_pairs_vs_sampling_size_bar,
                                                      df=pairs_num_per_sampling_size_selected_genome_df,
                                                      sampling_size=self.sample_sizes_slider,
                                                      is_all_regions=is_all_regions)
            pairs_bar_plot_pane = pn.pane.HoloViews(pairs_vs_sampling_size_bar_plot, width=530, sizing_mode="fixed")

            # Create a markdown for the pairs lost percent (binded to the slider)
            binded_text = pn.bind(widgets.create_pairs_lost_text, pairs_num_per_sampling_size_selected_genome_df,
                                  self.sample_sizes_slider)

            pairs_plot_column = pn.Column(pn.pane.Markdown(refs=binded_text, align='center'), pairs_bar_plot_pane,
                                          styles={'background-color': "white"})

            # Create the number of samples vs. subsampled regions bar plot
            samples_vs_sampling_size_bar_plot = pn.bind(ps.plot_samples_vs_sampling_size_bar,
                                                        df=pairs_num_per_sampling_size_selected_genome_df,
                                                        sampling_size=self.sample_sizes_slider,
                                                        is_all_regions=is_all_regions)
            samples_bar_plot_pane = pn.pane.HoloViews(samples_vs_sampling_size_bar_plot, width=530, sizing_mode="fixed")

            # Create a markdown for the number of samples (binded to the slider)
            binded_text = pn.bind(widgets.create_samples_num_text, pairs_num_per_sampling_size_selected_genome_df,
                                  self.sample_sizes_slider)

            samples_plot_column = pn.Column(pn.pane.Markdown(refs=binded_text, align='center'), samples_bar_plot_pane,
                                            styles={'background-color': "white"})

            plots_row = pn.Row(pairs_plot_column, pn.Spacer(width=20), samples_plot_column)
            slider_row = pn.Row(self.sample_sizes_slider, align='center')
            button_row = pn.Row(self.show_single_plots_button, align='center')

            initial_plots_column = pn.Column(
                plots_row,
                pn.Spacer(height=20),
                slider_row,
                button_row,
                self.plots_by_size_single_column,
            )
            self.APSS_analyses_single_column.append(initial_plots_column)

            self.synteny_single_tabs.clear()
            self.synteny_single_tabs.append(('APSS-based analyses', self.APSS_analyses_single_column))
            self.synteny_single_tabs.append(('Synteny per position', self.synteny_per_pos_plot_column))

            self.synteny_single_initial_plots_column.append(self.synteny_single_tabs)

            self.synteny_single_tabs.param.watch(self.changed_single_tabs, 'active')

            after = time.time()
            duration = after - before
            print("\ncreate_single_genome_column_syntracker_mode took " + str(duration) + " seconds.\n")

        # The ref-genome is not found in the SynTracker file - display a message
        else:
            message = "Species " + self.ref_genome + " is not found in the SynTracker input file."
            self.synteny_single_initial_plots_column.clear()
            self.synteny_single_initial_plots_column.append(pn.pane.Markdown(message,
                                                                             styles={'font-size': "20px",
                                                                                     'color': config.title_red_color,
                                                                                     'padding': "20px",
                                                                                     }))

    def create_single_genome_column_ANI_mode(self):
        self.jitter_card_ani.clear()
        self.metadata_jitter_card_ani.clear()
        self.clustermap_card_ani.clear()
        self.metadata_clustermap_card_ani.clear()
        self.network_card_ani.clear()
        self.metadata_colorby_card_ani.clear()
        self.layout_parameters_card_ani.clear()
        self.network_iterations_ani.value = config.network_iterations_options[0]
        self.highlight_sample_chkbox_ani.value = False
        self.highlight_nodes_by_feature_ani.value = False
        self.color_edges_by_feature_ani.value = False

        before = time.time()
        print("\n\nStart create_single_genome_column_ANI_mode in another thread.")

        # Unwatch ANI plots related watchers (if it's not the first time that this function is called
        # and only for watchers that have been defined before)
        if self.visited_ANI_tab:
            if not isinstance(self.threshold_select_ani_watcher, str):
                self.network_threshold_select_ani.param.unwatch(self.threshold_select_ani_watcher)
            if not isinstance(self.threshold_input_ani_watcher, str):
                self.network_threshold_input_ani.param.unwatch(self.threshold_input_ani_watcher)
            if not isinstance(self.highlight_sample_ani_watcher, str):
                self.highlight_sample_input_ani.param.unwatch(self.highlight_sample_ani_watcher)
            self.network_threshold_input_ani.value = config.ANI_connections_threshold_default
            
            if self.is_metadata:
                if not isinstance(self.feature_colormap_ani_watcher, str):
                    self.feature_colormap_ani.param.unwatch(self.feature_colormap_ani_watcher)
                if not isinstance(self.custom_colormap_clustermap_ani_watcher, str):
                    self.custom_colormap_input_clustermap_ani.param.unwatch(self.custom_colormap_clustermap_ani_watcher)
                if not isinstance(self.color_by_feature_clustermap_ani_watcher, str):
                    self.color_by_feature_ani.param.unwatch(self.color_by_feature_clustermap_ani_watcher)
                if not isinstance(self.continuous_clustermap_ani_watcher, str):
                    self.is_continuous_clustermap_ani.param.unwatch(self.continuous_clustermap_ani_watcher)
                if not isinstance(self.jitter_feature_ani_watcher, str):
                    self.jitter_feature_select_ani.param.unwatch(self.jitter_feature_ani_watcher)
                if not isinstance(self.continuous_network_ani_watcher, str):
                    self.is_continuous_network_ani.param.unwatch(self.continuous_network_ani_watcher)
                if not isinstance(self.colormap_ani_watcher, str):
                    self.nodes_colormap_ani.param.unwatch(self.colormap_ani_watcher)
                if not isinstance(self.custom_colormap_ani_watcher, str):
                    self.custom_colormap_input_ani.param.unwatch(self.custom_colormap_ani_watcher)
                if not isinstance(self.nodes_colorby_ani_watcher, str):
                    self.nodes_color_by_ani.param.unwatch(self.nodes_colorby_ani_watcher)
                if not isinstance(self.nodes_highlight_by_ani_watcher, str):
                    self.nodes_highlight_by_ani.param.unwatch(self.nodes_highlight_by_ani_watcher)

        self.is_continuous_clustermap_ani.value = False
        self.feature_colormap_ani.options = config.categorical_colormap_dict
        self.feature_colormap_ani.value = config.categorical_colormap_dict['cet_glasbey']

        # Verify that the current ref-genome is found in the ANI input
        if self.ref_genome in self.ref_genomes_list_ani:

            # Get the ANI scores table for the selected genome only
            self.ani_scores_selected_genome_df = ds.return_selected_genome_ani_table(self.ani_scores_all_genomes_df,
                                                                                     self.ref_genome)
            #print("\nANI df selected genome:")
            #print(self.ani_scores_selected_genome_df)

            # Display the number of samples and the number of compared pairs
            num_samples = pd.unique(self.ani_scores_selected_genome_df[["Sample1", "Sample2"]].values.ravel()).shape[0]
            num_pairs = len(self.ani_scores_selected_genome_df)
            samples_title = "Number of samples: " + str(num_samples)
            pairs_title = "Number of compared sample pairs: " + str(num_pairs)
            combined_title = samples_title + "\n" + pairs_title
            self.ani_single_plots_column.append(pn.pane.Markdown(combined_title, hard_line_break=True,
                                                                 styles={'font-size': "18px",
                                                                         'margin': "0 0 10px 0",
                                                                         'padding': "0"}))

            # Add the plots to the layout
            self.create_jitter_pane_ani(self.ani_scores_selected_genome_df)
            self.create_clustermap_pane_ani(self.ani_scores_selected_genome_df)
            self.create_network_pane_ani(self.ani_scores_selected_genome_df)
            self.visited_ANI_tab = 1

            plots_column = pn.Column(
                self.jitter_card_ani,
                pn.Spacer(height=20),
                self.clustermap_card_ani,
                pn.Spacer(height=20),
                self.network_card_ani
            )
            self.ani_single_plots_column.append(plots_column)

        # Display a message that the species is not found in the ANI input file
        else:
            message = "Species " + self.ref_genome + " is not found in the ANI input file."
            self.ani_single_plots_column.append(pn.pane.Markdown(message,
                                                                 styles={'font-size': "20px",
                                                                         'color': config.title_red_color,
                                                                         }))

        after = time.time()
        duration = after - before
        print("\ncreate_single_genome_column_ANI_mode took " + str(duration) + " seconds.\n")

    def create_single_genome_column_combined_mode(self):

        before = time.time()

        # Verify that the current ref-genome is found in the ANI input. If not, display an error message
        if self.ref_genome not in self.ref_genomes_list_ani:
            message = "Species " + self.ref_genome + " is not found in the ANI input file - " \
                                                     "cannot display the combined plot."
            self.combined_single_plots_column.append(pn.pane.Markdown(message,
                                                                      styles={'font-size': "20px",
                                                                              'color': config.title_red_color,
                                                                              }))
        # Verify that the current ref-genome is found in the SynTracker input. If not, display an error message
        elif self.ref_genome not in self.ref_genomes_list_syntracker:
            message = "Species " + self.ref_genome + " is not found in the SynTracker input file - " \
                                                     "cannot display the combined plot."
            self.combined_single_plots_column.append(pn.pane.Markdown(message,
                                                                      styles={'font-size': "20px",
                                                                              'color': config.title_red_color,
                                                                              }))

        else:
            # Verify that there is a specific subsampling size has already been selected and calculated
            if self.sampling_size == "":
                message = "No subsampling size for the APSS-based analyses was selected and calculated.\n" \
                          "Please return to the 'SynTracker' tab, select a subsampling size (using the slider) " \
                          "and click the button 'Display plots using the selected number of regions'."

                self.combined_single_plots_column.append(pn.pane.Markdown(message, hard_line_break=True,
                                                                          styles={'font-size': "20px",
                                                                                  'color': config.title_red_color
                                                                                  }))

        after = time.time()
        duration = after - before
        print("\ncreate_single_genome_column_combined_mode took " + str(duration) + " seconds.\n")

    def create_initial_synteny_per_pos_plot_tab(self):

        # Get the sorted contigs lists
        print("\n\nStart create_initial_synteny_per_pos_plot_tab in another thread.")
        before = time.time()
        self.contigs_list_by_name, self.contigs_list_by_length = \
            ds.return_sorted_contigs_lists(self.score_per_region_selected_genome_df)
        after = time.time()
        duration = after - before
        #print("return_sorted_contigs_lists took " + str(duration) + " seconds.\n")
        #print(self.score_per_region_selected_genome_df)

        # Calculate the average score and std for the whole genome (all contigs)
        self.avg_score_genome = self.score_per_region_selected_genome_df['Synteny_score'].mean()
        self.std_score_genome = self.score_per_region_selected_genome_df['Synteny_score'].std()
        median = self.score_per_region_selected_genome_df['Synteny_score'].median()
        print("\nAverage score for the genome = " + str(self.avg_score_genome))
        print("Std of score for the genome = " + str(self.std_score_genome))
        print("Median score for the genome = " + str(median))

        # Calculate the avg score per region for all the regions of all the contigs
        avg_score_per_region = self.score_per_region_selected_genome_df[['Contig_name', 'Position', 'Synteny_score']]. \
            groupby(['Contig_name', 'Position']).\
            agg(Count=('Synteny_score', 'size'), Avg_synteny_score=('Synteny_score', 'mean')).\
            sort_values(['Count'], ascending=False).reset_index()
        #print("\nAvg. score per region df for all the regions of all the contigs:")
        #print(avg_score_per_region)

        self.median_counts = avg_score_per_region['Count'].median()
        self.bottom_percentile_counts = avg_score_per_region['Count'].quantile(0.1)
        print("\nMedian of pairs per region counts is: " + str(self.median_counts))
        print("Percentile 10 of pairs per region counts is: " + str(self.bottom_percentile_counts))

        self.top_percentile = avg_score_per_region['Avg_synteny_score'].quantile(config.top_percentile)
        self.bottom_percentile = avg_score_per_region['Avg_synteny_score'].quantile(config.bottom_percentile)
        print("\nPercentile " + str(config.top_percentile) + " score for the genome = " + str(self.top_percentile))
        print("Percentile " + str(config.bottom_percentile) + " score for the genome = " + str(self.bottom_percentile))

        print("\nTotal number of pairs for the genome = " + str(self.total_pairs_genome))

        contigs_num = len(self.contigs_list_by_name)
        print("\nNumber of contigs: " + str(contigs_num))

        # If there is more than one contig - display a select widget to select the contig
        if contigs_num > 1:

            # Create a drop-down menu of the contigs
            self.contig_select.options = self.contigs_list_by_length
            self.contig_select.value = self.contigs_list_by_length[0]
            # Initialize the sorting method
            self.sorting_select.options = config.contig_sorting_options
            self.sorting_select.value = config.contig_sorting_options[0]

            contig_select_row = pn.Row(self.contig_select, pn.Spacer(width=20), self.sorting_select)

            self.synteny_per_pos_plot_column.clear()
            self.synteny_per_pos_plot_column.append(contig_select_row)
            self.synteny_per_pos_plot_column.append(self.selected_contig_column)

        # There is only one contig
        else:
            self.synteny_per_pos_plot_column.clear()
            self.synteny_per_pos_plot_column.append(self.selected_contig_column)

        self.finished_initial_synteny_per_pos_plot = 1
        print("\nFinished creating initial synteny_per_pos display")

    def changed_single_tabs(self, event):

        # The synteny_per_pos plots tab is selected for the first time for the current reference genome
        if self.synteny_single_tabs.active == 1 and self.visited_synteny_per_pos_tab == 0:

            # In case calculating the initial tab (running in a different thread) hasn't finished yet,
            # wait for it to finish
            while self.finished_initial_synteny_per_pos_plot == 0:
                time.sleep(1)

            # Building the initial synteny_per_pos plots tab has finished
            contigs_num = len(self.contigs_list_by_name)

            # If there is more than one contig -
            # trigger changed_contig to create the synteny_per_pos plot for the selected contig
            if contigs_num > 1:
                # When the selection of contig is changed, a new column is created for the selected contig
                self.contig_select_watcher = self.contig_select.param.watch(partial(self.changed_contig,
                                                                                    self.contig_select),
                                                                            'value', onlychanged=True)

                # When the selection of sorting method is changed, sort the contigs in the contigs-select drop-down
                # menu accordingly
                self.sorting_select_watcher = self.sorting_select.param.watch(partial(self.changed_contig_sorting,
                                                                                      self.sorting_select,
                                                                                      self.contig_select),
                                                                              'value', onlychanged=True)
                self.contig_select.param.trigger('value')

            # Create the synteny_per_pos plot for the Ref-genome
            else:
                self.contig_name = self.contigs_list_by_name[0]
                self.create_selected_contig_column()

            self.visited_synteny_per_pos_tab = 1

    def create_single_genome_plots_by_APSS(self, event):

        self.sampling_size = self.sample_sizes_slider.value
        print("\nSingle species visualization. Selected subsampling size = " + self.sampling_size)

        self.plots_by_size_single_column.clear()
        self.jitter_card.clear()
        self.clustermap_card.clear()
        self.metadata_clustermap_card.clear()
        self.network_card.clear()
        self.layout_parameters_card.clear()
        self.metadata_colorby_card.clear()
        self.metadata_jitter_card.clear()
        self.network_iterations.value = config.network_iterations_options[0]
        self.highlight_sample_chkbox.value = False
        self.highlight_nodes_by_feature.value = False
        self.color_edges_by_feature.value = False

        # Unwatch watchers (if it's not the first time that this function is called)
        if self.clicked_button_display_APSS:
            if not isinstance(self.threshold_select_watcher, str):
                self.network_threshold_select.param.unwatch(self.threshold_select_watcher)
            if not isinstance(self.threshold_input_watcher, str):
                self.network_threshold_input.param.unwatch(self.threshold_input_watcher)
            if not isinstance(self.highlight_sample_watcher, str):
                self.highlight_sample_input.param.unwatch(self.highlight_sample_watcher)
            self.network_threshold_input.value = config.APSS_connections_threshold_default
            if self.is_metadata:
                if not isinstance(self.continuous_network_watcher, str):
                    self.is_continuous_network.param.unwatch(self.continuous_network_watcher)
                if not isinstance(self.colormap_watcher, str):
                    self.nodes_colormap.param.unwatch(self.colormap_watcher)
                if not isinstance(self.custom_colormap_watcher, str):
                    self.custom_colormap_input.param.unwatch(self.custom_colormap_watcher)
                if not isinstance(self.nodes_colorby_watcher, str):
                    self.nodes_color_by.param.unwatch(self.nodes_colorby_watcher)
                if not isinstance(self.nodes_highlight_by_watcher, str):
                    self.nodes_highlight_by.param.unwatch(self.nodes_highlight_by_watcher)
                if not isinstance(self.feature_colormap_watcher, str):
                    self.feature_colormap.param.unwatch(self.feature_colormap_watcher)
                if not isinstance(self.custom_colormap_clustermap_watcher, str):
                    self.custom_colormap_input_clustermap.param.unwatch(self.custom_colormap_clustermap_watcher)
                if not isinstance(self.color_by_feature_clustermap_watcher, str):
                    self.color_by_feature.param.unwatch(self.color_by_feature_clustermap_watcher)
                if not isinstance(self.continuous_clustermap_watcher, str):
                    self.is_continuous_clustermap.param.unwatch(self.continuous_clustermap_watcher)
                if not isinstance(self.jitter_feature_watcher, str):
                    self.jitter_feature_select.param.unwatch(self.jitter_feature_watcher)
        self.clicked_button_display_APSS = 1

        self.is_continuous_network.value = False
        self.is_continuous_clustermap.value = False
        self.nodes_colormap.options = config.categorical_colormap_dict
        self.nodes_colormap.value = config.categorical_colormap_dict['cet_glasbey']
        self.feature_colormap.options = config.categorical_colormap_dict
        self.feature_colormap.value = config.categorical_colormap_dict['cet_glasbey']

        # Check if the requested genome and size have already been calculated. If so, fetch the specific dataframe
        if self.calculated_APSS_genome_size_dict[self.sampling_size]:
            print("\nThe selected size (" + self.sampling_size + ") has already been calculated - retrieve it.")
            selected_genome_and_size_avg_df = \
                self.APSS_by_genome_all_sizes_dict[self.sampling_size]

        else:
            # Calculate and return the dataframe with average scores for the selected genome and sampling size
            print("\nThe selected size (" + self.sampling_size + ") has not been calculated yet - calculate it.")
            selected_genome_and_size_avg_df = ds.calculate_avg_scores_selected_genome_size(
                self.score_per_region_selected_genome_df, self.ref_genome, self.sampling_size)
            # Save the dataframe in the main dictionary
            self.APSS_by_genome_all_sizes_dict[self.sampling_size] = selected_genome_and_size_avg_df
            self.calculated_APSS_genome_size_dict[self.sampling_size] = 1

        # No data at the selected sampling size
        if selected_genome_and_size_avg_df.empty:
            text = "The data obtained using " + self.sampling_size + " subsampled regions is not sufficient for " \
                                                                     "further processing"
            self.plots_by_size_single_column.append(pn.pane.Markdown(text, styles={'font-size': "18px",
                                                                                   'color': config.title_red_color,
                                                                                   'margin': "0"}))

        # Enough data -> creating plots
        else:
            if self.sampling_size == 'All':
                size_title = "Presenting plots using APSS from all available regions:"
            else:
                size_title = "Presenting plots using APSS from " + self.sampling_size + " subsampled regions:"

            self.plots_by_size_single_column.append(pn.pane.Markdown(size_title, styles={'font-size': "20px",
                                                                                         'color': config.title_purple_color,
                                                                                         'margin': "0 0 10px 0",
                                                                                         'padding': "0"}))

            # Add the plots to the layout
            self.create_jitter_pane(selected_genome_and_size_avg_df)
            self.create_clustermap_pane(selected_genome_and_size_avg_df)
            self.create_network_pane(selected_genome_and_size_avg_df)

            # In case of combined mode, verify that the ref-genome is found in both input files
            # and then prepare the combined plot for the selected sampling size
            if self.input_mode == 'both' and self.ref_genome in self.ref_genomes_list_ani:

                APSS_ANI_selected_genome_df = \
                    selected_genome_and_size_avg_df.loc[:, ['Sample1', 'Sample2', 'APSS']].copy()

                # Merge the ANI df and the APSS df, based on 'Sample1' and 'Sample2' - keep only common pairs
                APSS_ANI_selected_genome_df = APSS_ANI_selected_genome_df.merge(
                    self.ani_scores_selected_genome_df[['Sample1', 'Sample2', 'ANI']],
                    on=['Sample1', 'Sample2'],
                    how='inner'  # Only keep rows where the sample pair exists in both DataFrames
                )
                #print("\nCombined APSS and ANI df:")
                #print(APSS_ANI_selected_genome_df)
                self.create_combined_scatter_pane(APSS_ANI_selected_genome_df)

            plots_column = pn.Column(self.jitter_card, pn.Spacer(height=20), self.clustermap_card, pn.Spacer(height=20),
                                     self.network_card)
            self.plots_by_size_single_column.append(plots_column)

    def create_jitter_pane(self, selected_genome_and_size_avg_df):
        styling_title = "Plot styling options:"
        metadata_colors_row = pn.Row(self.jitter_same_color, pn.Spacer(width=3), self.jitter_different_color)
        metadata_col = pn.Column(self.jitter_feature_select,
                                 metadata_colors_row,
                                 styles={'padding': "10x"})
        self.metadata_jitter_card.append(metadata_col)
        options_row = pn.Row(self.jitter_type_select, self.jitter_color)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                options_row,
                                pn.Spacer(height=5),
                                self.use_metadata_jitter,
                                self.metadata_jitter_card
        )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_jitter)

        jitter_file = "Dist_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_jitter_file_path.placeholder = jitter_file

        self.download_jitter_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                 styles={'font-size': "15px",
                                                                         'font-weight': "bold",
                                                                         'color': config.title_blue_color,
                                                                         'margin': "0"}),
                                                self.jitter_image_format, self.save_jitter_file_path,
                                                download_button, pn.pane.Markdown())

        # Use metadata in plot
        if self.is_metadata:
            # Update the color nodes by- drop-down menu with the available metadata features
            self.jitter_feature_select.options = self.metadata_features_list
            self.jitter_feature_select.value = self.metadata_features_list[0]

            # Define a watcher for the jitter_feature_select widget
            self.jitter_feature_watcher = self.jitter_feature_select.param.watch(partial(self.change_jitter_feature,
                                                                                         selected_genome_and_size_avg_df),
                                                                                 'value', onlychanged=True)

        # No metadata
        else:
            self.use_metadata_jitter.disabled = True

        # Create the jitter plot
        self.jitter_plot = pn.bind(self.create_jitter_plot_syntracker, avg_df=selected_genome_and_size_avg_df,
                                   type=self.jitter_type_select, color=self.jitter_color,
                                   use_metadata=self.use_metadata_jitter, feature=self.jitter_feature_select.value,
                                   same_color=self.jitter_same_color, different_color=self.jitter_different_color)

        self.jitter_pane = pn.pane.Matplotlib(self.jitter_plot, width=520, dpi=300, tight=True, format='png')

        jitter_table_file = "Data_for_dist_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_jitter_table_path.placeholder = jitter_table_file
        download_table_button = pn.widgets.Button(name='Download data table in tsv format', button_type='primary')
        download_table_button.on_click(self.download_jitter_table)

        self.download_jitter_table_column = pn.Column(self.save_jitter_table_path,
                                                      download_table_button,
                                                      pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=25), self.download_jitter_column,
                                 self.download_jitter_table_column)

        jitter_row = pn.Row(controls_col, pn.Spacer(width=120), self.jitter_pane, styles={'padding': "15px"})
        self.jitter_card.append(jitter_row)

    def change_jitter_feature(self, selected_genome_and_size_avg_df, event):
        self.update_jitter_plot(selected_genome_and_size_avg_df)

    def update_jitter_plot(self, selected_genome_and_size_avg_df):
        # Update the jitter plot
        self.jitter_plot = pn.bind(self.create_jitter_plot_syntracker, avg_df=selected_genome_and_size_avg_df,
                                   type=self.jitter_type_select, color=self.jitter_color,
                                   use_metadata=self.use_metadata_jitter, feature=self.jitter_feature_select.value,
                                   same_color=self.jitter_same_color, different_color=self.jitter_different_color)
        self.jitter_pane.object = self.jitter_plot

    def create_jitter_pane_ani(self, ani_df):
        styling_title = "Plot styling options:"
        metadata_colors_row = pn.Row(self.jitter_same_color_ani, pn.Spacer(width=3), self.jitter_different_color_ani)
        metadata_col = pn.Column(self.jitter_feature_select_ani,
                                 metadata_colors_row,
                                 styles={'padding': "10x"})
        self.metadata_jitter_card_ani.append(metadata_col)
        options_row = pn.Row(self.jitter_type_select_ani, self.jitter_color_ani)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                options_row,
                                pn.Spacer(height=5),
                                self.use_metadata_jitter_ani,
                                self.metadata_jitter_card_ani
        )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_jitter_ani)

        jitter_file = "ANI_Dist_plot_" + self.ref_genome
        self.save_jitter_file_path_ani.placeholder = jitter_file

        self.download_jitter_column_ani = pn.Column(pn.pane.Markdown(save_file_title,
                                                                     styles={'font-size': "15px",
                                                                             'font-weight': "bold",
                                                                             'color': config.title_blue_color,
                                                                             'margin': "0"}),
                                                    self.jitter_image_format_ani, self.save_jitter_file_path_ani,
                                                    download_button, pn.pane.Markdown())

        # Use metadata in plot
        if self.is_metadata:
            # Update the color nodes by- drop-down menu with the available metadata features
            self.jitter_feature_select_ani.options = self.metadata_features_list
            self.jitter_feature_select_ani.value = self.metadata_features_list[0]

            # Define a watcher for the jitter_feature_select widget
            self.jitter_feature_ani_watcher = self.jitter_feature_select_ani.param.watch(
                partial(self.change_jitter_feature_ani, ani_df), 'value', onlychanged=True)

        # No metadata
        else:
            self.use_metadata_jitter_ani.disabled = True

        # Create the jitter plot
        self.jitter_plot_ani = pn.bind(self.create_jitter_plot_ani, ani_df=ani_df,
                                       type=self.jitter_type_select_ani, color=self.jitter_color_ani,
                                       use_metadata=self.use_metadata_jitter_ani,
                                       feature=self.jitter_feature_select_ani.value,
                                       same_color=self.jitter_same_color_ani,
                                       different_color=self.jitter_different_color_ani)

        self.jitter_pane_ani = pn.pane.Matplotlib(self.jitter_plot_ani, width=520, dpi=300, tight=True, format='png')

        jitter_table_file = "Data_for_ANI_dist_plot_" + self.ref_genome
        self.save_jitter_table_path_ani.placeholder = jitter_table_file
        download_table_button = pn.widgets.Button(name='Download data table in tsv format', button_type='primary')
        download_table_button.on_click(self.download_jitter_table_ani)

        self.download_jitter_table_column_ani = pn.Column(self.save_jitter_table_path_ani,
                                                          download_table_button, pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=25), self.download_jitter_column_ani,
                                 self.download_jitter_table_column_ani)

        jitter_row = pn.Row(controls_col, pn.Spacer(width=120), self.jitter_pane_ani, styles={'padding': "15px"})
        self.jitter_card_ani.append(jitter_row)

    def category_by_feature(self, row, feature):
        if is_missing_value(self.metadata_dict[feature][row['Sample1']]) \
                or is_missing_value(self.metadata_dict[feature][row['Sample2']]):
            return 'Unknown'
        else:
            if self.metadata_dict[feature][row['Sample1']] == self.metadata_dict[feature][row['Sample2']]:
                return 'Same ' + feature
            else:
                return 'Different ' + feature

    def change_jitter_feature_ani(self, ani_df, event):
        self.update_jitter_plot_ani(ani_df)

    def update_jitter_plot_ani(self, ani_df):
        # Update the jitter plot
        self.jitter_plot_ani = pn.bind(self.create_jitter_plot_ani, ani_df=ani_df,
                                       type=self.jitter_type_select_ani, color=self.jitter_color_ani,
                                       use_metadata=self.use_metadata_jitter_ani,
                                       feature=self.jitter_feature_select_ani.value,
                                       same_color=self.jitter_same_color_ani,
                                       different_color=self.jitter_different_color_ani)
        self.jitter_pane_ani.object = self.jitter_plot_ani

    def create_jitter_plot_syntracker(self, avg_df, type, color, use_metadata, feature, same_color, different_color):

        self.df_for_jitter = avg_df.loc[:, ['Sample1', 'Sample2', 'APSS']].copy()

        # Use metadata to separate plot to same/different feature
        if use_metadata:
            self.df_for_jitter['Category'] = self.df_for_jitter.apply(lambda row: self.category_by_feature(row, feature),
                                                                      axis=1)
            same_feature = 'Same ' + feature
            diff_feature = 'Different ' + feature

            # The selected feature contains 'Unknown' values - present them as a third category
            if (self.df_for_jitter['Category'] == 'Unknown').any():
                if type == 'Boxplot':
                    plot = sns.catplot(data=self.df_for_jitter, kind='box', x="Category", y="APSS",
                                       order=[same_feature, diff_feature, 'Unknown'],
                                       hue="Category", hue_order=[same_feature, diff_feature, 'Unknown'],
                                       palette=[same_color, different_color, 'gray'], width=0.5)

                else:
                    plot = sns.catplot(data=self.df_for_jitter, x="Category", y="APSS",
                                       order=[same_feature, diff_feature, 'Unknown'],
                                       hue="Category", hue_order=[same_feature, diff_feature, 'Unknown'],
                                       palette=[same_color, different_color, 'gray'], linewidth=0.1)
            # No Unknowns - present only two categories
            else:
                if type == 'Boxplot':
                    plot = sns.catplot(data=self.df_for_jitter, kind='box', x="Category", y="APSS",
                                       order=[same_feature, diff_feature],
                                       hue="Category", hue_order=[same_feature, diff_feature],
                                       palette=[same_color, different_color], width=0.5)

                else:
                    plot = sns.catplot(data=self.df_for_jitter, x="Category", y="APSS",
                                       order=[same_feature, diff_feature],
                                       hue="Category", hue_order=[same_feature, diff_feature],
                                       palette=[same_color, different_color], linewidth=0.1)

            # Calculate the P-value of the comparison
            same_array = self.df_for_jitter[self.df_for_jitter['Category'] == same_feature]['APSS']
            diff_array = self.df_for_jitter[self.df_for_jitter['Category'] == diff_feature]['APSS']

            # Sample size is enough for P-value calculation
            if len(same_array) >= 1 and len(diff_array) >= 1:
                before = time.time()
                u, p_val = return_p_value_mannwhitneyu(same_array, diff_array)
                effect_size = abs(1 - (2 * u) / (len(same_array) * len(diff_array)))
                after = time.time()
                duration = after - before
                print("\nP-value for " + feature + " comparison = " + str(p_val) + ", Effect size = " + str(effect_size))
                #print("P-value calculation took " + str(duration) + " seconds.\n")

                # Pvalue is valid and significant
                if str(p_val) != "nan" and p_val <= 0.05:
                    ax = plot.ax  # get the underlying matplotlib axis

                    # Place the p-value text between the two boxes, slightly above the max APSS
                    y_max = self.df_for_jitter["APSS"].max()
                    ax.text(
                        0.5, y_max * 1.02,  # x = middle of the two boxes, y = above max
                        f"p <= {p_val: .1e}\nEffect size = {effect_size: .2}",  # format p-value in scientific notation
                        ha="center", va="bottom", fontsize=9
                    )
            # Sample size is not enough for P-value calculation
            else:
                print("\nCannot calculate P-value for feature " + feature + ": Sample size is not enough.")

            # Rotate the x-axis labels in 45 degrees
            ax = plot.ax
            plt.setp(ax.get_xticklabels(), rotation=45)

        # Do not use metadata in plot - show all the comparisons together
        else:
            self.df_for_jitter['Category'] = 'All Comparisons'
            if type == 'Boxplot':
                plot = sns.catplot(data=self.df_for_jitter, kind='box', x="Category", y="APSS", color=color, width=0.5)
            else:
                plot = sns.catplot(data=self.df_for_jitter, x="Category", y="APSS", color=color, linewidth=0.1)

        # Remove the x-axis label
        plot.set_axis_labels("", "APSS")  # Sets x-label to empty string

        #print("\nDF for jitter plot:")
        #print(self.df_for_jitter)

        plt.close(plot.figure)

        return plot.figure

    def create_jitter_plot_ani(self, ani_df, type, color, use_metadata, feature, same_color, different_color):

        self.df_for_jitter_ani = ani_df.loc[:, ['Sample1', 'Sample2', 'ANI']].copy()

        # Use metadata to separate plot to same/different feature
        if use_metadata:
            self.df_for_jitter_ani['Category'] = self.df_for_jitter_ani.apply(
                lambda row: self.category_by_feature(row, feature),
                axis=1)
            print("\ndf_for_jitter_ani:")
            print(self.df_for_jitter_ani)

            same_feature = 'Same ' + feature
            diff_feature = 'Different ' + feature

            # The selected feature contains 'Unknown' values - present them as a third category
            if (self.df_for_jitter_ani['Category'] == 'Unknown').any():
                # Boxplot
                if type == 'Boxplot':
                    plot = sns.catplot(data=self.df_for_jitter_ani, kind='box', x="Category", y="ANI",
                                       order=[same_feature, diff_feature, 'Unknown'],
                                       hue="Category", hue_order=[same_feature, diff_feature, 'Unknown'],
                                       palette=[same_color, different_color, 'gray'], width=0.5)

                # Jitterplot
                else:
                    plot = sns.catplot(data=self.df_for_jitter_ani, x="Category", y="ANI",
                                       order=[same_feature, diff_feature, 'Unknown'],
                                       hue="Category", hue_order=[same_feature, diff_feature, 'Unknown'],
                                       palette=[same_color, different_color, 'gray'], linewidth=0.1)

            # No Unknowns - present only two categories
            else:
                # Boxplot
                if type == 'Boxplot':
                    plot = sns.catplot(data=self.df_for_jitter_ani, kind='box', x="Category", y="ANI",
                                       order=[same_feature, diff_feature],
                                       hue="Category", hue_order=[same_feature, diff_feature],
                                       palette=[same_color, different_color], width=0.5)

                # Jitterplot
                else:
                    plot = sns.catplot(data=self.df_for_jitter_ani, x="Category", y="ANI",
                                       order=[same_feature, diff_feature],
                                       hue="Category", hue_order=[same_feature, diff_feature],
                                       palette=[same_color, different_color], linewidth=0.1)

            # Calculate the P-value of the comparison
            same_array = self.df_for_jitter_ani[self.df_for_jitter_ani['Category'] == same_feature]['ANI']
            diff_array = self.df_for_jitter_ani[self.df_for_jitter_ani['Category'] == diff_feature]['ANI']

            # Sample size is enough for P-value calculation
            if len(same_array) >= 1 and len(diff_array) >= 1:
                before = time.time()
                u, p_val = return_p_value_mannwhitneyu(same_array, diff_array)
                effect_size = abs(1 - (2 * u) / (len(same_array) * len(diff_array)))
                after = time.time()
                duration = after - before
                print("\nP-value for " + feature + " comparison = " + str(p_val) + ", Effect size = " + str(effect_size))
                #print("P-value calculation took " + str(duration) + " seconds.\n")

                # P-value is valid and significant
                if str(p_val) != "nan" and p_val <= 0.05:
                    ax = plot.ax  # get the underlying matplotlib axis

                    # Place the p-value text between the two boxes, slightly above the max APSS
                    y_max = self.df_for_jitter_ani["ANI"].max()
                    ax.text(
                        0.5, y_max * 1.002,  # x = middle of the two boxes, y = above max
                        f"p <= {p_val: .1e}\nEffect size = {effect_size: .2}",  # format p-value in scientific notation
                        ha="center", va="bottom", fontsize=9
                    )
            # Sample size is not enough for P-value calculation
            else:
                print("\nCannot calculate P-value for feature " + feature + ": Sample size is not enough.")

            # Rotate the x-axis labels in case of more than one category
            ax = plot.ax
            plt.setp(ax.get_xticklabels(), rotation=45)

        # Do not use metadata in plot - show all the comparisons together
        else:
            self.df_for_jitter_ani['Category'] = 'All Comparisons'
            if type == 'Boxplot':
                plot = sns.catplot(data=self.df_for_jitter_ani, kind='box', x="Category", y="ANI", color=color,
                                   width=0.5)
            else:
                plot = sns.catplot(data=self.df_for_jitter_ani, x="Category", y="ANI", color=color, linewidth=0.1)

        # Remove the x-axis label
        plot.set_axis_labels("", "ANI")  # Sets x-label to empty string

        # print("\nDF for jitter plot:")
        # print(self.df_for_jitter_ani)

        plt.close(plot.figure)

        return plot.figure

    def download_jitter(self, event):
        fformat = self.jitter_image_format.value
        plot_type = self.jitter_type_select.value
        if plot_type != 'Boxplot':
            plot_type = 'Jitter_plot'

        # Update the placeholder of the filename for download.
        jitter_file = plot_type + "_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        if self.use_metadata_jitter.value:
            jitter_file += "_separate_by_" + self.jitter_feature_select.value
        self.save_jitter_file_path.placeholder = jitter_file

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

        download_message = "The image is successfully saved under:\n" + jitter_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_jitter_column.pop(4)
        self.download_jitter_column.append(download_floatpanel)

    def download_jitter_ani(self, event):
        fformat = self.jitter_image_format_ani.value
        plot_type = self.jitter_type_select_ani.value
        if plot_type != 'Boxplot':
            plot_type = 'Jitter_plot'

        # Update the placeholder of the filename for download.
        jitter_file = "ANI_" + plot_type + "_" + self.ref_genome
        if self.use_metadata_jitter_ani.value:
            jitter_file += "_separate_by_" + self.jitter_feature_select_ani.value
        self.save_jitter_file_path_ani.placeholder = jitter_file

        # Set the directory for saving
        if self.save_jitter_file_path_ani.value == "":
            jitter_file_path = self.downloads_dir_path + self.save_jitter_file_path_ani.placeholder + "." + fformat
        else:
            jitter_file_path = self.save_jitter_file_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, jitter_file_path, re.IGNORECASE):
                jitter_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(jitter_file_path):
                jitter_file_path = self.downloads_dir_path + jitter_file_path

        self.jitter_plot_ani().savefig(jitter_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + jitter_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_jitter_column_ani.pop(4)
        self.download_jitter_column_ani.append(download_floatpanel)

    def download_jitter_table(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_jitter_table_path.value == "":
            jitter_table_path = self.downloads_dir_path + self.save_jitter_table_path.placeholder + "." + fformat
        else:
            jitter_table_path = self.save_jitter_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, jitter_table_path, re.IGNORECASE):
                jitter_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(jitter_table_path):
                jitter_table_path = self.downloads_dir_path + jitter_table_path

        self.df_for_jitter.to_csv(jitter_table_path, sep='\t', index=False)

        download_message = "The table is successfully saved under:\n" + jitter_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_jitter_table_column.pop(2)
        self.download_jitter_table_column.append(download_floatpanel)

    def download_jitter_table_ani(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_jitter_table_path_ani.value == "":
            jitter_table_path = self.downloads_dir_path + self.save_jitter_table_path_ani.placeholder + "." + fformat
        else:
            jitter_table_path = self.save_jitter_table_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, jitter_table_path, re.IGNORECASE):
                jitter_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(jitter_table_path):
                jitter_table_path = self.downloads_dir_path + jitter_table_path

        self.df_for_jitter_ani.to_csv(jitter_table_path, sep='\t', index=False)

        download_message = "The table is successfully saved under:\n" + jitter_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_jitter_table_column_ani.pop(2)
        self.download_jitter_table_column_ani.append(download_floatpanel)

    def create_clustermap_pane(self, selected_genome_and_size_avg_df):

        self.clustermap_method.value = config.clustering_methods[0]

        styling_title = "Heatmap customization options:"
        continuous_col = pn.Column(pn.Spacer(height=20), self.is_continuous_clustermap)
        color_by_feature_row = pn.Row(self.color_by_feature, continuous_col)
        metadata_coloring_col = pn.Column(color_by_feature_row,
                                          self.feature_colormap,
                                          self.custom_colormap_input_clustermap,
                                          styles={'padding': "10x"})
        self.metadata_clustermap_card.append(metadata_coloring_col)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.clustermap_cmap,
                                self.clustermap_method,
                                pn.Spacer(height=5),
                                self.use_metadata_clustermap,
                                self.metadata_clustermap_card)

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_clustermap)

        clustermap_file = "Clustermap_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_clustermap_file_path.placeholder = clustermap_file

        self.download_clustermap_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                     styles={'font-size': "15px",
                                                                             'font-weight': "bold",
                                                                             'color': config.title_blue_color,
                                                                             'margin': "0"}),
                                                    self.clustermap_image_format, self.save_clustermap_file_path,
                                                    download_button, pn.pane.Markdown())
        matrix_file = "Matrix_for_clustermap_" + self.ref_genome + "_" + self.sampling_size + "_regions"
        self.save_matrix_file_path.placeholder = matrix_file
        download_table_button = pn.widgets.Button(name='Download matrix', button_type='primary')
        download_table_button.on_click(self.download_matrix)

        self.download_matrix_column = pn.Column(self.save_matrix_file_path,
                                                download_table_button,
                                                pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_clustermap_column,
                                 self.download_matrix_column)

        # Transform the data into a scoring matrix
        pivoted_df = selected_genome_and_size_avg_df.pivot(columns='Sample1', index='Sample2', values='APSS')
        self.scores_matrix = pivoted_df.combine_first(pivoted_df.T)
        np.fill_diagonal(self.scores_matrix.values, 1.0)
        #print("\nScores matrix:")
        #print(self.scores_matrix)

        # Check the number of columns in the matrix
        col_num = len(self.scores_matrix.columns)
        #print("\ncreate_clustermap_pane: number of columns = " + str(col_num))

        # If the num of columns exceeds the defined maximum, do not create the clustermap plot
        # and display a message + a possibility to download the matrix
        if col_num > config.max_clustermap_cols:
            error = "The number of samples exceeds the limit of " + str(config.max_clustermap_cols) + \
                    " so the heatmap plot cannot be well presented."
            suggestion = "It ia possible to download the scoring matrix and display the heatmap in another program."
            clustermap_col = pn.Column(
                                       pn.pane.Markdown(error, styles={'font-size': "16px",
                                                                       'color': config.title_red_color,
                                                                       'margin-bottom': "0"}),
                                       pn.pane.Markdown(suggestion,
                                                        styles={'font-size': "14px", 'margin-top': "0"}),
                                       self.download_matrix_column,
                                       styles={'padding': "15px"})

            self.clustermap_card.append(clustermap_col)

        # Display the full clustermap pane, including the plot
        else:
            # There is metadata
            if self.is_metadata:
                # Update the color rows by- drop-down menu with the available metadata features
                self.color_by_feature.options = self.metadata_features_list
                self.color_by_feature.value = self.metadata_features_list[0]

                # Define watchers for the metadata widgets
                self.color_by_feature_clustermap_watcher = self.color_by_feature.param.watch(
                    self.set_not_continuous_clustermap, 'value', onlychanged=True)
                self.continuous_clustermap_watcher = self.is_continuous_clustermap.param.watch(
                    self.change_continuous_state_clustermap, 'value', onlychanged=True)
                self.feature_colormap_watcher = self.feature_colormap.param.watch(self.change_colormap_clustermap,
                                                                                  'value', onlychanged=True)
                self.custom_colormap_clustermap_watcher = \
                    self.custom_colormap_input_clustermap.param.watch(self.get_custom_colormap_clustermap,
                                                                      'value', onlychanged=True)

            # No metadata
            else:
                self.use_metadata_clustermap.disabled = True

            self.clustermap_plot = pn.bind(ps.create_clustermap, matrix=self.scores_matrix, type="APSS",
                                           cmap=self.clustermap_cmap, method=self.clustermap_method,
                                           is_metadata=self.use_metadata_clustermap,
                                           feature=self.color_by_feature.value,
                                           is_continuous=self.is_continuous_clustermap.value,
                                           cmap_metadata=self.feature_colormap.value_name,
                                           custom_cmap=self.custom_colormap_input_clustermap.value,
                                           metadata_dict=self.metadata_dict)

            self.clustermap_pane = pn.pane.Matplotlib(self.clustermap_plot, height=600, dpi=300, tight=True,
                                                      format='png')
            clustermap_row = pn.Row(controls_col, pn.Spacer(width=30), self.clustermap_pane, styles={'padding': "15px"})

            self.clustermap_card.append(clustermap_row)

    def change_colormap_clustermap(self, event):
        self.update_clustermap_plot()

    def get_custom_colormap_clustermap(self, event):
        self.update_clustermap_plot()

    def change_continuous_state_clustermap(self, event):
        # Continuous feature
        if self.is_continuous_clustermap.value:
            #print("\nIn change_continuous_state. Continuous feature")

            # Verify that the feature is indeed continuous
            feature = self.color_by_feature.value
            unique_groups = sorted(
                list(set([self.metadata_dict[feature][sample] for sample in self.scores_matrix.iloc[:, 0].index])))
            str_type = 0
            for group in unique_groups:
                if isinstance(group, str):
                    str_type = 1

            # Feature is not really continuous, treat as categorical
            if str_type == 1:
                #print("The feature is not really continuous - uncheck...")
                self.is_continuous_clustermap.value = False

            # Feature is indeed really continuous
            else:
                self.feature_colormap.options = config.continuous_colormap_dict
                self.feature_colormap.value = config.continuous_colormap_dict['plasma_r']

        # Categorical feature
        else:
            #print("\nIn change_continuous_state. Categorical feature")
            self.feature_colormap.options = config.categorical_colormap_dict
            self.feature_colormap.value = config.categorical_colormap_dict['cet_glasbey']

    def set_not_continuous_clustermap(self, event):
        #print("\nIn set_not_continuous")
        self.is_continuous_clustermap.value = False
        self.update_clustermap_plot()

    # Update the clustermap plot
    def update_clustermap_plot(self):
        #print("\nIn update_clustermap_plot")
        self.clustermap_plot = pn.bind(ps.create_clustermap, matrix=self.scores_matrix, type="APSS",
                                       cmap=self.clustermap_cmap, method=self.clustermap_method,
                                       is_metadata=self.use_metadata_clustermap,
                                       feature=self.color_by_feature.value,
                                       is_continuous=self.is_continuous_clustermap.value,
                                       cmap_metadata=self.feature_colormap.value_name,
                                       custom_cmap=self.custom_colormap_input_clustermap.value,
                                       metadata_dict=self.metadata_dict)

        self.clustermap_pane.object = self.clustermap_plot

    def download_clustermap(self, event):
        fformat = self.clustermap_image_format.value

        # Update the placeholder of the filename for download.
        clustermap_file = "Clustermap_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                          self.clustermap_cmap.value + "_" + self.clustermap_method.value
        if self.use_metadata_clustermap.value:
            if self.feature_colormap.value_name == 'Define custom colormap':
                clustermap_file += "_colorby_" + self.color_by_feature.value + "_custom_colormap"
            else:
                clustermap_file += "_colorby_" + self.color_by_feature.value + "_" + self.feature_colormap.value_name
        self.save_clustermap_file_path.placeholder = clustermap_file

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

        download_message = "The image is successfully saved under:\n" + clustermap_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_clustermap_column.pop(4)
        self.download_clustermap_column.append(download_floatpanel)

    def download_matrix(self, event):
        fformat = "txt"

        # Set the directory for saving
        if self.save_matrix_file_path.value == "":
            matrix_path = self.downloads_dir_path + self.save_matrix_file_path.placeholder + "." + fformat
        else:
            matrix_path = self.save_matrix_file_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, matrix_path, re.IGNORECASE):
                matrix_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(matrix_path):
                matrix_path = self.downloads_dir_path + matrix_path

        self.scores_matrix.to_csv(matrix_path, sep='\t')

        download_message = "The table is successfully saved under:\n" + matrix_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_matrix_column.pop(2)
        self.download_matrix_column.append(download_floatpanel)

    def create_clustermap_pane_ani(self, selected_genome_ani_df):

        self.clustermap_method_ani.value = config.clustering_methods[0]

        styling_title = "Heatmap customization options:"
        continuous_col = pn.Column(pn.Spacer(height=20), self.is_continuous_clustermap_ani)
        color_by_feature_row = pn.Row(self.color_by_feature_ani, continuous_col)
        metadata_coloring_col = pn.Column(color_by_feature_row,
                                          self.feature_colormap_ani,
                                          self.custom_colormap_input_clustermap_ani,
                                          styles={'padding': "10x"})
        self.metadata_clustermap_card_ani.append(metadata_coloring_col)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.clustermap_cmap_ani,
                                self.clustermap_method_ani,
                                pn.Spacer(height=5),
                                self.use_metadata_clustermap_ani,
                                self.metadata_clustermap_card_ani)

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_clustermap_ani)

        clustermap_file = "ANI_Clustermap_" + self.ref_genome
        self.save_clustermap_file_path_ani.placeholder = clustermap_file

        self.download_clustermap_column_ani = pn.Column(pn.pane.Markdown(save_file_title,
                                                                         styles={'font-size': "15px",
                                                                                 'font-weight': "bold",
                                                                                 'color': config.title_blue_color,
                                                                                 'margin': "0"}),
                                                        self.clustermap_image_format_ani,
                                                        self.save_clustermap_file_path_ani,
                                                        download_button, pn.pane.Markdown())
        matrix_file = "Matrix_for_ANI_Clustermap_" + self.ref_genome
        self.save_matrix_file_path_ani.placeholder = matrix_file
        download_table_button = pn.widgets.Button(name='Download matrix', button_type='primary')
        download_table_button.on_click(self.download_matrix_ani)

        self.download_matrix_column_ani = pn.Column(self.save_matrix_file_path_ani,
                                                    download_table_button, pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_clustermap_column_ani,
                                 self.download_matrix_column_ani)

        # Transform the data into a scoring matrix
        pivoted_df = selected_genome_ani_df.pivot(columns='Sample1', index='Sample2', values='ANI')
        self.scores_matrix_ani = pivoted_df.combine_first(pivoted_df.T)
        np.fill_diagonal(self.scores_matrix_ani.values, 1.0)
        #print("\nScores matrix:")
        #print(self.scores_matrix_ani)

        # Check the number of columns in the matrix
        col_num = len(self.scores_matrix_ani.columns)
        #print("\ncreate_clustermap_pane: number of columns = " + str(col_num))

        # If the num of columns exceeds the defined maximum, do not create the clustermap plot
        # and display a message + a possibility to download the matrix
        if col_num > config.max_clustermap_cols:
            error = "The number of samples exceeds the limit of " + str(config.max_clustermap_cols) + \
                    " so the heatmap plot cannot be well presented."
            suggestion = "It ia possible to download the scoring matrix and display the heatmap in another program."
            clustermap_col = pn.Column(
                                       pn.pane.Markdown(error, styles={'font-size': "16px",
                                                                       'color': config.title_red_color,
                                                                       'margin-bottom': "0"}),
                                       pn.pane.Markdown(suggestion,
                                                        styles={'font-size': "14px", 'margin-top': "0"}),
                                       self.download_matrix_column_ani,
                                       styles={'padding': "15px"})

            self.clustermap_card_ani.append(clustermap_col)

        # Display the full clustermap pane, including the plot
        else:
            # There is metadata
            if self.is_metadata:
                # Update the color rows by- drop-down menu with the available metadata features
                self.color_by_feature_ani.options = self.metadata_features_list
                self.color_by_feature_ani.value = self.metadata_features_list[0]

                # Define watchers for the metadata widgets
                self.color_by_feature_clustermap_ani_watcher = self.color_by_feature_ani.param.watch(
                    self.set_not_continuous_clustermap_ani, 'value', onlychanged=True)
                self.continuous_clustermap_ani_watcher = self.is_continuous_clustermap_ani.param.watch(
                    self.change_continuous_state_clustermap_ani, 'value', onlychanged=True)
                self.feature_colormap_ani_watcher = \
                    self.feature_colormap_ani.param.watch(self.change_colormap_clustermap_ani,
                                                          'value', onlychanged=True)
                self.custom_colormap_clustermap_ani_watcher = \
                    self.custom_colormap_input_clustermap_ani.param.watch(self.get_custom_colormap_clustermap_ani,
                                                                          'value', onlychanged=True)

            # No metadata
            else:
                self.use_metadata_clustermap_ani.disabled = True

            self.clustermap_plot_ani = pn.bind(ps.create_clustermap, matrix=self.scores_matrix_ani, type="ANI",
                                               cmap=self.clustermap_cmap_ani, method=self.clustermap_method_ani,
                                               is_metadata=self.use_metadata_clustermap_ani,
                                               feature=self.color_by_feature_ani.value,
                                               is_continuous=self.is_continuous_clustermap_ani.value,
                                               cmap_metadata=self.feature_colormap_ani.value_name,
                                               custom_cmap=self.custom_colormap_input_clustermap_ani.value,
                                               metadata_dict=self.metadata_dict)

            self.clustermap_pane_ani = pn.pane.Matplotlib(self.clustermap_plot_ani, height=600, dpi=300, tight=True,
                                                          format='png')
            clustermap_row = pn.Row(controls_col, pn.Spacer(width=30), self.clustermap_pane_ani,
                                    styles={'padding': "15px"})

            self.clustermap_card_ani.append(clustermap_row)

    def change_colormap_clustermap_ani(self, event):
        self.update_clustermap_plot_ani()

    def get_custom_colormap_clustermap_ani(self, event):
        self.update_clustermap_plot_ani()

    def change_continuous_state_clustermap_ani(self, event):
        # Continuous feature
        if self.is_continuous_clustermap_ani.value:
            #print("\nIn change_continuous_state. Continuous feature")

            # Verify that the feature is indeed continuous
            feature = self.color_by_feature_ani.value
            unique_groups = sorted(
                list(set([self.metadata_dict[feature][sample] for sample in self.scores_matrix_ani.iloc[:, 0].index])))
            str_type = 0
            for group in unique_groups:
                if isinstance(group, str):
                    str_type = 1

            # Feature is not really continuous, treat as categorical
            if str_type == 1:
                #print("The feature is not really continuous - uncheck...")
                self.is_continuous_clustermap_ani.value = False

            # Feature is indeed really continuous
            else:
                self.feature_colormap_ani.options = config.continuous_colormap_dict
                self.feature_colormap_ani.value = config.continuous_colormap_dict['plasma_r']

        # Categorical feature
        else:
            #print("\nIn change_continuous_state. Categorical feature")
            self.feature_colormap_ani.options = config.categorical_colormap_dict
            self.feature_colormap_ani.value = config.categorical_colormap_dict['cet_glasbey']

    def set_not_continuous_clustermap_ani(self, event):
        #print("\nIn set_not_continuous")
        self.is_continuous_clustermap_ani.value = False
        self.update_clustermap_plot_ani()

    # Update the clustermap plot
    def update_clustermap_plot_ani(self):
        print("\nIn update_clustermap_plot_ani")
        self.clustermap_plot_ani = pn.bind(ps.create_clustermap, matrix=self.scores_matrix_ani, type="ANI",
                                           cmap=self.clustermap_cmap_ani, method=self.clustermap_method_ani,
                                           is_metadata=self.use_metadata_clustermap_ani,
                                           feature=self.color_by_feature_ani.value,
                                           is_continuous=self.is_continuous_clustermap_ani.value,
                                           cmap_metadata=self.feature_colormap_ani.value_name,
                                           custom_cmap=self.custom_colormap_input_clustermap_ani.value,
                                           metadata_dict=self.metadata_dict)

        self.clustermap_pane_ani.object = self.clustermap_plot_ani

    def download_clustermap_ani(self, event):
        fformat = self.clustermap_image_format_ani.value

        # Update the placeholder of the filename for download.
        clustermap_file = "ANI_Clustermap_" + self.ref_genome + "_" + self.clustermap_cmap_ani.value + "_" + \
                          self.clustermap_method_ani.value
        if self.use_metadata_clustermap_ani.value:
            if self.feature_colormap_ani.value_name == 'Define custom colormap':
                clustermap_file += "_colorby_" + self.color_by_feature_ani.value + "_custom_colormap"
            else:
                clustermap_file += "_colorby_" + self.color_by_feature_ani.value + "_" + \
                                   self.feature_colormap_ani.value_name
        self.save_clustermap_file_path_ani.placeholder = clustermap_file

        # Set the directory for saving
        if self.save_clustermap_file_path_ani.value == "":
            clustermap_file_path = self.downloads_dir_path + self.save_clustermap_file_path_ani.placeholder + "." \
                                   + fformat
        else:
            clustermap_file_path = self.save_clustermap_file_path_ani.value

            # Add a .png suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, clustermap_file_path, re.IGNORECASE):
                clustermap_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(clustermap_file_path):
                clustermap_file_path = self.downloads_dir_path + clustermap_file_path

        self.clustermap_plot_ani().savefig(clustermap_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + clustermap_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_clustermap_column_ani.pop(4)
        self.download_clustermap_column_ani.append(download_floatpanel)

    def download_matrix_ani(self, event):
        fformat = "txt"

        # Set the directory for saving
        if self.save_matrix_file_path_ani.value == "":
            matrix_path = self.downloads_dir_path + self.save_matrix_file_path_ani.placeholder + "." + fformat
        else:
            matrix_path = self.save_matrix_file_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, matrix_path, re.IGNORECASE):
                matrix_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(matrix_path):
                matrix_path = self.downloads_dir_path + matrix_path

        self.scores_matrix_ani.to_csv(matrix_path, sep='\t')

        download_message = "The table is successfully saved under:\n" + matrix_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_matrix_column_ani.pop(2)
        self.download_matrix_column_ani.append(download_floatpanel)

    def change_continuous_state_network(self, event):
        # Continuous feature
        if self.is_continuous_network.value:
            #print("\nIn change_continuous_state. Continuous feature")

            # Verify that the feature is indeed continuous
            nodes_feature = self.nodes_color_by.value
            unique_groups = list(set([self.network.nodes[node][nodes_feature] for node in self.network.nodes()]))
            str_type = 0
            for group in unique_groups:
                if isinstance(group, str):
                    str_type = 1

            # Feature is not really continuous, treat as categorical
            if str_type == 1:
                #print("The feature is not really continuous - uncheck...")
                self.is_continuous_network.value = False

            # Feature is indeed really continuous
            else:
                self.nodes_colormap.options = config.continuous_colormap_dict
                self.nodes_colormap.value = config.continuous_colormap_dict['plasma_r']

        # Categorical feature
        else:
            #print("\nIn change_continuous_state. Categorical feature")
            self.nodes_colormap.options = config.categorical_colormap_dict
            self.nodes_colormap.value = config.categorical_colormap_dict['cet_glasbey']

    def change_colormap_network(self, event):
        #print("\nIn change_colormap. Continuous state = " + str(self.is_continuous.value))
        self.update_network_plot()

    def get_custom_colormap_network(self, event):
        #print("\nIn change_colormap. Continuous state = " + str(self.is_continuous.value))
        self.update_network_plot()

    def set_not_continuous_network(self, event):
        #print("\nIn set_not_continuous")
        self.is_continuous_network.value = False
        self.update_network_plot()

    def fill_feature_groups(self, event):
        feature = self.nodes_highlight_by.value
        unique_groups = sorted(list(set([str(self.network.nodes[node][feature]) for node in self.network.nodes()])))
        # print("\nUnique groups:")
        # print(unique_groups)
        if 'nan' in unique_groups:
            unique_groups.remove('nan')
            unique_groups.append('nan')
        self.nodes_highlight_group_select.options = unique_groups
        self.nodes_highlight_group_select.value = unique_groups[0]

        self.update_network_plot()

    def create_network_pane(self, selected_genome_and_size_avg_df):
        mean_only = 0
        mean_std = 0

        init_button = pn.widgets.Button(name='Initialize nodes positions', button_type='primary',
                                        button_style='outline')
        init_button.on_click(self.init_positions)
        init_button_row = pn.Row(init_button, align='center')

        styling_title = "Network customization options:"
        nodes_cust_title = pn.pane.Markdown("Nodes customizations:",
                                            styles={'font-size': "15px", 'font-weight': "bold", 'margin': "0 0 0 15px"})
        edges_cust_title = pn.pane.Markdown("Edges customizations:",
                                            styles={'font-size': "15px", 'font-weight': "bold", 'margin': "0 0 0 15px"})
        no_metadata_colors_row = pn.Row(self.network_node_color, pn.Spacer(width=10), self.network_edge_color)
        continuous_col = pn.Column(pn.Spacer(height=20), self.is_continuous_network)
        nodes_color_by_row = pn.Row(self.nodes_color_by, continuous_col, styles={'margin': "0 0 0 8px"})
        edges_color_by_row = pn.Row(self.edges_color_by, self.network_within_color, pn.Spacer(width=3),
                                    self.network_between_color, styles={'margin': "0 0 10px 5px"})
        highlight_nodes_row = pn.Row(self.nodes_highlight_by, self.nodes_highlight_group_select,
                                     styles={'margin': "0 0 0 5px"})
        metadata_coloring_col = pn.Column(nodes_cust_title,
                                          nodes_color_by_row,
                                          self.nodes_colormap,
                                          self.custom_colormap_input,
                                          pn.Spacer(height=15),
                                          self.highlight_nodes_by_feature,
                                          highlight_nodes_row,
                                          pn.Spacer(height=10),
                                          edges_cust_title,
                                          self.color_edges_by_feature,
                                          edges_color_by_row)
        self.metadata_colorby_card.append(metadata_coloring_col)
        network_threshold_row = pn.Row(self.network_threshold_select, self.network_threshold_input)
        params_col = pn.Column(network_threshold_row,
                               self.network_iterations,
                               init_button_row)
        self.layout_parameters_card.append(params_col)
        highlight_sample_row = pn.Row(self.highlight_sample_chkbox, self.highlight_sample_input)
        show_labels_col = pn.Column(pn.Spacer(height=5), self.show_labels_chkbox)
        show_labels_row = pn.Row(show_labels_col, self.all_or_highlighted_select)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.layout_parameters_card,
                                pn.Spacer(height=5),
                                no_metadata_colors_row,
                                pn.Spacer(height=5),
                                self.use_metadata_network,
                                self.metadata_colorby_card,
                                pn.Spacer(height=5),
                                highlight_sample_row,
                                show_labels_row
                                )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_network)

        self.download_network_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                  styles={'font-size': "15px",
                                                                          'font-weight': "bold",
                                                                          'color': config.title_blue_color,
                                                                          'margin': "0"}),
                                                 self.network_image_format, self.save_network_plot_path,
                                                 download_button, pn.pane.Markdown())
        download_table_button = pn.widgets.Button(name='Download network data in tsv format', button_type='primary')
        download_table_button.on_click(self.download_network_table)

        self.download_network_table_column = pn.Column(self.save_network_table_path, download_table_button,
                                                       pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_network_column,
                                 self.download_network_table_column)

        ########################################################
        # Create a table for the network
        self.df_for_network = selected_genome_and_size_avg_df.loc[:, ['Sample1', 'Sample2', 'APSS']].copy()

        eps = 1e-6
        self.df_for_network.loc[(self.df_for_network['APSS'] < 0), 'APSS'] = 0
        self.df_for_network.loc[(self.df_for_network['APSS'] == 1), 'APSS'] = 1 - eps

        # Set a score threshold for the connections (below it zero the weight).
        # The default threshold is the mean APSS + 1std
        mean_APSS = self.df_for_network.loc[:, 'APSS'].mean().round(2)
        std_APSS = self.df_for_network.loc[:, 'APSS'].std().round(2)

        print("\ncreate_network_pane:")
        print("Mean APSS: " + str(mean_APSS))
        print("Standard deviation APSS: " + str(std_APSS) + "\n")

        if mean_APSS <= 0.99:
            self.APSS_connections_threshold = mean_APSS
        else:
            self.APSS_connections_threshold = 0.99

        # First stage: zero the scores below the threshold
        self.df_for_network['filtered_score'] = self.df_for_network['APSS'].where(self.df_for_network['APSS'] >=
                                                                                  self.APSS_connections_threshold, 0)

        # Second stage: perform min-max normalization on the scores above the threshold
        mask = self.df_for_network['filtered_score'] > 0
        min_val = self.df_for_network.loc[mask, 'filtered_score'].min()
        max_val = self.df_for_network.loc[mask, 'filtered_score'].max()
        print("min_value = " + str(min_val) + ", max_value = " + str(max_val))
        self.df_for_network['norm_score'] = 0.0
        # avoid division by zero if all remaining values are equal
        if max_val > min_val:
            self.df_for_network.loc[mask, 'norm_score'] = (self.df_for_network.loc[mask, 'filtered_score'] - min_val) / \
                                                          (max_val - min_val)
            # replace exact zeros and ones to avoid later problems in the log function
            self.df_for_network.loc[mask, 'norm_score'] = \
                self.df_for_network.loc[mask, 'norm_score'].clip(eps, 1 - eps)
        else:
            self.df_for_network.loc[mask, 'norm_score'] = max_val

        # Third stage: apply the log function to transform the scores into weight values
        self.df_for_network['weight'] = 0.0
        self.df_for_network.loc[mask, 'weight'] = np.negative(np.log(1 - self.df_for_network.loc[mask, 'norm_score']))

        # Forth stage: calculate line width from weights (normalize values between 0.5 and 4)
        self.df_for_network['width'] = 0.0
        min_width = self.df_for_network.loc[mask, 'weight'].min()
        max_width = self.df_for_network.loc[mask, 'weight'].max()
        new_min_width = config.min_edge_width
        new_max_width = config.max_edge_width
        self.df_for_network.loc[mask, 'width'] = new_min_width + \
                                                 (self.df_for_network.loc[mask, 'weight'] - min_width) * \
                                                 (new_max_width - new_min_width) / (max_width - min_width)

        #print("\nDF for network with weights and widths:")
        #print(self.df_for_network)

        # Update the placeholder of the filenames for download with the default threshold.
        network_file = "Network_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                       self.network_iterations.value + "_iterations_threshold_" + str(self.APSS_connections_threshold)
        self.save_network_plot_path.placeholder = network_file
        table_file = "Network_" + self.ref_genome + "_" + self.sampling_size + "_regions_threshold_" + \
                     str(self.APSS_connections_threshold)
        self.save_network_table_path.placeholder = table_file

        # Create a network using networkx
        self.network = nx.from_pandas_edgelist(self.df_for_network, source='Sample1', target='Sample2',
                                               edge_attr=['weight', 'width'])
        self.nodes_list = list(self.network.nodes)
        nodes_num = len(self.nodes_list)
        #print("\nNumber of nodes in the network = " + str(nodes_num))

        # If the number of nodes in the network exceeds the defined maximum, do not create the plot
        # and display only a message + a possibility to download the network data in tsv format
        if nodes_num > config.max_network_nodes:
            error = "The number of samples exceeds the limit of 300, so the network cannot be presented with all its " \
                    "interactive features."
            suggestion = "It ia possible to download the network data in tsv format and visualize it using another " \
                         "program."
            network_col = pn.Column(
                pn.pane.Markdown(error, styles={'font-size': "16px",
                                                'color': config.title_red_color,
                                                'margin-bottom': "0"}),
                pn.pane.Markdown(suggestion,
                                 styles={'font-size': "14px", 'margin-top': "0"}),
                self.download_network_table_column,
                styles={'padding': "15px"})

            self.network_card.append(network_col)

        # Display the full network pane, including the plot
        else:
            # Initialize the network positions
            self.generate_rand_positions()

            # Add the actual threshold value to the network_threshold_select widget
            self.network_threshold_input.disabled = True
            self.network_threshold_select.options = []

            if mean_APSS > 0.99:
                mean_APSS = 0.99
                mean_only = 1

            else:
                mean_std = round((mean_APSS + std_APSS), 2)
                if mean_std > 0.99:
                    mean_std = 0.99

            str_mean = config.network_thresholds_options[0] + " (APSS=" + str(mean_APSS) + ")"
            self.network_threshold_select.options.append(str_mean)
            if not mean_only:
                str_mean_std = config.network_thresholds_options[1] + " (APSS=" + str(mean_std) + ")"
                self.network_threshold_select.options.append(str_mean_std)

            self.network_threshold_select.options.append(config.network_thresholds_options[2])
            # Set the mean as the default threshold
            self.network_threshold_select.value = self.network_threshold_select.options[0]

            # Set watchers for the threshold widgets
            self.threshold_select_watcher = self.network_threshold_select.param.watch(partial(
                self.changed_threshold_select, mean_APSS, mean_std, mean_only), 'value',
                onlychanged=True)
            self.threshold_input_watcher = self.network_threshold_input.param.watch(self.changed_threshold_input,
                                                                                    'value', onlychanged=True)
            self.highlight_sample_watcher = self.highlight_sample_input.param.watch(self.change_highlighted_sample,
                                                                                    'value', onlychanged=True)

            # There is metadata
            if self.is_metadata:

                # Insert the features information as nodes attributes
                for node in self.nodes_list:
                    for feature in self.metadata_features_list:
                        self.network.nodes[node][feature] = self.metadata_dict[feature][node]

                    # Add node attribute 'SampleID' for the hover tooltip
                    self.network.nodes[node]['SampleID'] = node

                # Update the color by- drop-down menus with the available metadata features
                first_feature = self.metadata_features_list[0]
                self.nodes_color_by.options = self.metadata_features_list
                self.nodes_color_by.value = first_feature
                self.edges_color_by.options = self.metadata_features_list
                self.edges_color_by.value = first_feature
                self.nodes_highlight_by.options = self.metadata_features_list
                self.nodes_highlight_by.value = first_feature
                # Fill the groups for the first feature
                unique_groups = sorted(list(set([str(self.network.nodes[node][first_feature]) for node in
                                                 self.network.nodes()])))
                if 'nan' in unique_groups:
                    unique_groups.remove('nan')
                    unique_groups.append('nan')
                self.nodes_highlight_group_select.options = unique_groups
                self.nodes_highlight_group_select.value = unique_groups[0]

                self.nodes_colorby_watcher = self.nodes_color_by.param.watch(self.set_not_continuous_network, 'value',
                                                                             onlychanged=True)
                self.continuous_network_watcher = self.is_continuous_network.param.watch(
                    self.change_continuous_state_network, 'value', onlychanged=True)
                self.colormap_watcher = self.nodes_colormap.param.watch(self.change_colormap_network, 'value',
                                                                        onlychanged=True)
                self.custom_colormap_watcher = self.custom_colormap_input.param.watch(self.get_custom_colormap_network,
                                                                                      'value', onlychanged=True)
                self.nodes_highlight_by_watcher = self.nodes_highlight_by.param.watch(self.fill_feature_groups,
                                                                                      'value', onlychanged=True)

            # No metadata
            else:
                self.use_metadata_network.disabled = True

                # Add node attribute 'SampleID' for the hover tooltip
                for node in self.nodes_list:
                    self.network.nodes[node]['SampleID'] = node

            # Create the network plot using the selected parameters
            self.network_plot_hv = pn.bind(ps.cretae_network_plot, network=self.network,
                                           is_metadata=self.use_metadata_network,
                                           nodes_feature=self.nodes_color_by.value,
                                           is_continuous=self.is_continuous_network.value,
                                           cmap=self.nodes_colormap.value,
                                           custom_cmap=self.custom_colormap_input.value,
                                           node_color=self.network_node_color,
                                           edge_color=self.network_edge_color,
                                           is_highlight_group=self.highlight_nodes_by_feature,
                                           highlight_feature=self.nodes_highlight_by.value,
                                           highlight_group=self.nodes_highlight_group_select,
                                           is_edge_colorby=self.color_edges_by_feature,
                                           edges_feature=self.edges_color_by,
                                           within_edge_color=self.network_within_color,
                                           between_edge_color=self.network_between_color,
                                           iterations=self.network_iterations, pos_dict=self.pos_dict,
                                           show_labels=self.show_labels_chkbox,
                                           all_or_highlighted=self.all_or_highlighted_select,
                                           is_highlight_samples=self.highlight_sample_chkbox,
                                           samples_to_highlight=self.highlight_sample_input.value,
                                           metadata_dict=self.metadata_dict)
            self.network_pane = pn.pane.HoloViews(self.network_plot_hv, height=600, width=700, sizing_mode="fixed")

            network_row = pn.Row(controls_col, pn.Spacer(width=15), self.network_pane, styles={'padding': "15px"})
            self.network_card.append(network_row)

    def download_network(self, event):
        fformat = self.network_image_format.value

        # Update the placeholder of the filenames for download with the default threshold.
        if self.use_metadata_network.value:
            if self.nodes_colormap.value_name == 'Define custom colormap':
                network_file = "Network_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                               self.network_iterations.value + "_iterations_threshold_" + \
                               str(self.APSS_connections_threshold) + "_colorby_" + self.nodes_color_by.value + \
                               "_custom_colormap"
            else:
                network_file = "Network_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                               self.network_iterations.value + "_iterations_threshold_" + \
                               str(self.APSS_connections_threshold) + "_colorby_" + self.nodes_color_by.value + "_" + \
                               self.nodes_colormap.value_name
        else:
            network_file = "Network_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                           self.network_iterations.value + "_iterations_threshold_" + \
                           str(self.APSS_connections_threshold)
        if self.show_labels_chkbox.value:
            network_file += "_with_labels"

        self.save_network_plot_path.placeholder = network_file
        table_file = "Network_" + self.ref_genome + "_" + self.sampling_size + "_regions_threshold_" + \
                     str(self.APSS_connections_threshold)
        self.save_network_table_path.placeholder = table_file

        # Set the directory for saving
        if self.save_network_plot_path.value == "":
            network_file_path = self.downloads_dir_path + self.save_network_plot_path.placeholder + "." + fformat
        else:
            network_file_path = self.save_network_plot_path.value

            # Add a .png suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, network_file_path, re.IGNORECASE):
                network_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(network_file_path):
                network_file_path = self.downloads_dir_path + network_file_path

        # Get the updated matplotlib plot
        self.network_plot_matplotlib = pn.bind(ps.cretae_network_plot_matplotlib, network=self.network,
                                               is_metadata=self.use_metadata_network,
                                               nodes_feature=self.nodes_color_by.value,
                                               is_continuous=self.is_continuous_network.value,
                                               cmap=self.nodes_colormap.value_name,
                                               custom_cmap=self.custom_colormap_input.value,
                                               node_color=self.network_node_color,
                                               edge_color=self.network_edge_color,
                                               is_highlight_group=self.highlight_nodes_by_feature,
                                               highlight_feature=self.nodes_highlight_by.value,
                                               highlight_group=self.nodes_highlight_group_select,
                                               is_edge_colorby=self.color_edges_by_feature,
                                               edges_feature=self.edges_color_by,
                                               within_edge_color=self.network_within_color,
                                               between_edge_color=self.network_between_color,
                                               iterations=self.network_iterations, pos_dict=self.pos_dict,
                                               show_labels=self.show_labels_chkbox,
                                               all_or_highlighted=self.all_or_highlighted_select,
                                               is_highlight_samples=self.highlight_sample_chkbox,
                                               samples_to_highlight=self.highlight_sample_input.value,
                                               metadata_dict=self.metadata_dict)

        # Save the network plot in the requested format using matplotlib
        self.network_plot_matplotlib().savefig(network_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + network_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_network_column.pop(4)
        self.download_network_column.append(download_floatpanel)

    def download_network_table(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_network_table_path.value == "":
            network_table_path = self.downloads_dir_path + self.save_network_table_path.placeholder + "." + fformat
        else:
            network_table_path = self.save_network_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, network_table_path, re.IGNORECASE):
                network_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(network_table_path):
                network_table_path = self.downloads_dir_path + network_table_path

        # Export just the necessary columns to csv
        self.df_for_network[['Sample1', 'Sample2', 'APSS', 'weight']].to_csv(network_table_path, index=False, sep='\t')

        download_message = "The table is successfully saved under:\n" + network_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_network_table_column.pop(2)
        self.download_network_table_column.append(download_floatpanel)

    def init_positions(self, event):
        self.generate_rand_positions()
        self.update_network_plot()

    def generate_rand_positions(self):
        self.pos_dict = {}

        for node in self.nodes_list:
            pos_x = generate_rand_pos()
            pos_y = generate_rand_pos()
            pos_tuple = (pos_x, pos_y)
            self.pos_dict[node] = pos_tuple

    def changed_threshold_select(self, mean, mean_std, mean_only, event):
        #print("\nchanged_threshold_select:")
        #print("Current mean: " + str(mean))

        # The mean option
        if self.network_threshold_select.value == self.network_threshold_select.options[0]:
            self.APSS_connections_threshold = mean
            self.network_threshold_input.disabled = True

        elif self.network_threshold_select.value == self.network_threshold_select.options[1]:
            # The second option is the custom threshold
            if mean_only:
                self.APSS_connections_threshold = self.network_threshold_input.value
                self.network_threshold_input.disabled = False
            # The second option is the mean + 1std
            else:
                self.APSS_connections_threshold = mean_std
                self.network_threshold_input.disabled = True

        # The last option is custom threshold
        else:
            self.APSS_connections_threshold = self.network_threshold_input.value
            self.network_threshold_input.disabled = False

        self.change_weight_attribute()

    def changed_threshold_input(self, event):
        #print("\nIn changed_threshold_input")
        self.APSS_connections_threshold = self.network_threshold_input.value
        self.change_weight_attribute()

    def change_weight_attribute(self):

        self.APSS_connections_threshold = round(self.APSS_connections_threshold, 2)
        print("\nchange_weight_attribute:")
        print("APSS_connections_threshold = " + str(self.APSS_connections_threshold))

        # Update the threshold in the deafult filenames for download
        network_file = "Network_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions_" + \
                       self.network_iterations.value + "_iterations_threshold_" + str(self.APSS_connections_threshold)
        self.save_network_plot_path.placeholder = network_file
        table_file = "Network_" + self.ref_genome + "_" + self.sampling_size + "_regions_threshold_" + \
                     str(self.APSS_connections_threshold)
        self.save_network_table_path.placeholder = table_file

        ## Recalculate the weights
        # First stage: zero the scores below the threshold
        self.df_for_network['filtered_score'] = self.df_for_network['APSS'].where(self.df_for_network['APSS'] >=
                                                                                  self.APSS_connections_threshold, 0)

        # Second stage: perform min-max normalization on the scores above the threshold
        mask = self.df_for_network['filtered_score'] > 0
        min_val = self.df_for_network.loc[mask, 'filtered_score'].min()
        max_val = self.df_for_network.loc[mask, 'filtered_score'].max()
        self.df_for_network['norm_score'] = 0.0
        # avoid division by zero if all remaining values are equal
        eps = 1e-6
        if max_val > min_val:
            self.df_for_network.loc[mask, 'norm_score'] = (self.df_for_network.loc[mask, 'filtered_score'] - min_val) / \
                                                          (max_val - min_val)
            # replace exact zeros (min_val) with 10-6 (so that they will have a very low weight but not zero)
            self.df_for_network.loc[mask, 'norm_score'] = \
                self.df_for_network.loc[mask, 'norm_score'].clip(eps, 1 - eps)
        else:
            self.df_for_network.loc[mask, 'norm_score'] = max_val

        # Third stage: apply the log function to transform the scores into weight values
        self.df_for_network['weight'] = 0.0
        self.df_for_network.loc[mask, 'weight'] = np.negative(np.log(1 - self.df_for_network.loc[mask, 'norm_score']))

        # Forth stage: calculate line width from weights (normalize values between 0.5 and 4)
        self.df_for_network['width'] = 0.0
        min_width = self.df_for_network.loc[mask, 'weight'].min()
        max_width = self.df_for_network.loc[mask, 'weight'].max()
        new_min_width = config.min_edge_width
        new_max_width = config.max_edge_width
        self.df_for_network.loc[mask, 'width'] = new_min_width + \
                                                 (self.df_for_network.loc[mask, 'weight'] - min_width) * \
                                                 (new_max_width - new_min_width) / (max_width - min_width)

        #print("\nDF for network with weights and widths:")
        #print(self.df_for_network)

        before = time.time()
        # Create a new network from the updated df using networkx
        self.network = nx.from_pandas_edgelist(self.df_for_network, source='Sample1', target='Sample2',
                                               edge_attr=['weight', 'width'])
        self.nodes_list = list(self.network.nodes)

        if self.is_metadata:
            # Insert the features information as nodes attributes to the new network
            for node in self.nodes_list:
                for feature in self.metadata_features_list:
                    self.network.nodes[node][feature] = self.metadata_dict[feature][node]

                # Add node attribute 'SampleID' for the hover tooltip
                self.network.nodes[node]['SampleID'] = node

        after = time.time()
        duration = after - before
        #print("Creating a new network took " + str(duration) + " seconds.\n")

        before = time.time()
        self.update_network_plot()
        after = time.time()
        duration = after - before
        print("Updating the network plot took " + str(duration) + " seconds.\n")

    def change_highlighted_sample(self, event):
        print("\nIn change_highlighted_sample")
        print("Samples to highlight: " + self.highlight_sample_input.value)
        self.update_network_plot()

    # Update the network plot using the selected parameters and the new positions dict
    def update_network_plot(self):
        #print("\nIn update_network_plot")
        self.network_plot_hv = pn.bind(ps.cretae_network_plot, network=self.network,
                                       is_metadata=self.use_metadata_network, nodes_feature=self.nodes_color_by.value,
                                       is_continuous=self.is_continuous_network.value, cmap=self.nodes_colormap.value,
                                       custom_cmap=self.custom_colormap_input.value,
                                       node_color=self.network_node_color, edge_color=self.network_edge_color,
                                       is_highlight_group=self.highlight_nodes_by_feature,
                                       highlight_feature=self.nodes_highlight_by.value,
                                       highlight_group=self.nodes_highlight_group_select,
                                       is_edge_colorby=self.color_edges_by_feature, edges_feature=self.edges_color_by,
                                       within_edge_color=self.network_within_color,
                                       between_edge_color=self.network_between_color,
                                       iterations=self.network_iterations, pos_dict=self.pos_dict,
                                       show_labels=self.show_labels_chkbox,
                                       all_or_highlighted=self.all_or_highlighted_select,
                                       is_highlight_samples=self.highlight_sample_chkbox,
                                       samples_to_highlight=self.highlight_sample_input.value,
                                       metadata_dict=self.metadata_dict)

        self.network_pane.object = self.network_plot_hv
    def change_continuous_state_network_ani(self, event):
        # Continuous feature
        if self.is_continuous_network_ani.value:
            #print("\nIn change_continuous_state. Continuous feature")

            # Verify that the feature is indeed continuous
            nodes_feature = self.nodes_color_by_ani.value
            unique_groups = list(set([self.network_ani.nodes[node][nodes_feature] for node in self.network_ani.nodes()]))
            str_type = 0
            for group in unique_groups:
                if isinstance(group, str):
                    str_type = 1

            # Feature is not really continuous, treat as categorical
            if str_type == 1:
                #print("The feature is not really continuous - uncheck...")
                self.is_continuous_network_ani.value = False

            # Feature is indeed really continuous
            else:
                self.nodes_colormap_ani.options = config.continuous_colormap_dict
                self.nodes_colormap_ani.value = config.continuous_colormap_dict['plasma_r']

        # Categorical feature
        else:
            #print("\nIn change_continuous_state. Categorical feature")
            self.nodes_colormap_ani.options = config.categorical_colormap_dict
            self.nodes_colormap_ani.value = config.categorical_colormap_dict['cet_glasbey']

    def change_colormap_network_ani(self, event):
        #print("\nIn change_colormap. Continuous state = " + str(self.is_continuous.value))
        self.update_network_plot_ani()

    def get_custom_colormap_network_ani(self, event):
        #print("\nIn change_colormap. Continuous state = " + str(self.is_continuous.value))
        self.update_network_plot_ani()

    def set_not_continuous_network_ani(self, event):
        #print("\nIn set_not_continuous")
        self.is_continuous_network_ani.value = False
        self.update_network_plot_ani()

    def fill_feature_groups_ani(self, event):
        feature = self.nodes_highlight_by_ani.value
        unique_groups = sorted(list(set([str(self.network_ani.nodes[node][feature]) for node in self.network_ani.nodes()])))
        # print("\nUnique groups:")
        # print(unique_groups)
        if 'nan' in unique_groups:
            unique_groups.remove('nan')
            unique_groups.append('nan')
        self.nodes_highlight_group_select_ani.options = unique_groups
        self.nodes_highlight_group_select_ani.value = unique_groups[0]

        self.update_network_plot_ani()

    def create_network_pane_ani(self, ani_scores_selected_genome_df):
        mean_only = 0
        mean_std = 0

        init_button = pn.widgets.Button(name='Initialize nodes positions', button_type='primary',
                                        button_style='outline')
        init_button.on_click(self.init_positions_ani)
        init_button_row = pn.Row(init_button, align='center')

        styling_title = "Network customization options:"
        nodes_cust_title = pn.pane.Markdown("Nodes customizations:",
                                            styles={'font-size': "15px", 'font-weight': "bold", 'margin': "0 0 0 15px"})
        edges_cust_title = pn.pane.Markdown("Edges customizations:",
                                            styles={'font-size': "15px", 'font-weight': "bold", 'margin': "0 0 0 15px"})
        no_metadata_colors_row = pn.Row(self.network_node_color_ani, pn.Spacer(width=10), self.network_edge_color_ani)
        continuous_col = pn.Column(pn.Spacer(height=20), self.is_continuous_network_ani)
        nodes_color_by_row = pn.Row(self.nodes_color_by_ani, continuous_col, styles={'margin': "0 0 0 8px"})
        edges_color_by_row = pn.Row(self.edges_color_by_ani, self.network_within_color_ani, pn.Spacer(width=3),
                                    self.network_between_color_ani, styles={'margin': "0 0 10px 5px"})
        highlight_nodes_row = pn.Row(self.nodes_highlight_by_ani, self.nodes_highlight_group_select_ani,
                                     styles={'margin': "0 0 0 5px"})
        metadata_coloring_col = pn.Column(nodes_cust_title,
                                          nodes_color_by_row,
                                          self.nodes_colormap_ani,
                                          self.custom_colormap_input_ani,
                                          pn.Spacer(height=15),
                                          self.highlight_nodes_by_feature_ani,
                                          highlight_nodes_row,
                                          pn.Spacer(height=10),
                                          edges_cust_title,
                                          self.color_edges_by_feature_ani,
                                          edges_color_by_row)
        self.metadata_colorby_card_ani.append(metadata_coloring_col)
        network_threshold_row = pn.Row(self.network_threshold_select_ani, self.network_threshold_input_ani)
        params_col = pn.Column(network_threshold_row,
                               self.network_iterations_ani,
                               init_button_row)
        self.layout_parameters_card_ani.append(params_col)
        highlight_sample_row = pn.Row(self.highlight_sample_chkbox_ani, self.highlight_sample_input_ani)
        show_labels_col = pn.Column(pn.Spacer(height=5), self.show_labels_chkbox_ani)
        show_labels_row = pn.Row(show_labels_col, self.all_or_highlighted_select_ani)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.layout_parameters_card_ani,
                                pn.Spacer(height=5),
                                no_metadata_colors_row,
                                pn.Spacer(height=5),
                                self.use_metadata_network_ani,
                                self.metadata_colorby_card_ani,
                                pn.Spacer(height=5),
                                highlight_sample_row,
                                show_labels_row
                                )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_network_ani)

        self.download_network_column_ani = pn.Column(pn.pane.Markdown(save_file_title,
                                                                      styles={'font-size': "15px",
                                                                              'font-weight': "bold",
                                                                              'color': config.title_blue_color,
                                                                              'margin': "0"}),
                                                     self.network_image_format_ani, self.save_network_plot_path_ani,
                                                     download_button, pn.pane.Markdown())
        download_table_button = pn.widgets.Button(name='Download network data in tsv format', button_type='primary')
        download_table_button.on_click(self.download_network_table_ani)

        self.download_network_table_column_ani = pn.Column(self.save_network_table_path_ani, download_table_button,
                                                           pn.pane.Markdown())

        controls_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_network_column_ani,
                                 self.download_network_table_column_ani)

        ########################################################
        # Create a table for the network
        self.df_for_network_ani = ani_scores_selected_genome_df.loc[:, ['Sample1', 'Sample2', 'ANI']].copy()

        eps = 1e-6
        self.df_for_network_ani.loc[(self.df_for_network_ani['ANI'] < 0), 'ANI'] = 0
        self.df_for_network_ani.loc[(self.df_for_network_ani['ANI'] == 1), 'ANI'] = 1 - eps

        # Set a score threshold for the connections (below it zero the weight).
        # The default threshold is the mean APSS + 1std
        mean_ANI = self.df_for_network_ani.loc[:, 'ANI'].mean().round(3)
        std_ANI = self.df_for_network_ani.loc[:, 'ANI'].std().round(3)

        # Set the default connections threshold
        if mean_ANI <= 0.999:
            self.ani_connections_threshold = mean_ANI
        else:
            self.ani_connections_threshold = 0.999

        # First stage: zero the scores below the threshold
        self.df_for_network_ani['filtered_score'] = \
            self.df_for_network_ani['ANI'].where(self.df_for_network_ani['ANI'] >= self.ani_connections_threshold, 0)

        # Second stage: perform min-max normalization on the scores above the threshold
        mask = self.df_for_network_ani['filtered_score'] > 0
        min_val = self.df_for_network_ani.loc[mask, 'filtered_score'].min()
        max_val = self.df_for_network_ani.loc[mask, 'filtered_score'].max()
        self.df_for_network_ani['norm_score'] = 0.0
        # avoid division by zero if all remaining values are equal
        if max_val > min_val:
            self.df_for_network_ani.loc[mask, 'norm_score'] = \
                (self.df_for_network_ani.loc[mask, 'filtered_score'] - min_val) / (max_val - min_val)
            # replace exact zeros (min_val) with 10-6 (so that they will have a very low weight but not zero)
            self.df_for_network_ani.loc[mask, 'norm_score'] = \
                self.df_for_network_ani.loc[mask, 'norm_score'].clip(eps, 1 - eps)
        else:
            self.df_for_network_ani.loc[mask, 'norm_score'] = max_val

        # Third stage: apply the log function to transform the scores into weight values
        self.df_for_network_ani['weight'] = 0.0
        self.df_for_network_ani.loc[mask, 'weight'] = \
            np.negative(np.log(1 - self.df_for_network_ani.loc[mask, 'norm_score']))

        # Forth stage: calculate line width from weights (normalize values between 0.5 and 4)
        self.df_for_network_ani['width'] = 0.0
        min_width = self.df_for_network_ani.loc[mask, 'weight'].min()
        max_width = self.df_for_network_ani.loc[mask, 'weight'].max()
        new_min_width = config.min_edge_width
        new_max_width = config.max_edge_width
        self.df_for_network_ani.loc[mask, 'width'] = new_min_width + \
                                                 (self.df_for_network_ani.loc[mask, 'weight'] - min_width) * \
                                                 (new_max_width - new_min_width) / (max_width - min_width)

        #print("\nDF for network with weights and widths:")
        #print(self.df_for_network_ani)
        print("\ncreate_network_pane_ani:")
        print("Mean ANI: " + str(mean_ANI))
        print("Standard deviation ANI: " + str(std_ANI) + "\n")

        # Update the placeholder of the filenames for download with the default threshold.
        network_file = "Network_plot_ANI_" + self.ref_genome + "_" + self.network_iterations_ani.value + \
                       "_iterations_threshold_" + str(self.ani_connections_threshold)
        self.save_network_plot_path_ani.placeholder = network_file
        table_file = "Network_" + self.ref_genome + "_threshold_" + str(self.ani_connections_threshold)
        self.save_network_table_path_ani.placeholder = table_file

        # Create a network using networkx
        self.network_ani = nx.from_pandas_edgelist(self.df_for_network_ani, source='Sample1', target='Sample2',
                                                   edge_attr=['weight', 'width'])
        self.nodes_list_ani = list(self.network_ani.nodes)
        nodes_num = len(self.nodes_list_ani)
        #print("\nNumber of nodes in the network = " + str(nodes_num))

        # If the number of nodes in the network exceeds the defined maximum, do not create the plot
        # and display only a message + a possibility to download the network data in tsv format
        if nodes_num > config.max_network_nodes:
            error = "The number of samples exceeds the limit of 300, so the network cannot be presented with all its " \
                    "interactive features."
            suggestion = "It ia possible to download the network data in tsv format and visualize it using another " \
                         "program."
            network_col = pn.Column(
                pn.pane.Markdown(error, styles={'font-size': "16px",
                                                'color': config.title_red_color,
                                                'margin-bottom': "0"}),
                pn.pane.Markdown(suggestion,
                                 styles={'font-size': "14px", 'margin-top': "0"}),
                self.download_network_table_column_ani,
                styles={'padding': "15px"})

            self.network_card_ani.append(network_col)

        # Display the full network pane, including the plot
        else:
            # Initialize the network positions
            self.generate_rand_positions_ani()

            # Add the actual threshold value to the network_threshold_select widget
            self.network_threshold_input_ani.disabled = True
            self.network_threshold_select_ani.options = []

            if mean_ANI >= 0.999:
                mean_ANI = 0.999
                str_mean = config.network_thresholds_options_ani[0] + " (ANI=0.999)"
                self.network_threshold_select_ani.options.append(str_mean)
                mean_only = 1

            else:
                str_mean = config.network_thresholds_options_ani[0] + " (ANI=" + str(mean_ANI) + ")"
                self.network_threshold_select_ani.options.append(str_mean)

                mean_std = round((mean_ANI + std_ANI), 3)
                if mean_std >= 0.999:
                    mean_std = 0.999
                    str_mean_std = config.network_thresholds_options_ani[1] + " (ANI=0.999)"
                    self.network_threshold_select_ani.options.append(str_mean_std)

                else:
                    str_mean_std = config.network_thresholds_options_ani[1] + " (ANI=" + str(mean_std) + ")"
                    self.network_threshold_select_ani.options.append(str_mean_std)

            self.network_threshold_select_ani.options.append(config.network_thresholds_options_ani[2])
            # Set the mean as the default threshold
            self.network_threshold_select_ani.value = self.network_threshold_select_ani.options[0]

            # Set watchers for the threshold widgets
            self.threshold_select_ani_watcher = self.network_threshold_select_ani.param.watch(partial(
                self.changed_threshold_select_ani, mean_ANI, mean_std, mean_only), 'value',
                onlychanged=True)
            self.is_defined_threshold_select_ani_watcher = 1
            self.threshold_input_ani_watcher = \
                self.network_threshold_input_ani.param.watch(self.changed_threshold_input_ani,
                                                             'value', onlychanged=True)
            self.highlight_sample_ani_watcher = \
                self.highlight_sample_input_ani.param.watch(self.change_highlighted_sample_ani,
                                                            'value', onlychanged=True)

            # There is metadata
            if self.is_metadata:

                # Insert the features information as nodes attributes
                for node in self.nodes_list_ani:
                    for feature in self.metadata_features_list:
                        self.network_ani.nodes[node][feature] = self.metadata_dict[feature][node]

                    # Add node attribute 'SampleID' for the hover tooltip
                    self.network_ani.nodes[node]['SampleID'] = node

                # Update the color by- drop-down menus with the available metadata features
                first_feature = self.metadata_features_list[0]
                self.nodes_color_by_ani.options = self.metadata_features_list
                self.nodes_color_by_ani.value = first_feature
                self.edges_color_by_ani.options = self.metadata_features_list
                self.edges_color_by_ani.value = first_feature
                self.nodes_highlight_by_ani.options = self.metadata_features_list
                self.nodes_highlight_by_ani.value = first_feature
                # Fill the groups for the first feature
                unique_groups = sorted(list(set([str(self.network_ani.nodes[node][first_feature])
                                                 for node in self.network_ani.nodes()])))
                if 'nan' in unique_groups:
                    unique_groups.remove('nan')
                    unique_groups.append('nan')
                self.nodes_highlight_group_select_ani.options = unique_groups
                self.nodes_highlight_group_select_ani.value = unique_groups[0]

                self.nodes_colorby_ani_watcher = self.nodes_color_by_ani.param.watch(
                    self.set_not_continuous_network_ani, 'value', onlychanged=True)
                self.continuous_network_ani_watcher = self.is_continuous_network_ani.param.watch(
                    self.change_continuous_state_network_ani, 'value', onlychanged=True)
                self.colormap_ani_watcher = self.nodes_colormap_ani.param.watch(self.change_colormap_network_ani,
                                                                                'value', onlychanged=True)
                self.custom_colormap_ani_watcher = self.custom_colormap_input_ani.param.watch(
                    self.get_custom_colormap_network_ani, 'value', onlychanged=True)
                self.nodes_highlight_by_ani_watcher = self.nodes_highlight_by_ani.param.watch(
                    self.fill_feature_groups_ani, 'value', onlychanged=True)

            # No metadata
            else:
                self.use_metadata_network_ani.disabled = True

                # Add node attribute 'SampleID' for the hover tooltip
                for node in self.nodes_list_ani:
                    self.network_ani.nodes[node]['SampleID'] = node

            # Create the network plot using the selected parameters
            self.network_plot_hv_ani = pn.bind(ps.cretae_network_plot, network=self.network_ani,
                                               is_metadata=self.use_metadata_network_ani,
                                               nodes_feature=self.nodes_color_by_ani.value,
                                               is_continuous=self.is_continuous_network_ani.value,
                                               cmap=self.nodes_colormap_ani.value,
                                               custom_cmap=self.custom_colormap_input_ani.value,
                                               node_color=self.network_node_color_ani,
                                               edge_color=self.network_edge_color_ani,
                                               is_highlight_group=self.highlight_nodes_by_feature_ani,
                                               highlight_feature=self.nodes_highlight_by_ani.value,
                                               highlight_group=self.nodes_highlight_group_select_ani,
                                               is_edge_colorby=self.color_edges_by_feature_ani,
                                               edges_feature=self.edges_color_by_ani,
                                               within_edge_color=self.network_within_color_ani,
                                               between_edge_color=self.network_between_color_ani,
                                               iterations=self.network_iterations_ani, pos_dict=self.pos_dict_ani,
                                               show_labels=self.show_labels_chkbox_ani,
                                               all_or_highlighted=self.all_or_highlighted_select_ani,
                                               is_highlight_samples=self.highlight_sample_chkbox_ani,
                                               samples_to_highlight=self.highlight_sample_input_ani.value,
                                               metadata_dict=self.metadata_dict)
            self.network_pane_ani = pn.pane.HoloViews(self.network_plot_hv_ani, height=600, width=700,
                                                      sizing_mode="fixed")

            network_row = pn.Row(controls_col, pn.Spacer(width=15), self.network_pane_ani, styles={'padding': "15px"})
            self.network_card_ani.append(network_row)

    def download_network_ani(self, event):
        fformat = self.network_image_format_ani.value

        # Update the placeholder of the filenames for download with the default threshold.
        if self.use_metadata_network_ani.value:
            if self.nodes_colormap_ani.value_name == 'Define custom colormap':
                network_file = "Network_plot_ANI_" + self.ref_genome + "_" + \
                               self.network_iterations_ani.value + "_iterations_threshold_" + \
                               str(self.ani_connections_threshold) + "_colorby_" + self.nodes_color_by_ani.value + \
                               "_custom_colormap"
            else:
                network_file = "Network_plot_ANI_" + self.ref_genome + "_" + \
                               self.network_iterations.value + "_iterations_threshold_" + \
                               str(self.ani_connections_threshold) + "_colorby_" + self.nodes_color_by_ani.value + "_" + \
                               self.nodes_colormap_ani.value_name
        else:
            network_file = "Network_plot_ANI_" + self.ref_genome + "_" + \
                           self.network_iterations_ani.value + "_iterations_threshold_" + \
                           str(self.ani_connections_threshold)
        if self.show_labels_chkbox_ani.value:
            network_file += "_with_labels"

        self.save_network_plot_path_ani.placeholder = network_file
        table_file = "Network_ANI_" + self.ref_genome + "_threshold_" + str(self.ani_connections_threshold)
        self.save_network_table_path_ani.placeholder = table_file

        # Set the directory for saving
        if self.save_network_plot_path_ani.value == "":
            network_file_path = self.downloads_dir_path + self.save_network_plot_path_ani.placeholder + "." + fformat
        else:
            network_file_path = self.save_network_plot_path_ani.value

            # Add a .png suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, network_file_path, re.IGNORECASE):
                network_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(network_file_path):
                network_file_path = self.downloads_dir_path + network_file_path

        # Get the updated matplotlib plot
        self.network_plot_matplotlib_ani = pn.bind(ps.cretae_network_plot_matplotlib, network=self.network_ani,
                                                   is_metadata=self.use_metadata_network_ani,
                                                   nodes_feature=self.nodes_color_by_ani.value,
                                                   is_continuous=self.is_continuous_network_ani.value,
                                                   cmap=self.nodes_colormap_ani.value_name,
                                                   custom_cmap=self.custom_colormap_input_ani.value,
                                                   node_color=self.network_node_color_ani,
                                                   edge_color=self.network_edge_color_ani,
                                                   is_highlight_group=self.highlight_nodes_by_feature_ani,
                                                   highlight_feature=self.nodes_highlight_by_ani.value,
                                                   highlight_group=self.nodes_highlight_group_select_ani,
                                                   is_edge_colorby=self.color_edges_by_feature_ani,
                                                   edges_feature=self.edges_color_by_ani,
                                                   within_edge_color=self.network_within_color_ani,
                                                   between_edge_color=self.network_between_color_ani,
                                                   iterations=self.network_iterations_ani, pos_dict=self.pos_dict_ani,
                                                   show_labels=self.show_labels_chkbox_ani,
                                                   all_or_highlighted=self.all_or_highlighted_select_ani,
                                                   is_highlight_samples=self.highlight_sample_chkbox_ani,
                                                   samples_to_highlight=self.highlight_sample_input_ani.value,
                                                   metadata_dict=self.metadata_dict)

        # Save the network plot in the requested format using matplotlib
        self.network_plot_matplotlib_ani().savefig(network_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + network_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_network_column_ani.pop(4)
        self.download_network_column_ani.append(download_floatpanel)

    def download_network_table_ani(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_network_table_path_ani.value == "":
            network_table_path = self.downloads_dir_path + self.save_network_table_path_ani.placeholder + "." + fformat
        else:
            network_table_path = self.save_network_table_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, network_table_path, re.IGNORECASE):
                network_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(network_table_path):
                network_table_path = self.downloads_dir_path + network_table_path

        # Export just the necessary columns to csv
        self.df_for_network_ani[['Sample1', 'Sample2', 'ANI', 'weight']].to_csv(network_table_path, index=False,
                                                                                sep='\t')

        download_message = "The table is successfully saved under:\n" + network_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_network_table_column_ani.pop(2)
        self.download_network_table_column_ani.append(download_floatpanel)

    def init_positions_ani(self, event):
        self.generate_rand_positions_ani()
        self.update_network_plot_ani()

    def generate_rand_positions_ani(self):
        self.pos_dict_ani = {}

        for node in self.nodes_list_ani:
            pos_x = generate_rand_pos()
            pos_y = generate_rand_pos()
            pos_tuple = (pos_x, pos_y)
            self.pos_dict_ani[node] = pos_tuple

    def changed_threshold_select_ani(self, mean, mean_std, mean_only, event):
        #print("\nchanged_threshold_select:")
        #print("Current mean: " + str(mean))

        if self.network_threshold_select_ani.value == self.network_threshold_select_ani.options[0]:
            self.ani_connections_threshold = mean
            self.network_threshold_input_ani.disabled = True

        elif self.network_threshold_select_ani.value == self.network_threshold_select_ani.options[1]:
            # The second option is the custom threshold
            if mean_only:
                self.ani_connections_threshold = self.network_threshold_input_ani.value
                self.network_threshold_input_ani.disabled = False
            else:
                self.ani_connections_threshold = mean_std
                self.network_threshold_input_ani.disabled = True

        else:
            self.ani_connections_threshold = self.network_threshold_input_ani.value
            self.network_threshold_input_ani.disabled = False

        self.change_weight_attribute_ani()

    def changed_threshold_input_ani(self, event):
        #print("\nIn changed_threshold_input")
        self.ani_connections_threshold = self.network_threshold_input_ani.value
        self.change_weight_attribute_ani()

    def change_weight_attribute_ani(self):

        self.ani_connections_threshold = round(self.ani_connections_threshold, 3)
        print("\nchange_weight_attribute_ani:")
        print("ANI_connections_threshold = " + str(self.ani_connections_threshold))

        # Update the threshold in the default filenames for download
        network_file = "Network_plot_ANI_" + self.ref_genome + "_" + \
                       self.network_iterations_ani.value + "_iterations_threshold_" + str(self.ani_connections_threshold)
        self.save_network_plot_path_ani.placeholder = network_file
        table_file = "Network_" + self.ref_genome + "_threshold_" + str(self.ani_connections_threshold)
        self.save_network_table_path_ani.placeholder = table_file

        # Recalculate the weights
        #self.df_for_network_ani['weight'] = np.where(self.df_for_network_ani['ANI'] >= self.ani_connections_threshold,
        #                                             np.negative(np.log(1 - self.df_for_network_ani['ANI'])), 0)
        # First stage: zero the scores below the threshold
        self.df_for_network_ani['filtered_score'] = \
            self.df_for_network_ani['ANI'].where(self.df_for_network_ani['ANI'] >= self.ani_connections_threshold, 0)

        # Second stage: perform min-max normalization on the scores above the threshold
        mask = self.df_for_network_ani['filtered_score'] > 0
        min_val = self.df_for_network_ani.loc[mask, 'filtered_score'].min()
        max_val = self.df_for_network_ani.loc[mask, 'filtered_score'].max()
        self.df_for_network_ani['norm_score'] = 0.0
        # avoid division by zero if all remaining values are equal
        eps = 1e-6
        if max_val > min_val:
            self.df_for_network_ani.loc[mask, 'norm_score'] = \
                (self.df_for_network_ani.loc[mask, 'filtered_score'] - min_val) / (max_val - min_val)
            # replace exact zeros (min_val) with 10-6 (so that they will have a very low weight but not zero)
            self.df_for_network_ani.loc[mask, 'norm_score'] = \
                self.df_for_network_ani.loc[mask, 'norm_score'].clip(eps, 1 - eps)
        else:
            self.df_for_network_ani.loc[mask, 'norm_score'] = max_val

        # Third stage: apply the log function to transform the scores into weight values
        self.df_for_network_ani['weight'] = 0.0
        self.df_for_network_ani.loc[mask, 'weight'] = \
            np.negative(np.log(1 - self.df_for_network_ani.loc[mask, 'norm_score']))

        # Forth stage: calculate line width from weights (normalize values between 0.5 and 4)
        self.df_for_network_ani['width'] = 0.0
        min_width = self.df_for_network_ani.loc[mask, 'weight'].min()
        max_width = self.df_for_network_ani.loc[mask, 'weight'].max()
        new_min_width = config.min_edge_width
        new_max_width = config.max_edge_width
        self.df_for_network_ani.loc[mask, 'width'] = new_min_width + \
                                                     (self.df_for_network_ani.loc[mask, 'weight'] - min_width) * \
                                                     (new_max_width - new_min_width) / (max_width - min_width)

        #print("\nDF for network with weights and widths:")
        #print(self.df_for_network_ani)

        before = time.time()
        # Create a new network from the updated df using networkx
        self.network_ani = nx.from_pandas_edgelist(self.df_for_network_ani, source='Sample1', target='Sample2',
                                                   edge_attr=['weight', 'width'])
        self.nodes_list_ani = list(self.network_ani.nodes)

        if self.is_metadata:
            # Insert the features information as nodes attributes to the new network
            for node in self.nodes_list_ani:
                for feature in self.metadata_features_list:
                    self.network_ani.nodes[node][feature] = self.metadata_dict[feature][node]

                # Add node attribute 'SampleID' for the hover tooltip
                self.network_ani.nodes[node]['SampleID'] = node

        after = time.time()
        duration = after - before
        #print("Creating a new network took " + str(duration) + " seconds.\n")

        before = time.time()
        self.update_network_plot_ani()
        after = time.time()
        duration = after - before
        print("Updating the network plot took " + str(duration) + " seconds.\n")

    def change_highlighted_sample_ani(self, event):
        print("\nIn change_highlighted_sample_ani")
        print("Samples to highlight: " + self.highlight_sample_input_ani.value)
        self.update_network_plot_ani()

    # Update the network plot using the selected parameters and the new positions dict
    def update_network_plot_ani(self):
        #print("\nIn update_network_plot")
        self.network_plot_hv_ani = pn.bind(ps.cretae_network_plot, network=self.network_ani,
                                           is_metadata=self.use_metadata_network_ani,
                                           nodes_feature=self.nodes_color_by_ani.value,
                                           is_continuous=self.is_continuous_network_ani.value,
                                           cmap=self.nodes_colormap_ani.value,
                                           custom_cmap=self.custom_colormap_input_ani.value,
                                           node_color=self.network_node_color_ani,
                                           edge_color=self.network_edge_color_ani,
                                           is_highlight_group=self.highlight_nodes_by_feature_ani,
                                           highlight_feature=self.nodes_highlight_by_ani.value,
                                           highlight_group=self.nodes_highlight_group_select_ani,
                                           is_edge_colorby=self.color_edges_by_feature_ani,
                                           edges_feature=self.edges_color_by_ani,
                                           within_edge_color=self.network_within_color_ani,
                                           between_edge_color=self.network_between_color_ani,
                                           iterations=self.network_iterations_ani, pos_dict=self.pos_dict_ani,
                                           show_labels=self.show_labels_chkbox_ani,
                                           all_or_highlighted=self.all_or_highlighted_select_ani,
                                           is_highlight_samples=self.highlight_sample_chkbox_ani,
                                           samples_to_highlight=self.highlight_sample_input_ani.value,
                                           metadata_dict=self.metadata_dict)

        self.network_pane_ani.object = self.network_plot_hv_ani

    def create_combined_scatter_pane(self, APSS_ANI_selected_genome_df):

        self.combined_single_plots_column.clear()
        self.combined_scatter_card.clear()
        self.metadata_combined_scatter_card.clear()

        size_title = "Presenting plot using APSS from " + self.sampling_size + " subsampled regions"

        self.combined_single_plots_column.append(pn.pane.Markdown(size_title, styles={'font-size': "20px",
                                                                                      'color': config.title_purple_color,
                                                                                      'margin': "0",
                                                                                      'padding': "0"}))

        num_pairs = len(APSS_ANI_selected_genome_df)
        print("\nNumber of common sample pairs: " + str(num_pairs))

        # Check that there is at least 1 common pair. If not, print a message and don't display the plot-card.
        if num_pairs == 0:
            pairs_title = "There are no common compared sample-pairs for the selected number of subsampled regions - " \
                          "cannot display the combined plot"
            self.combined_single_plots_column.append(pn.pane.Markdown(pairs_title, styles={'font-size': "18px",
                                                                                           'color': config.title_red_color,
                                                                                           'margin': "0 0 10px 0",
                                                                                           'padding': "0"}))

        else:
            pairs_title = "Number of common compared sample pairs: " + str(num_pairs)
            self.combined_single_plots_column.append(pn.pane.Markdown(pairs_title, styles={'font-size': "18px",
                                                                                           'margin': "0 0 10px 0",
                                                                                           'padding': "0"}))

            styling_title = "Plot styling options:"
            metadata_colors_row = pn.Row(self.combined_scatter_same_color, pn.Spacer(width=3),
                                         self.combined_scatter_different_color)
            metadata_col = pn.Column(self.combined_scatter_feature_select,
                                    metadata_colors_row,
                                    styles={'padding': "10x"})
            self.metadata_combined_scatter_card.append(metadata_col)
            styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                            'color': config.title_blue_color,
                                                                            'margin': "0"}),
                                    self.combined_scatter_color,
                                    pn.Spacer(height=5),
                                    self.use_metadata_combined_scatter,
                                    self.metadata_combined_scatter_card
                                    )

            save_file_title = "Plot download options:"
            download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
            download_button.on_click(self.download_combined_scatter)

            combined_scatter_file = "ANI_vs_APSS_plot_" + self.ref_genome + "_" + self.sampling_size + "_regions"
            self.save_combined_scatter_file_path.placeholder = combined_scatter_file

            self.download_combined_scatter_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                              styles={'font-size': "15px",
                                                                      'font-weight': "bold",
                                                                      'color': config.title_blue_color,
                                                                      'margin': "0"}),
                                                              self.combined_scatter_image_format,
                                                              self.save_combined_scatter_file_path,
                                                              download_button, pn.pane.Markdown())

            # Use metadata in plot
            if self.is_metadata:
                # Update the color nodes by- drop-down menu with the available metadata features
                self.combined_scatter_feature_select.options = self.metadata_features_list
                self.combined_scatter_feature_select.value = self.metadata_features_list[0]
            # No metadata
            else:
                self.use_metadata_combined_scatter.disabled = True

            # Create the scatter plot
            self.combined_scatter_plot = pn.bind(self.create_combined_scatter_plot,
                                                 combined_df=APSS_ANI_selected_genome_df,
                                                 color=self.combined_scatter_color,
                                                 use_metadata=self.use_metadata_combined_scatter,
                                                 feature=self.combined_scatter_feature_select,
                                                 same_color=self.combined_scatter_same_color,
                                                 different_color=self.combined_scatter_different_color)

            combined_scatter_pane = pn.pane.Matplotlib(self.combined_scatter_plot, height=550, dpi=300, tight=True,
                                                       format='png')

            combined_scatter_table_file = "Data_for_ANI_vs_APSS_plot_" + self.ref_genome + "_" + self.sampling_size + \
                                          "_regions"
            self.save_combined_scatter_table_path.placeholder = combined_scatter_table_file
            download_table_button = pn.widgets.Button(name='Download data table in tsv format', button_type='primary')
            download_table_button.on_click(self.download_combined_scatter_table)

            self.download_combined_scatter_table_column = pn.Column(self.save_combined_scatter_table_path,
                                                                    download_table_button,
                                                                    pn.pane.Markdown())

            controls_col = pn.Column(styling_col, pn.Spacer(height=25), self.download_combined_scatter_column,
                                     self.download_combined_scatter_table_column)

            combined_scatter_row = pn.Row(controls_col, pn.Spacer(width=120), combined_scatter_pane,
                                          styles={'padding': "15px"})
            self.combined_scatter_card.append(combined_scatter_row)

            self.combined_single_plots_column.append(self.combined_scatter_card)

    def download_combined_scatter(self, event):
        fformat = self.combined_scatter_image_format.value

        # Set the directory for saving
        if self.save_combined_scatter_file_path.value == "":
            combined_scatter_file_path = self.downloads_dir_path + self.save_combined_scatter_file_path.placeholder \
                                         + "." + fformat
        else:
            combined_scatter_file_path = self.save_combined_scatter_file_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, combined_scatter_file_path, re.IGNORECASE):
                combined_scatter_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(combined_scatter_file_path):
                combined_scatter_file_path = self.downloads_dir_path + combined_scatter_file_path

        self.combined_scatter_plot().savefig(combined_scatter_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + combined_scatter_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_combined_scatter_column.pop(4)
        self.download_combined_scatter_column.append(download_floatpanel)

    def download_combined_scatter_table(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_combined_scatter_table_path.value == "":
            combined_scatter_table_path = self.downloads_dir_path + self.save_combined_scatter_table_path.placeholder \
                                          + "." + fformat
        else:
            combined_scatter_table_path = self.save_combined_scatter_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, combined_scatter_table_path, re.IGNORECASE):
                combined_scatter_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(combined_scatter_table_path):
                combined_scatter_table_path = self.downloads_dir_path + combined_scatter_table_path

        self.df_for_combined_scatter.to_csv(combined_scatter_table_path, sep='\t', index=False)

        download_message = "The table is successfully saved under:\n" + combined_scatter_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_combined_scatter_table_column.pop(2)
        self.download_combined_scatter_table_column.append(download_floatpanel)

    def category_by_feature_scatter(self, row, feature, same_color, diff_color):
        if is_missing_value(self.metadata_dict[feature][row['Sample1']]) \
                or is_missing_value(self.metadata_dict[feature][row['Sample2']]):
            return 'gray'
        else:
            if self.metadata_dict[feature][row['Sample1']] == self.metadata_dict[feature][row['Sample2']]:
                return same_color
            else:
                return diff_color

    def create_combined_scatter_plot(self, combined_df, color, use_metadata, feature, same_color, different_color):

        self.df_for_combined_scatter = combined_df.copy()

        fig, ax1 = plt.subplots(figsize=(8, 7))

        # Use metadata to separate plot to same/different feature
        if use_metadata:
            self.df_for_combined_scatter['Color'] = self.df_for_combined_scatter.apply(
                lambda row: self.category_by_feature_scatter(row, feature, same_color, different_color), axis=1)
            #print(self.df_for_combined_scatter)

            ax1.scatter(self.df_for_combined_scatter['APSS'], self.df_for_combined_scatter['ANI'],
                        c=self.df_for_combined_scatter['Color'], linewidths=0.1, edgecolors="gray")

            # Define the color-to-category (same/different) mapping
            color_to_category = {
                same_color: 'Same ' + feature,
                different_color: 'Different ' + feature,
            }
            # The selected feature contains 'Unknown' values - add this category to the legend
            if (self.df_for_combined_scatter['Color'] == 'gray').any():
                color_to_category['gray'] = 'Unknown'

            # Create legend handles
            handles = [
                mpatches.Patch(color=color, label=label)
                for color, label in color_to_category.items()
            ]

            # Add legend to the plot
            ax1.legend(handles=handles, loc='lower left', bbox_to_anchor=(0, 1), fontsize=12)

        # Do not use metadata in plot - show all the comparisons together
        else:
            ax1.scatter(self.df_for_combined_scatter['APSS'], self.df_for_combined_scatter['ANI'],
                        c=color, linewidths=0.1, edgecolors="gray")

        plt.xlabel("APSS", labelpad=5, fontsize=12)
        plt.ylabel("ANI", labelpad=5, fontsize=12)

        # Calculate the Spearman correlation for the plot
        corr, p_value = spearmanr(self.df_for_combined_scatter['APSS'], self.df_for_combined_scatter['ANI'])

        # Add to plot (adjust x, y text position as needed)
        ax1.text(0.05, 0.95, f"r = {corr:.2f}", transform=ax1.transAxes, fontsize=12, verticalalignment='top')

        plt.close(fig)

        return fig

    def changed_contig_sorting(self, sorting_select, contig_select, event):
        sorting_method = sorting_select.value
        print("\nChanged contigs sorting method. Sort by: " + sorting_method)

        # Change the order of the contigs selection and present the first one in the new order
        if sorting_method == "Contig name":
            contig_select.options = self.contigs_list_by_name
            contig_select.value = self.contigs_list_by_name[0]
        else:
            contig_select.options = self.contigs_list_by_length
            contig_select.value = self.contigs_list_by_length[0]

    def changed_contig(self, contig, event):
        self.contig_name = contig.value
        print("\nChanged_contig, contig name: " + self.contig_name)

        # Unwatch all contig-specific widgets (if it's not the first time that the contig has changed)
        if self.visited_synteny_per_pos_tab:
            print("changed_contig: remove contig-related watchers")
            self.avg_plot_chkbox.param.unwatch(self.avg_plot_chkbox_watcher)
            self.avg_plot_color.param.unwatch(self.avg_plot_color_watcher)
            self.coverage_plot_chkbox.param.unwatch(self.coverage_plot_chkbox_watcher)
            self.coverage_plot_color.param.unwatch(self.coverage_plot_color_watcher)
            self.hypervar_chkbox.param.unwatch(self.hypervar_chkbox_watcher)
            self.hypervar_color.param.unwatch(self.hypervar_color_watcher)
            self.hypercons_chkbox.param.unwatch(self.hypercons_chkbox_watcher)
            self.hypercons_color.param.unwatch(self.hypercons_color_watcher)
            self.alpha_slider.param.unwatch(self.alpha_slider_watcher)
            self.show_annotations_chkbox.param.unwatch(self.show_annotations_watcher)
            if self.is_metadata:
                self.synteny_per_pos_feature_select.param.unwatch(self.synteny_per_pos_feature_select_watcher)

        self.coverage_plot = ""
        self.line_avg_plot = ""
        self.hypervar_bars = ""
        self.hypercons_bars = ""
        self.filter_plot_by_metadata = 0
        self.show_annotations = False

        self.create_selected_contig_column()

    def create_selected_contig_column(self):

        self.score_per_pos_contig = self.score_per_region_selected_genome_df[
            self.score_per_region_selected_genome_df['Contig_name'] == self.contig_name].copy()

        self.selected_contig_column.clear()
        self.annotation_upload_col.clear()
        self.show_annotations_col.clear()
        self.show_annotations_row.clear()

        contig_name_title = "Contig name: " + self.contig_name
        self.selected_contig_column.append(pn.pane.Markdown(contig_name_title,
                                                            styles={'font-size': "17px",
                                                                    'color': config.title_purple_color,
                                                                    'margin': "5px 5px 0px 5px",
                                                                    'padding-bottom': "0px"}))

        # Find contig length by the last position
        self.score_per_pos_contig['Position'] = self.score_per_pos_contig['Position'].astype(int)
        self.score_per_pos_contig = self.score_per_pos_contig.sort_values('Position')
        print("\nLast position: " + str(self.score_per_pos_contig.iloc[-1]['Position']))
        self.contig_length = self.score_per_pos_contig.iloc[-1]['Position'] + config.region_length
        contig_length_title = "Contig length: " + str(self.contig_length) + " bp"
        self.selected_contig_column.append(pn.pane.Markdown(contig_length_title,
                                                            styles={'font-size': "16px", 'margin': "0px 3px 3px 5px",
                                                                    'padding-top': "0px"}))

        ###################################
        # Create the customization column
        # All the binded widgets have to be defined again to prevent wrong calls to create_synteny_per_pos_plot
        # with previous contigs

        # Set the length range
        start_pos = '0'
        end_pos = str(self.contig_length)
        self.start_pos_input.placeholder = start_pos
        self.start_pos_input.value = start_pos
        self.end_pos_input.placeholder = end_pos
        self.end_pos_input.value = end_pos

        reset_range_button = pn.widgets.Button(name='Reset range', button_type='primary',
                                               styles={'margin-top': "22px"})
        reset_range_button.on_click(self.reset_range)
        change_range_button = pn.widgets.Button(name='Set new range', button_type='primary',
                                                styles={'margin-top': "22px"})
        change_range_button.on_click(self.change_range)

        pos_range_cust_row = pn.Row(pn.pane.Markdown("Set contig length range:",
                                                     styles={'font-size': "14px", 'margin': "10px 5px 5px 5px"}),
                                    pn.Spacer(width=10), self.start_pos_input, pn.Spacer(width=5), self.end_pos_input,
                                    pn.Spacer(width=5), change_range_button, pn.Spacer(width=5), reset_range_button,
                                    styles={'margin-left': "5px"})

        avg_plot_chkbox_col = pn.Column(pn.Spacer(height=25), self.avg_plot_chkbox)
        coverage_plot_chkbox_col = pn.Column(pn.Spacer(height=25), self.coverage_plot_chkbox)

        coverage_plot_cust_row = pn.Row(avg_plot_chkbox_col, pn.Spacer(width=5), self.avg_plot_color,
                                        pn.Spacer(width=30),
                                        coverage_plot_chkbox_col, pn.Spacer(width=5), self.coverage_plot_color,
                                        styles={'margin-left': "5px"})

        hypervar_chkbox_col = pn.Column(pn.Spacer(height=25), self.hypervar_chkbox)
        hypercons_chkbox_col = pn.Column(pn.Spacer(height=25), self.hypercons_chkbox)

        hypercons_cust_row = pn.Row(hypervar_chkbox_col, pn.Spacer(width=5), self.hypervar_color,
                                    pn.Spacer(width=30),
                                    hypercons_chkbox_col, pn.Spacer(width=5), self.hypercons_color,
                                    pn.Spacer(width=30),
                                    self.alpha_slider,
                                    styles={'margin-left': "5px"})

        # There is metadata
        if self.is_metadata:

            self.filter_by_metadata_card.clear()

            # Update the synteny_per_pos_feature_select drop-down menu with the available metadata features
            feature = self.metadata_features_list[0]
            self.synteny_per_pos_feature_select.options = self.metadata_features_list
            self.synteny_per_pos_feature_select.value = feature

            # Fill the groups for the first feature
            unique_groups = sorted(self.groups_per_feature_dict[feature], key=str)
            if 'nan' in unique_groups:
                unique_groups.remove('nan')
                unique_groups.append('nan')
            self.synteny_per_pos_groups_select.options = unique_groups

            # Define a watcher for the feature-selection change
            self.synteny_per_pos_feature_select_watcher = self.synteny_per_pos_feature_select.param.watch(
                self.fill_groups_for_multiselect, 'value', onlychanged=True)

            buttons_row = pn.Row(self.filter_synteny_per_pos_button, pn.Spacer(width=10),
                                 self.reset_filter_synteny_per_pos_button)
            buttons_col = pn.Column(pn.Spacer(height=20), buttons_row)
            filter_metadata_row = pn.Row(self.synteny_per_pos_feature_select, pn.Spacer(width=10),
                                         self.synteny_per_pos_groups_select, pn.Spacer(width=10),
                                         buttons_col,
                                         styles={'padding': "10px"})
            self.filter_by_metadata_card.append(filter_metadata_row)

        # No metadata
        else:
            self.filter_by_metadata_card.collapsible = False
            self.filter_by_metadata_card.hide_header = True

        self.add_annotation_card.clear()

        # No annotation file for the current genome has already been uploaded - add the FileInput widget
        if self.annotation_per_ref_genome_dict[self.ref_genome] is None:
            annotation_upload_title = "Upload an annotation file for the reference genome (in gff format):"
            self.upload_annotation_file.disabled = False
            self.submit_annotation_file_button.disabled = False
            self.show_annotations_chkbox.disabled = True

            self.annotation_upload_col.append(pn.pane.Markdown(annotation_upload_title,
                                                               styles={'font-size': "14px",
                                                                       'margin-bottom': "0",
                                                                       'margin-top': "0",
                                                                       'padding-bottom': "0"}))

        # An annotation file for the current genome has already been uploaded - write a message
        else:
            annotation_upload_title = "Found annotation file for the current reference genome: " + self.gff_filename
            self.upload_annotation_file.disabled = True
            self.submit_annotation_file_button.disabled = True
            self.show_annotations_chkbox.value = False
            self.show_annotations_chkbox.disabled = False

            self.annotation_upload_col.append(pn.pane.Markdown(annotation_upload_title,
                                                               styles={'font-size': "14px",
                                                                       'color': config.title_red_color,
                                                                       'margin-bottom': "0",
                                                                       'margin-top': "0",
                                                                       'padding-bottom': "0"}))

        annotation_upload_row = pn.Row(self.upload_annotation_file, self.submit_annotation_file_button)
        self.annotation_upload_col.append(annotation_upload_row)

        self.show_annotations_col.append(pn.Spacer(height=50))
        self.show_annotations_row = pn.Row(self.show_annotations_chkbox, self.annotation_help_button,
                                           pn.Spacer(width=1))
        self.show_annotations_col.append(self.show_annotations_row)
        add_annotation_row = pn.Row(self.annotation_upload_col, pn.Spacer(width=45), self.show_annotations_col,
                                    styles={'padding': "10px"})

        self.add_annotation_card.append(add_annotation_row)

        styling_title = "Customization options:"
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "10px 5px 5px 5px"}),
                                pos_range_cust_row,
                                coverage_plot_cust_row,
                                hypercons_cust_row,
                                self.filter_by_metadata_card,
                                self.add_annotation_card)

        # Define watchers for visual widgets
        self.avg_plot_chkbox_watcher = self.avg_plot_chkbox.param.watch(self.show_hide_avg_plot, 'value',
                                                                        onlychanged=True)
        self.avg_plot_color_watcher = self.avg_plot_color.param.watch(self.change_avg_plot_color, 'value',
                                                                      onlychanged=True)
        self.coverage_plot_chkbox_watcher = self.coverage_plot_chkbox.param.watch(self.show_hide_coverage_plot, 'value',
                                                                                  onlychanged=True)
        self.coverage_plot_color_watcher = self.coverage_plot_color.param.watch(self.change_coverage_plot_color,
                                                                                'value', onlychanged=True)
        self.hypervar_chkbox_watcher = self.hypervar_chkbox.param.watch(self.show_hide_hypervar_plot, 'value',
                                                                        onlychanged=True)
        self.hypervar_color_watcher = self.hypervar_color.param.watch(self.change_hypervar_color, 'value',
                                                                      onlychanged=True)
        self.hypercons_chkbox_watcher = self.hypercons_chkbox.param.watch(self.show_hide_hypercons_plot,
                                                                          'value', onlychanged=True)
        self.hypercons_color_watcher = self.hypercons_color.param.watch(self.change_hypercons_color,
                                                                        'value', onlychanged=True)
        self.alpha_slider_watcher = self.alpha_slider.param.watch(self.change_alpha, 'value', onlychanged=True)
        self.show_annotations_watcher = self.show_annotations_chkbox.param.watch(self.show_hide_annotations_plot,
                                                                                 'value', onlychanged=True)

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_synteny_per_pos_plot)

        synteny_per_pos_file = "Synteny_per_position_plot_" + self.ref_genome + "_" + self.contig_name
        self.save_synteny_per_pos_plot_path.placeholder = synteny_per_pos_file

        self.download_synteny_per_pos_plot_column = pn.Column(self.synteny_per_pos_image_format,
                                                              self.save_synteny_per_pos_plot_path,
                                                              pn.Spacer(height=10), download_button, pn.pane.Markdown(),
                                                              align='end')

        synteny_per_pos_table = "Data_for_synteny_per_position_plot_" + self.ref_genome + "_" + self.contig_name
        self.save_synteny_per_pos_table_path.placeholder = synteny_per_pos_table

        download_table_button = pn.widgets.Button(name='Download underlying data in tsv format', button_type='primary')
        download_table_button.on_click(self.download_synteny_per_pos_table)

        self.download_synteny_per_pos_table_column = pn.Column(self.save_synteny_per_pos_table_path,
                                                               pn.Spacer(height=10),
                                                               download_table_button, pn.pane.Markdown(), align='end')

        download_files_row = pn.Row(self.download_synteny_per_pos_plot_column, pn.Spacer(width=100),
                                    self.download_synteny_per_pos_table_column)

        download_synteny_per_pos_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                    styles={'font-size': "15px", 'font-weight': "bold",
                                                            'color': config.title_blue_color,
                                                            'margin': "10px 5px 5px 5px"}), download_files_row)

        # Prepare the data structures necessary for the synteny_per_pos plots of a specific contig
        # Calculate the average synteny scores for each position
        self.avg_score_per_pos_contig = self.score_per_pos_contig[['Contig_name', 'Position', 'Synteny_score']]. \
            sort_values(['Position']).groupby('Position') \
            .agg(Count=('Synteny_score', 'size'), Avg_synteny_score=('Synteny_score', 'mean')) \
            .reset_index()
        #print("\nAverage + count df:")
        #print(avg_score_per_pos_contig)

        # Fill the missing positions with score=0 (default jump=5000)
        self.avg_score_per_pos_contig = self.avg_score_per_pos_contig. \
            merge(how='right', on='Position',
                  right=pd.DataFrame(
                      {'Position': np.arange(self.avg_score_per_pos_contig.iloc[0]['Position'],
                                             self.avg_score_per_pos_contig.iloc[-1]['Position'] + 2 * config.region_length,
                                             config.region_length)
                       })).sort_values(by='Position').reset_index(). \
            drop(['index'], axis=1)
        self.avg_score_per_pos_contig['Position'] = self.avg_score_per_pos_contig['Position'].astype(int)

        #print("\nAfter filling missing positions:")
        #print(avg_score_per_pos_contig)

        hypervar_threshold = self.bottom_percentile
        hypercons_threshold = self.top_percentile

        self.avg_score_per_pos_contig['Hypervariable'] = np.where(
            (self.avg_score_per_pos_contig['Avg_synteny_score'] <= hypervar_threshold) &
            (self.avg_score_per_pos_contig['Count'] >= self.bottom_percentile_counts), 1, 0)
        self.avg_score_per_pos_contig['Hyperconserved'] = np.where(
            (self.avg_score_per_pos_contig['Avg_synteny_score'] >= hypercons_threshold) &
            (self.avg_score_per_pos_contig['Count'] >= self.median_counts), 1, 0)

        # Create the synteny_per_pos plot for the selected contig of the Ref-genome with the customized parameters
        self.synteny_per_pos_plot = self.create_synteny_per_pos_plot()

        self.synteny_per_pos_pane = pn.pane.Matplotlib(self.synteny_per_pos_plot, dpi=300, tight=True, format='png',
                                                       width=1105)

        synteny_per_pos_plot_row = pn.Row(styles={'padding': "0", 'margin': "0"})
        synteny_per_pos_plot_row.append(self.synteny_per_pos_pane)

        self.selected_contig_column.append(synteny_per_pos_plot_row)
        self.selected_contig_column.append(pn.Spacer(width=20))
        self.selected_contig_column.append(styling_col)
        self.selected_contig_column.append(download_synteny_per_pos_column)

    def read_annotation_file(self, event):
        gff_file = io.BytesIO(self.upload_annotation_file.value)
        gff_text = gff_file.read().decode("utf-8")
        self.gff_filename = self.upload_annotation_file.filename
        print("\nUploaded annotation file: " + self.gff_filename)

        # Read the gff file and get a dict of features-per-contig for the current ref-genome
        self.annotation_per_ref_genome_dict[self.ref_genome] = load_genes_from_gff(gff_text)

        annotation_upload_title = "Uploaded annotation file: " + self.gff_filename
        self.upload_annotation_file.disabled = True
        self.submit_annotation_file_button.disabled = True
        self.show_annotations_chkbox.disabled = False

        self.annotation_upload_col.clear()
        self.annotation_upload_col.append(pn.pane.Markdown(annotation_upload_title,
                                                           styles={'font-size': "14px",
                                                                   'color': config.title_red_color,
                                                                   'margin-bottom': "0",
                                                                   'margin-top': "0",
                                                                   'padding-bottom': "0"}))
        annotation_upload_row = pn.Row(self.upload_annotation_file, self.submit_annotation_file_button)
        self.annotation_upload_col.append(annotation_upload_row)

    def show_annotations_help_float_panel(self, event):
        message = "Plot annotated genes at the bottom of the main plot.\n" \
                  "The genes plot can only be displayed if the contig's length does not exceed 100,000 bp."
        markdown = pn.pane.Markdown(message, styles={'font-size': "14px", 'color': config.title_red_color})
        floatpanel = pn.layout.FloatPanel(markdown, name='genes plot help', margin=10, width=400)
        self.show_annotations_row[2] = floatpanel

    def fill_groups_for_multiselect(self, event):
        feature = self.synteny_per_pos_feature_select.value
        unique_groups = sorted(self.groups_per_feature_dict[feature], key=str)

        # Move the 'nan' group (if any) to the end of the list
        if 'nan' in unique_groups:
            unique_groups.remove('nan')
            unique_groups.append('nan')
        #print(unique_groups)

        self.synteny_per_pos_groups_select.options = unique_groups

        # Reset filteration if needed
        if self.filter_plot_by_metadata:
            self.filter_plot_by_metadata = 0
            self.update_synteny_per_pos_plot()

    def create_synteny_per_pos_plot(self):
        before = time.time()
        #print("\ncreate_synteny_per_pos_plot:\nContig name: " + self.contig_name)

        # Set the requested positions range
        start_pos = self.start_pos_input.value
        end_pos = self.end_pos_input.value
        #print("Start position: " + start_pos)
        #print("End position: " + end_pos)

        # The user requested to filter the data by a metadata feature - use the filtered tables
        if self.filter_plot_by_metadata:
            # Consider only the positions within the requested range
            score_per_pos_contig = self.score_per_pos_contig_filtered[
                self.score_per_pos_contig_filtered['Position'] >= int(start_pos)]
            score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] < int(end_pos)]
            avg_score_per_pos_contig = self.avg_score_per_pos_contig_filtered[
                self.avg_score_per_pos_contig_filtered['Position'] >= int(start_pos)]
            avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

        # No filtering by metadata - use the full tables
        else:
            # Consider only the positions within the requested range
            score_per_pos_contig = self.score_per_pos_contig[self.score_per_pos_contig['Position'] >= int(start_pos)]
            score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] < int(end_pos)]
            avg_score_per_pos_contig = self.avg_score_per_pos_contig[
                self.avg_score_per_pos_contig['Position'] >= int(start_pos)]
            avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

        #print("\nscore_per_pos_contig table:")
        #print(score_per_pos_contig)
        #print("\navg_score_per_pos_contig table:")
        #print(avg_score_per_pos_contig)

        # Prepare data for plotting the avg scores as lines
        avg_score_per_pos_contig_end_pos = avg_score_per_pos_contig.copy()
        avg_score_per_pos_contig_end_pos['Position'] = avg_score_per_pos_contig['Position'] + config.region_length - 50
        # print(avg_score_per_pos_contig_end_pos)

        avg_score_per_pos_contig_for_line_plot = pd.concat([avg_score_per_pos_contig, avg_score_per_pos_contig_end_pos],
                                                           ignore_index=True).sort_values(by='Position')

        #print("\nAverage data for line plot after concatenating:")
        #print(avg_score_per_pos_contig_for_line_plot)

        pos_array = np.full((2, len(score_per_pos_contig.index)), 0)
        pos_array[1, :] = config.region_length - 50

        avg_pos_array = np.full((2, len(avg_score_per_pos_contig.index)), 0)
        avg_pos_array[1, :] = config.region_length

        # Create a figure with two subplots, one for the synteny information and one for the annotations
        if self.show_annotations:
            fig, (self.ax_for_synteny_per_pos_plot, self.ax_annotations) = plt.subplots(2, 1, figsize=(12, 8),
                                                                                        sharex=True,
                                                                                        gridspec_kw={'height_ratios': [3, 2]})
        # Create a figure with one subplot for the synteny only
        else:
            fig, self.ax_for_synteny_per_pos_plot = plt.subplots(figsize=(12, 5))

        self.coverage_plot = self.ax_for_synteny_per_pos_plot.errorbar(score_per_pos_contig['Position'],
                                                                       score_per_pos_contig['Synteny_score'],
                                                                       xerr=pos_array,
                                                                       color=self.coverage_plot_color.value, fmt='none',
                                                                       elinewidth=0.7, zorder=1, label='Synteny scores')
        # Hide the coverage plot if the checkbox is not checked
        if not self.coverage_plot_chkbox.value:
            for barline in self.coverage_plot[2]:
                barline.set_visible(False)
            self.coverage_plot.set_label('_Synteny scores')  # Hide the label in the legend

        self.line_avg_plot = self.ax_for_synteny_per_pos_plot.plot(avg_score_per_pos_contig_for_line_plot['Position'],
                                                                   avg_score_per_pos_contig_for_line_plot['Avg_synteny_score'],
                                                                   color=self.avg_plot_color.value, zorder=2,
                                                                   label='Average synteny scores')
        # Hide the average plot if the checkbox is not checked
        if not self.avg_plot_chkbox.value:
            self.line_avg_plot[0].set_visible(False)
            self.line_avg_plot[0].set_label('_Average synteny scores')

        height = 1
        bottom_val = 0
        min_score = self.score_per_pos_contig['Synteny_score'].min()
        if min_score < 0:
            height += abs(min_score) + 0.05
            bottom_val = min_score - 0.05
        #print("\nMin score = " + str(min_score))
        #print("Height = " + str(height))

        avg_score_per_pos_contig['Hypervariable'] = np.where(avg_score_per_pos_contig['Hypervariable'] == 0, 0, height)
        avg_score_per_pos_contig['Hyperconserved'] = np.where(avg_score_per_pos_contig['Hyperconserved'] == 0, 0, height)

        self.hypervar_bars = self.ax_for_synteny_per_pos_plot.bar(avg_score_per_pos_contig['Position'],
                                     avg_score_per_pos_contig['Hypervariable'], align='edge',
                                     width=config.region_length, bottom=bottom_val, color=self.hypervar_color.value,
                                     linewidth=0, alpha=self.alpha_slider.value, label='Hypervariable regions')

        # Hide the hyper variability plot if the checkbox is not checked
        if not self.hypervar_chkbox.value:
            for bar in self.hypervar_bars:
                bar.set_visible(False)
            self.hypervar_bars.set_label('_Hypervariable regions')

        self.hypercons_bars = self.ax_for_synteny_per_pos_plot.bar(avg_score_per_pos_contig['Position'],
                                      avg_score_per_pos_contig['Hyperconserved'], align='edge',
                                      width=config.region_length, bottom=bottom_val, color=self.hypercons_color.value,
                                      linewidth=0, alpha=self.alpha_slider.value, label='Hyperconserved regions')

        # Hide the hyper conservation plot if the checkbox is not checked
        if not self.hypercons_chkbox.value:
            for bar in self.hypercons_bars:
                bar.set_visible(False)
            self.hypercons_bars.set_label('_Hyperconserved regions')

        # Add a dummy row to the end of the df in order to include the end-position in the X-axis
        avg_score_per_pos_contig = avg_score_per_pos_contig.reset_index()
        avg_score_per_pos_contig.loc[len(avg_score_per_pos_contig)] = [len(avg_score_per_pos_contig), int(end_pos), 0,
                                                                       0, 0, 0]
        #print("\nFinal AVG score per position table:")
        #print(avg_score_per_pos_contig)

        # Set the X-ticks and labels of the main plot
        positions = avg_score_per_pos_contig['Position']
        self.ax_for_synteny_per_pos_plot.set_xlim(int(start_pos), int(end_pos))
        self.ax_for_synteny_per_pos_plot.set_xticks(positions)
        self.ax_for_synteny_per_pos_plot.tick_params(axis='x', labelsize=6, rotation=90)
        self.ax_for_synteny_per_pos_plot.set_xlabel("Position in reference genome/contig", labelpad=8, fontsize=10)

        # Set the Y-label of the synteny plot
        self.ax_for_synteny_per_pos_plot.set_ylabel("Synteny score", labelpad=8, fontsize=10)

        # Set the legend
        self.ax_for_synteny_per_pos_plot.legend(fontsize='small', loc=(0, 1.02))

        # Customize the axes in order to display the synteny + annotations plot
        if self.show_annotations:

            # Filter to features that overlap the visible x-range (in absolute coordinates)
            features_in_range = set_features_in_range(
                self.annotation_per_ref_genome_dict[self.ref_genome][self.contig_name], int(start_pos), int(end_pos))

            # Draw the annotations plot
            seq_length = int(end_pos) - int(start_pos)
            record = GraphicRecord(
                sequence_length=seq_length,
                features=features_in_range,
                first_index=int(start_pos)
            )
            record.plot(ax=self.ax_annotations, with_ruler=True)

            # Set the X-ticks of both plots according to the positions
            self.ax_annotations.set_xticks(positions)
            self.ax_annotations.tick_params(axis='x', labelsize=7, labelcolor='black')

            self.ax_for_synteny_per_pos_plot.tick_params(axis='x', labelsize=7, reset=True, top=False)
            self.ax_for_synteny_per_pos_plot.set_xlabel("")

            # ---- FIX Y-AXIS VISIBILITY ----
            # Re-enable the left spine (DFV turns it off)
            self.ax_annotations.spines['left'].set_visible(True)

            # Make sure Matplotlib knows it should show the y-axis
            self.ax_annotations.yaxis.set_visible(True)

            # Enable invisible ticks (even invisible ticks help trigger label placement)
            self.ax_annotations.tick_params(axis='y', left=False, labelleft=False)

            # Set the Y-label
            self.ax_annotations.set_ylabel("Annotated genes", labelpad=36, fontsize=10)

            # Add global X-label for the entire figure
            fig.supxlabel("Position in reference genome/contig", y=0.04, fontsize=10)

        after = time.time()
        duration = after - before
        print("Create/update the synteny_per_pos plot took " + str(duration) + " seconds")

        return fig

    def update_synteny_per_pos_plot(self):
        self.synteny_per_pos_plot = self.create_synteny_per_pos_plot()
        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_range(self, event):

        # Check if show_annotations is checked
        if self.show_annotations_chkbox.value:
            length_range = int(self.end_pos_input.value) - int(self.start_pos_input.value)

            # The requested length range is too big - Print a message that the length exceeds the allowed maximum
            if length_range > config.max_range_for_annotation:
                # Print a message that the length exceeds the allowed maximum
                print("\nchange_range: length range is too big")
                message = "Cannot plot annotated genes - set contig length to max. 100,000 bp"
                self.show_annotations_col[0] = pn.pane.Markdown(message,
                                                                styles={'font-size': "14px",
                                                                        'color': config.title_red_color,
                                                                        'margin-bottom': "0",
                                                                        'margin-top': "0",
                                                                        'padding-bottom': "0"})
                self.show_annotations = False

            # The requested length range is suitable for showing the annotation plot
            else:
                self.show_annotations_col[0] = pn.Spacer(height=50)
                self.show_annotations = True

        # show_annotations is not checked
        else:
            self.show_annotations = False

        self.update_synteny_per_pos_plot()

    def reset_range(self, event):
        start_pos = '0'
        end_pos = str(self.contig_length)

        self.start_pos_input.placeholder = start_pos
        self.start_pos_input.value = start_pos
        self.end_pos_input.placeholder = end_pos
        self.end_pos_input.value = end_pos
        print("\nIn reset_range")
        print("Start=" + start_pos + ", end=" + end_pos)

        # Check if show_annotations is checked
        if self.show_annotations_chkbox.value:
            length_range = int(self.end_pos_input.value) - int(self.start_pos_input.value)

            # The requested length range is too big - Print a message that the length exceeds the allowed maximum
            if length_range > config.max_range_for_annotation:
                print("\nreset_range: length range is too big")
                message = "Cannot plot annotated genes - set contig length to max. 100,000 bp"
                self.show_annotations_col[0] = pn.pane.Markdown(message,
                                                                styles={'font-size': "14px",
                                                                        'color': config.title_red_color,
                                                                        'margin-bottom': "0",
                                                                        'margin-top': "0",
                                                                        'padding-bottom': "0"})
                self.show_annotations = False

            # The requested length range is suitable for showing the annotation plot
            else:
                self.show_annotations_col[0] = pn.Spacer(height=50)
                self.show_annotations = True

        # show_annotations is not checked
        else:
            self.show_annotations = False

        self.update_synteny_per_pos_plot()

    def update_legend(self):
        self.ax_for_synteny_per_pos_plot.get_legend().remove()
        self.ax_for_synteny_per_pos_plot.legend(fontsize='small', loc=(0, 1.02))

    def show_hide_avg_plot(self, event):
        if self.avg_plot_chkbox.value:
            self.line_avg_plot[0].set_visible(True)
            self.line_avg_plot[0].set_label('Average synteny scores')
        else:
            self.line_avg_plot[0].set_visible(False)
            self.line_avg_plot[0].set_label('_Average synteny scores')

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_avg_plot_color(self, event):
        color = self.avg_plot_color.value
        print("\nIn change_avg_plot_color. New color: " + color)
        self.line_avg_plot[0].set_color(color)

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def show_hide_coverage_plot(self, event):
        if self.coverage_plot_chkbox.value:
            for barline in self.coverage_plot[2]:
                barline.set_visible(True)
            self.coverage_plot.set_label('Synteny scores')
        else:
            for barline in self.coverage_plot[2]:
                barline.set_visible(False)
            self.coverage_plot.set_label('_Synteny scores')  # Hide the label in the legend

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_coverage_plot_color(self, event):
        color = self.coverage_plot_color.value

        for barline in self.coverage_plot[2]:  # 2 = error bar lines
            barline.set_color(color)

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def show_hide_hypervar_plot(self, event):
        if self.hypervar_chkbox.value:
            for bar in self.hypervar_bars:
                bar.set_visible(True)
            self.hypervar_bars.set_label('Hypervariable regions')
        else:
            for bar in self.hypervar_bars:
                bar.set_visible(False)
            self.hypervar_bars.set_label('_Hypervariable regions')

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_hypervar_color(self, event):
        color = self.hypervar_color.value

        for bar in self.hypervar_bars:
            bar.set_color(color)

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def show_hide_hypercons_plot(self, event):
        if self.hypercons_chkbox.value:
            for bar in self.hypercons_bars:
                bar.set_visible(True)
            self.hypercons_bars.set_label('Hyperconserved regions')
        else:
            for bar in self.hypercons_bars:
                bar.set_visible(False)
            self.hypercons_bars.set_label('_Hyperconserved regions')

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_hypercons_color(self, event):
        color = self.hypercons_color.value

        for bar in self.hypercons_bars:
            bar.set_color(color)

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def change_alpha(self, event):
        alpha = self.alpha_slider.value

        for bar in self.hypervar_bars:
            bar.set_alpha(alpha)
        for bar in self.hypercons_bars:
            bar.set_alpha(alpha)

        self.update_legend()

        self.synteny_per_pos_pane.object = self.synteny_per_pos_plot

    def show_hide_annotations_plot(self, event):
        # Show annotations is checked
        if self.show_annotations_chkbox.value:
            # Check the current contig length range. I f it's longer than 50000, uncheck the box and print a message
            length_range = int(self.end_pos_input.value) - int(self.start_pos_input.value)
            if length_range > config.max_range_for_annotation:
                print("\nshow_hide_annotations_plot: length range is too big")
                message = "Cannot plot annotated genes - set contig length to max. 100,000 bp"
                self.show_annotations_col[0] = pn.pane.Markdown(message,
                                                                styles={'font-size': "14px",
                                                                        'color': config.title_red_color,
                                                                        'margin-bottom': "0",
                                                                        'margin-top': "0",
                                                                        'padding-bottom': "0"})
                #self.show_annotations_chkbox.value = False
                self.show_annotations = False
            else:
                self.show_annotations_col[0] = pn.Spacer(height=50)
                self.show_annotations = True
                self.update_synteny_per_pos_plot()
        else:
            self.show_annotations_col[0] = pn.Spacer(height=50)
            self.show_annotations = False
            self.update_synteny_per_pos_plot()

    def filter_synteny_per_pos_plot(self, event):
        #print("\nIn filter_synteny_per_pos_plot")
        self.filter_plot_by_metadata = 1

        feature = self.synteny_per_pos_feature_select.value
        groups = self.synteny_per_pos_groups_select.value

        # There are selected groups
        if len(groups) > 0:
            score_per_pos_contig = self.score_per_pos_contig.copy()
            score_per_pos_contig['Group1'] = score_per_pos_contig['Sample1'].map(self.metadata_dict[feature])
            score_per_pos_contig['Group2'] = score_per_pos_contig['Sample2'].map(self.metadata_dict[feature])

            self.score_per_pos_contig_filtered = score_per_pos_contig[score_per_pos_contig['Group1'].isin(groups) &
                                                                      score_per_pos_contig['Group2'].isin(groups)]

            # Calculate the average synteny scores for each position
            self.avg_score_per_pos_contig_filtered = \
                self.score_per_pos_contig_filtered[['Contig_name', 'Position', 'Synteny_score']]. \
                sort_values(['Position']).groupby('Position') \
                .agg(Count=('Synteny_score', 'size'), Avg_synteny_score=('Synteny_score', 'mean')).reset_index()

            # Fill the missing positions with score=0 (default jump=5000)
            self.avg_score_per_pos_contig_filtered = self.avg_score_per_pos_contig_filtered.\
                merge(how='right', on='Position', right=pd.DataFrame(
                      {'Position': np.arange(self.avg_score_per_pos_contig_filtered.iloc[0]['Position'],
                                             self.avg_score_per_pos_contig_filtered.iloc[-1]['Position'] + 2 *
                                             config.region_length, config.region_length)
                       })).sort_values(by='Position').reset_index().drop(['index'], axis=1)
            self.avg_score_per_pos_contig_filtered['Position'] = self.avg_score_per_pos_contig_filtered['Position'].astype(int)

            hypervar_threshold = self.bottom_percentile
            hypercons_threshold = self.top_percentile
            self.median_counts_filtered = self.avg_score_per_pos_contig_filtered['Count'].median()
            self.bottom_percentile_counts_filtered = self.avg_score_per_pos_contig_filtered['Count'].quantile(0.1)
            print("\nMedian of pairs per region counts (with filtering) is: " + str(self.median_counts_filtered))
            print("Percentile 10 of pairs per region counts (with filtering) is: " +
                  str(self.bottom_percentile_counts_filtered))

            self.avg_score_per_pos_contig_filtered['Hypervariable'] = np.where(
                (self.avg_score_per_pos_contig_filtered['Avg_synteny_score'] <= hypervar_threshold) &
                (self.avg_score_per_pos_contig_filtered['Count'] >= self.bottom_percentile_counts_filtered), 1, 0)
            self.avg_score_per_pos_contig_filtered['Hyperconserved'] = np.where(
                (self.avg_score_per_pos_contig_filtered['Avg_synteny_score'] >= hypercons_threshold) &
                (self.avg_score_per_pos_contig_filtered['Count'] >= self.median_counts_filtered), 1, 0)

        # No groups are selected for filteration -> use the full tables
        else:
            self.score_per_pos_contig_filtered = self.score_per_pos_contig
            self.avg_score_per_pos_contig_filtered = self.avg_score_per_pos_contig

        self.update_synteny_per_pos_plot()

    def reset_filter_synteny_per_pos_plot(self, event):
        self.filter_plot_by_metadata = 0
        self.synteny_per_pos_groups_select.value = []
        self.update_synteny_per_pos_plot()

    def download_synteny_per_pos_plot(self, event):
        fformat = self.synteny_per_pos_image_format.value

        # Update the placeholder of the filenames adding the positions range
        start_pos = self.start_pos_input.value
        end_pos = self.end_pos_input.value
        synteny_per_pos_file = "Synteny_per_position_plot_" + self.ref_genome + "_" + self.contig_name + "_" + \
                               str(start_pos) + "-" + str(end_pos)

        # Update the placeholder of the filename when filtering by metadata, adding the selected groups.
        if self.filter_plot_by_metadata:
            groups = self.synteny_per_pos_groups_select.value
            if len(groups) > 0:
                synteny_per_pos_file += "_filtered_by_" + self.synteny_per_pos_feature_select.value
                for group in self.synteny_per_pos_groups_select.value:
                    synteny_per_pos_file += "_" + str(group)

        self.save_synteny_per_pos_plot_path.placeholder = synteny_per_pos_file

        # Set the directory for saving
        if self.save_synteny_per_pos_plot_path.value == "":
            synteny_per_pos_file_path = self.downloads_dir_path + self.save_synteny_per_pos_plot_path.placeholder + \
                                        "." + fformat
        else:
            synteny_per_pos_file_path = self.save_synteny_per_pos_plot_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, synteny_per_pos_file_path, re.IGNORECASE):
                synteny_per_pos_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(synteny_per_pos_file_path):
                synteny_per_pos_file_path = self.downloads_dir_path + synteny_per_pos_file_path

        self.synteny_per_pos_plot.savefig(synteny_per_pos_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + synteny_per_pos_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_synteny_per_pos_plot_column.pop(4)
        self.download_synteny_per_pos_plot_column.append(download_floatpanel)

    def download_synteny_per_pos_table(self, event):
        fformat = "tsv"

        # Update the placeholder of the filenames adding the positions range
        start_pos = self.start_pos_input.value
        end_pos = self.end_pos_input.value
        synteny_per_pos_table = "Data_for_synteny_per_position_plot_" + self.ref_genome + "_" + self.contig_name + \
                                "_" + str(start_pos) + "-" + str(end_pos)

        # Update the placeholder of the filename when filtering by metadata, adding the selected groups.
        if self.filter_plot_by_metadata:
            groups = self.synteny_per_pos_groups_select.value
            if len(groups) > 0:
                synteny_per_pos_table += "_filtered_by_" + self.synteny_per_pos_feature_select.value
                for group in groups:
                    synteny_per_pos_table += "_" + str(group)

        self.save_synteny_per_pos_table_path.placeholder = synteny_per_pos_table

        # Set the directory for saving
        if self.save_synteny_per_pos_table_path.value == "":
            synteny_per_pos_file_path = self.downloads_dir_path + self.save_synteny_per_pos_table_path.placeholder + \
                                        "." + fformat
        else:
            synteny_per_pos_file_path = self.save_synteny_per_pos_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, synteny_per_pos_file_path, re.IGNORECASE):
                synteny_per_pos_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(synteny_per_pos_file_path):
                synteny_per_pos_file_path = self.downloads_dir_path + synteny_per_pos_file_path

        # The user requested to filter the data by a metadata feature - use the filtered tables
        if self.filter_plot_by_metadata:
            # Consider only the positions within the requested range
            avg_score_per_pos_contig = self.avg_score_per_pos_contig_filtered[
                self.avg_score_per_pos_contig_filtered['Position'] >= int(start_pos)]
            avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

        # No filtering by metadata - use the full tables
        else:
            # Consider only the positions within the requested range
            avg_score_per_pos_contig = self.avg_score_per_pos_contig[
                self.avg_score_per_pos_contig['Position'] >= int(start_pos)]
            avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

        avg_score_per_pos_contig.to_csv(synteny_per_pos_file_path, sep='\t', index=False)

        download_message = "The table is successfully saved under:\n" + synteny_per_pos_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_synteny_per_pos_table_column.pop(3)
        self.download_synteny_per_pos_table_column.append(download_floatpanel)

    def create_multi_genome_column(self):

        self.sample_sizes_slider_multi.value = config.sampling_sizes[0]

        self.genomes_subset_select.options = self.ref_genomes_list_by_pairs_num
        self.selected_genomes_subset = self.ref_genomes_list_by_pairs_num
        self.selected_subset_species_num = len(self.selected_genomes_subset)

        # Doesn't work
        #if self.number_of_genomes > 20:
        #    self.genomes_subset_select.size = 20
        #else:
        #    self.genomes_subset_select.size = self.number_of_genomes

        self.genomes_select_card.append(self.genomes_subset_select)

        all_or_subset_row = pn.Row(self.all_or_subset_radio, pn.Spacer(width=15), self.genomes_sort_select_multi)

        # Create the species selection part (which is common to all input modes)
        species_select_col = pn.Column(
            all_or_subset_row,
            self.genomes_select_card,
            self.update_genomes_selection_button,
            styles={'margin': "0 0 10px 0", 'padding': "0"}
        )
        self.main_multi_column.clear()
        self.main_multi_column.append(species_select_col)

        if self.input_mode == "SynTracker":
            self.create_multi_genomes_column_syntracker_mode()
            self.main_multi_column.append(self.synteny_multi_initial_plots_column)

        elif self.input_mode == "ANI":
            self.create_multi_genomes_column_ANI_mode()
            self.main_multi_column.append(self.ani_multi_plots_column)

        else:
            self.create_multi_genomes_column_both_mode()
            self.main_multi_column.append(self.synteny_ani_multi_tabs)

    def create_multi_genomes_column_both_mode(self):
        before = time.time()
        print("\nStart create_multi_genomes_column_both_mode")

        self.synteny_multi_initial_plots_column.styles = config.both_mode_other_style
        self.ani_multi_plots_column.styles = config.both_mode_other_style

        # Run the task of creating the ANI single view in another thread.
        thread_ani_multi = threading.Thread(target=self.create_multi_genomes_column_ANI_mode)
        thread_ani_multi.start()  # Start the thread

        self.create_multi_genomes_column_syntracker_mode()
        #self.create_multi_genomes_column_ANI_mode()

        self.synteny_ani_multi_tabs.clear()
        self.synteny_ani_multi_tabs.append(('SynTracker', self.synteny_multi_initial_plots_column))
        self.synteny_ani_multi_tabs.append(('ANI', self.ani_multi_plots_column))

        after = time.time()
        duration = after - before
        print("\ncreate_multi_genomes_column_both_mode took " + str(duration) + " seconds.\n")

    def update_genomes_selection(self, event):
        # Build the two bar-plots for subsampled regions based on the selected list of genomes
        if self.all_or_subset_radio.value == 'All species':
            self.selected_genomes_subset = self.genomes_subset_select.options
        else:
            self.selected_genomes_subset = self.genomes_subset_select.value

        self.selected_subset_species_num = len(self.selected_genomes_subset)

        # No selected species -> treat as all species
        if self.selected_subset_species_num == 0:
            self.selected_genomes_subset = self.genomes_subset_select.options

        print("\nupdate_genomes_selection:\nNumber of selected species: " + str(self.selected_subset_species_num))
        #print(self.selected_genomes_subset)

        if self.input_mode == "SynTracker":
            self.create_multi_genomes_column_syntracker_mode()

        elif self.input_mode == "ANI":
            self.create_multi_genomes_column_ANI_mode()

        else:
            self.create_multi_genomes_column_syntracker_mode()
            self.create_multi_genomes_column_ANI_mode()

    def create_multi_genomes_column_syntracker_mode(self):
        self.synteny_multi_initial_plots_column.clear()
        self.plots_by_size_multi_column.clear()

        # Get the score-per-region table for the selected genomes only
        self.score_per_region_genomes_subset_df = dm.return_genomes_subset_table(self.score_per_region_all_genomes_df,
                                                                                 self.selected_genomes_subset)
        # Get the number of genomes from the selected subset that actually appear in the SynTracker file
        presented_genomes_list = self.score_per_region_genomes_subset_df['Ref_genome'].unique()
        self.species_num_in_subset_syntracker = len(presented_genomes_list)
        print("\ncreate_multi_genomes_column_syntracker_mode: Number of available species: " +
              str(self.species_num_in_subset_syntracker))

        # No synteny data for the selected genomes subset
        if self.species_num_in_subset_syntracker == 0:
            text = "None of the species in the selected subset appears in the SynTracker input file."
            self.synteny_multi_initial_plots_column.append(pn.pane.Markdown(text, styles={'font-size': "18px",
                                                                                          'color': config.title_red_color,
                                                                                          'margin': "0"}))
        # At least one selected species appear in the SynTracker file
        else:
            # Create the df for plot presenting the number of pairs vs. subsampled regions
            self.pairs_num_per_sampling_size_multi_genomes_df = \
                dm.create_pairs_num_per_sampling_size(self.score_per_region_genomes_subset_df)

            # If the number of species with 40 sampled regions is smaller than the selected requested of species,
            # present also the 'All regions' bar.
            # If not, do not present this bar (and remove this option from the slider)
            species_at_40 = self.pairs_num_per_sampling_size_multi_genomes_df['Number_of_species'].iloc[1]
            if species_at_40 == self.species_num_in_subset_syntracker:
                self.sample_sizes_slider_multi.options = config.sampling_sizes_wo_all
                self.sample_sizes_slider_multi.value = config.sampling_sizes_wo_all[0]
                is_all_regions = 0
            else:
                self.sample_sizes_slider_multi.options = config.sampling_sizes
                self.sample_sizes_slider_multi.value = config.sampling_sizes[0]
                is_all_regions = 1

            # Create the number of pairs vs. subsampled regions bar plot
            pairs_vs_sampling_size_bar_plot = pn.bind(pm.plot_pairs_vs_sampling_size_bar,
                                                      df=self.pairs_num_per_sampling_size_multi_genomes_df,
                                                      sampling_size=self.sample_sizes_slider_multi,
                                                      is_all_regions=is_all_regions)
            pairs_bar_plot_pane = pn.pane.HoloViews(pairs_vs_sampling_size_bar_plot, width=530, sizing_mode="fixed")

            # Create a markdown for the pairs lost percent (binded to the slider)
            binded_text = pn.bind(widgets.create_pairs_lost_text, self.pairs_num_per_sampling_size_multi_genomes_df,
                                  self.sample_sizes_slider_multi)

            pairs_plot_column = pn.Column(pn.pane.Markdown(refs=binded_text, align='center'), pairs_bar_plot_pane,
                                          styles={'background-color': "white"})

            # Create the number of species vs. subsampled regions bar plot
            species_vs_sampling_size_bar_plot = pn.bind(pm.plot_species_vs_sampling_size_bar,
                                                        df=self.pairs_num_per_sampling_size_multi_genomes_df,
                                                        sampling_size=self.sample_sizes_slider_multi,
                                                        is_all_regions=is_all_regions)
            species_bar_plot_pane = pn.pane.HoloViews(species_vs_sampling_size_bar_plot, width=530, sizing_mode="fixed")

            # Create a markdown for the pairs lost percent (binded to the slider)
            binded_text = pn.bind(widgets.create_species_num_text, self.pairs_num_per_sampling_size_multi_genomes_df,
                                  self.sample_sizes_slider_multi)

            species_plot_column = pn.Column(pn.pane.Markdown(refs=binded_text, align='center'), species_bar_plot_pane,
                                            styles={'background-color': "white"})

            plots_row = pn.Row(pairs_plot_column, pn.Spacer(width=20), species_plot_column)
            slider_row = pn.Row(self.sample_sizes_slider_multi, align='center')
            button_row = pn.Row(self.show_box_plot_multi_button,  align='center')

            plots_column = pn.Column(
                plots_row,
                pn.Spacer(height=20),
                slider_row,
                button_row,
                self.plots_by_size_multi_column,
            )

            self.synteny_multi_initial_plots_column.append(plots_column)

    def create_multi_genomes_plots_by_APSS(self, event):

        # Unwatch watchers (if it's not the first time that this function is called)
        if self.clicked_button_display_APSS_multi and self.is_metadata:
            self.box_plot_feature_select.param.unwatch(self.feature_select_watcher)
        self.clicked_button_display_APSS_multi = 1

        self.sorted_selected_genomes_subset = []
        self.sampling_size_multi = self.sample_sizes_slider_multi.value
        print("\nMulti species visualization. Selected subsampling size = " + self.sampling_size_multi)

        self.species_num_at_sampling_size = self.pairs_num_per_sampling_size_multi_genomes_df.loc[
            self.pairs_num_per_sampling_size_multi_genomes_df['Subsampled_regions'] == self.sampling_size_multi,
            'Number_of_species'].values[0]
        print("Number of species at this sampling size: " + str(self.species_num_at_sampling_size))

        self.plots_by_size_multi_column.clear()
        self.box_plot_card.clear()
        self.metadata_box_plot_card.clear()

        # Np species at the selected sampling size
        if self.species_num_at_sampling_size == 0:
            text = "The data obtained using " + self.sampling_size + " subsampled regions is not sufficient for " \
                                                                     "further processing"
            self.plots_by_size_multi_column.append(pn.pane.Markdown(text, styles={'font-size': "18px",
                                                                                  'color': config.title_red_color,
                                                                                  'margin': "0"}))

        # Enough data to present the plot
        else:
            if self.sampling_size_multi == 'All':
                size_title = "Presenting plots using APSS from all available regions:"
            else:
                size_title = "Presenting plots using APSS from " + self.sampling_size_multi + \
                             " subsampled regions:"

            self.plots_by_size_multi_column.append(pn.pane.Markdown(size_title, styles={'font-size': "20px",
                                                                                        'color': config.title_purple_color,
                                                                                        'margin': "0 0 10px 0",
                                                                                        'padding': "0"}))

            # Check if the requested genome and size have already been calculated. If so, fetch the specific dataframe
            if self.calculated_APSS_all_genomes_size_dict[self.sampling_size_multi]:
                print("\nThe selected size (" + self.sampling_size_multi + ") has already been calculated - retrieve it.")
                all_genomes_selected_size_APSS_df = self.APSS_all_genomes_all_sizes_dict[self.sampling_size_multi]

            else:
                # Calculate and return the dataframe with average scores for the selected genome and sampling size
                print("\nThe selected size (" + self.sampling_size_multi + ") has not been calculated yet - calculate it...")
                before = time.time()
                all_genomes_selected_size_APSS_df = dm.calculate_APSS_all_genomes_sampling_size(
                    self.score_per_region_all_genomes_df, self.sampling_size_multi)
                after = time.time()
                duration = after - before
                print("Calculating APSS with " + str(self.sampling_size_multi ) + " regions for " +
                      str(self.number_of_genomes) + " species took " + str(duration) + " seconds.\n")
                # Save the dataframe in the main dictionary
                self.APSS_all_genomes_all_sizes_dict[self.sampling_size_multi] = all_genomes_selected_size_APSS_df
                self.calculated_APSS_all_genomes_size_dict[self.sampling_size_multi] = 1

            # Get the table containing the selected genomes only
            self.genomes_subset_selected_size_APSS_df = dm.return_genomes_subset_APSS_selected_size_table(
                all_genomes_selected_size_APSS_df, self.selected_genomes_subset)

            # Sort the present genomes list by the requested sorting method
            presented_genomes_list = self.genomes_subset_selected_size_APSS_df['Ref_genome'].unique()
            for genome in self.selected_genomes_subset:
                if genome in presented_genomes_list:
                    self.sorted_selected_genomes_subset.append(genome)

            # Add the plots to the layout
            self.create_box_plot_multi_pane()
            self.plots_by_size_multi_column.append(self.box_plot_card)

    def create_multi_genomes_column_ANI_mode(self):
        before = time.time()
        print("\n\nStart create_multi_genomes_column_ANI_mode in another thread.")

        # Unwatch watchers (if it's not the first time that this function is called)
        if self.visited_ANI_tab_multi and self.is_metadata:
            self.box_plot_feature_select_ani.param.unwatch(self.feature_select_ani_watcher)
        self.visited_ANI_tab_multi = 1

        self.sorted_selected_genomes_subset_ani = []
        self.ani_multi_plots_column.clear()
        self.box_plot_card_ani.clear()
        self.metadata_box_plot_card_ani.clear()

        self.ani_scores_genomes_subset_df = dm.return_genomes_subset_ani_table(self.ani_scores_all_genomes_df,
                                                                               self.selected_genomes_subset)

        presented_genomes_list = self.ani_scores_genomes_subset_df['Ref_genome'].unique()
        self.species_num_in_subset_ani = len(presented_genomes_list)
        print("\nNumber of available species for ANI analysis: " +
              str(self.species_num_in_subset_ani))

        # No data for the selected genomes subset
        if self.species_num_in_subset_ani == 0:
            text = "None of the species in the selected subset appears in the ANI input file."
            self.ani_multi_plots_column.append(pn.pane.Markdown(text, styles={'font-size': "18px",
                                                                              'color': config.title_red_color,
                                                                              'margin': "0"}))

        # Enough data -> creating plots
        else:
            # Sort the present genomes list by the requested sorting method
            for genome in self.selected_genomes_subset:
                if genome in presented_genomes_list:
                    self.sorted_selected_genomes_subset_ani.append(genome)

            # Add the plots to the layout
            self.create_box_plot_multi_pane_ani()
            self.ani_multi_plots_column.append(self.box_plot_card_ani)

        after = time.time()
        duration = after - before
        print("\ncreate_multi_genomes_column_ANI_mode took " + str(duration) + " seconds.\n")

    def create_box_plot_multi_pane(self):
        styling_title = "Plot styling options:"
        metadata_colors_row = pn.Row(self.box_plot_same_color, pn.Spacer(width=3), self.box_plot_different_color)
        metadata_col = pn.Column(self.box_plot_feature_select, metadata_colors_row, styles={'padding': "10x"})
        self.metadata_box_plot_card.append(metadata_col)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.box_plot_color,
                                pn.Spacer(height=5),
                                self.use_metadata_box_plot,
                                self.metadata_box_plot_card
                                )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_box_plot)

        if self.species_num_at_sampling_size == self.number_of_genomes:
            box_plot_file = "Boxplot_all_species_" + self.sampling_size_multi
            boxplot_table = "Data_for_boxplot_all_species_" + self.sampling_size_multi
            pvalues_table = "P-values_all_species_" + self.sampling_size_multi
        else:
            box_plot_file = "Boxplot_" + str(self.species_num_at_sampling_size) + "_species_" + \
                            self.sampling_size_multi + "_regions"
            boxplot_table = "Data_for_boxplot_" + str(self.species_num_at_sampling_size) + "_species_" + \
                            self.sampling_size_multi + "_regions"
            pvalues_table = "P-values_" + str(self.species_num_at_sampling_size) + "_species_" + \
                            self.sampling_size_multi + "_regions"
        self.save_box_plot_file_path.placeholder = box_plot_file
        self.save_boxplot_table_path.placeholder = boxplot_table
        self.save_pvalues_table_path.placeholder = pvalues_table

        self.download_box_plot_column = pn.Column(pn.pane.Markdown(save_file_title,
                                                                   styles={'font-size': "15px",
                                                                           'font-weight': "bold",
                                                                           'color': config.title_blue_color,
                                                                           'margin': "0"}),
                                                  self.box_plot_image_format, self.save_box_plot_file_path,
                                                  download_button, pn.pane.Markdown())

        download_table_button = pn.widgets.Button(name='Download data table in tsv format', button_type='primary')
        download_table_button.on_click(self.download_boxplot_table)

        self.download_boxplot_table_column = pn.Column(self.save_boxplot_table_path, download_table_button,
                                                       pn.pane.Markdown())

        self.download_pvalues_table_column = pn.Column(self.save_pvalues_table_path, self.download_pvalues_button,
                                                       pn.pane.Markdown())

        self.download_multi_col = pn.Column(styling_col, pn.Spacer(height=30), self.download_box_plot_column,
                                            self.download_boxplot_table_column)

        # There is metadata
        if self.is_metadata:

            self.box_plot_feature_select.options = self.metadata_features_list
            self.box_plot_feature_select.value = self.metadata_features_list[0]

            self.calculate_metadata_for_box_plot()

            # Set watcher for the feature-select widget
            self.feature_select_watcher = self.box_plot_feature_select.param.watch(self.update_feature_in_boxplot,
                                                                                   'value', onlychanged=True)

            # Add the p-values download column
            self.download_multi_col.append(self.download_pvalues_table_column)

        # No metadata
        else:
            self.use_metadata_box_plot.disabled = True

        self.box_plot = pn.bind(pm.create_box_plot, avg_df=self.genomes_subset_selected_size_APSS_df,
                                sorted_genomes_list=self.sorted_selected_genomes_subset,
                                pvalues_df=self.boxplot_p_values_df, color=self.box_plot_color,
                                use_metadata=self.use_metadata_box_plot, feature=self.box_plot_feature_select.value,
                                same_color=self.box_plot_same_color, different_color=self.box_plot_different_color)

        self.box_plot_pane = pn.pane.Matplotlib(self.box_plot, width=720, dpi=300, tight=True, format='png')

        box_plot_row = pn.Row(self.download_multi_col, pn.Spacer(width=30), self.box_plot_pane,
                              styles={'padding': "15px"})
        self.box_plot_card.append(box_plot_row)

    def calculate_metadata_for_box_plot(self):
        before = time.time()

        feature = self.box_plot_feature_select.value

        print("\ncalculate_metadata_for_box_plot:")
        print("Number of species to present: " + str(self.species_num_at_sampling_size))
        print("Compared feature: " + feature)

        if self.species_num_at_sampling_size == self.number_of_genomes:
            box_plot_file = "Boxplot_all_species_"
            boxplot_table = "Data_for_boxplot_all_species_"
            pvalues_table = "P-values_all_species_"
        else:
            box_plot_file = "Boxplot_" + str(self.species_num_at_sampling_size) + "_species_"
            boxplot_table = "Data_for_boxplot_" + str(self.species_num_at_sampling_size) + "_species_"
            pvalues_table = "P-values_" + str(self.species_num_at_sampling_size) + "_species_"

        self.save_box_plot_file_path.placeholder = box_plot_file + self.sampling_size_multi + "_regions_feature_" + \
                                                   feature
        self.save_boxplot_table_path.placeholder = boxplot_table + self.sampling_size_multi + "_regions_feature_" + \
                                                   feature
        self.save_pvalues_table_path.placeholder = pvalues_table + self.sampling_size_multi + "_regions_feature_" + \
                                                   feature

        self.genomes_subset_selected_size_APSS_df['Category'] = self.genomes_subset_selected_size_APSS_df.apply(
            lambda row: category_by_feature(row, feature, self.metadata_dict), axis=1)
        #print("\nDF for box_plot with category:")
        #print(self.genomes_subset_selected_size_APSS_df)

        same_feature = 'Same ' + feature
        diff_feature = 'Different ' + feature

        # Calculate P-values for each genome
        valid_pval_list = []
        genome_pval_dict = {}
        genome_effect_size_dict = {}
        pval_corrected = []
        for genome in self.sorted_selected_genomes_subset:
            #print("\nGenome: " + genome + ", Feature: " + feature)
            same_array = self.genomes_subset_selected_size_APSS_df[
                (self.genomes_subset_selected_size_APSS_df['Ref_genome'] == genome) &
                (self.genomes_subset_selected_size_APSS_df['Category'] == same_feature)]['APSS']
            diff_array = self.genomes_subset_selected_size_APSS_df[
                (self.genomes_subset_selected_size_APSS_df['Ref_genome'] == genome) &
                (self.genomes_subset_selected_size_APSS_df['Category'] == diff_feature)]['APSS']
        
            # Sample size is enough for P-value calculation
            if len(same_array) >= 1 and len(diff_array) >= 1:
                before = time.time()
                u, p_val = return_p_value_mannwhitneyu(same_array, diff_array)
                after = time.time()
                duration = after - before
                #print("\nP-value = " + str(p_val))
                #print("P-value calculation took " + str(duration) + " seconds.\n")
                if str(p_val) != "NaN":
                    valid_pval_list.append(p_val)
                    genome_pval_dict[genome] = p_val
                    effect_size = abs(1 - (2 * u) / (len(same_array) * len(diff_array)))
                    genome_effect_size_dict[genome] = effect_size
                    #print("Effect size = " + str(effect_size))
                # Sample size is not enough for P-value calculation
            else:
                genome_pval_dict[genome] = np.nan
                genome_effect_size_dict[genome] = np.nan
                #print("\nCannot calculate P-value for feature " + feature + ": Sample size is not enough.")

        #print("\nOriginal valid p-values:")
        #print(valid_pval_list)

        if len(valid_pval_list) >= 2:
            reject, pval_corrected, _, q_values = multipletests(valid_pval_list, method='fdr_bh')
            #print("Corrected p-values:")
            #print(pval_corrected)

        valid_counter = 0
        if len(pval_corrected) > 0:
            for genome in self.sorted_selected_genomes_subset:
                if not np.isnan(genome_pval_dict[genome]):
                    genome_pval_dict[genome] = pval_corrected[valid_counter]
                    valid_counter += 1

        updated_pval_list = genome_pval_dict.values()
        effect_size_list = genome_effect_size_dict.values()
        genomes_pvalues_dict = {'Ref_genome': self.sorted_selected_genomes_subset,
                                'P_value': updated_pval_list, 'Effect_size': effect_size_list}
        self.boxplot_p_values_df = pd.DataFrame(genomes_pvalues_dict)

        self.boxplot_p_values_df['Significance'] = self.boxplot_p_values_df.apply(lambda row: return_significance(row),
                                                                                  axis=1)
        self.boxplot_p_values_df['P_value'] = \
            self.boxplot_p_values_df['P_value'].map(lambda v: f"{v: .1e}" if pd.notna(v) else v)
        self.boxplot_p_values_df['Effect_size'] = \
            self.boxplot_p_values_df['Effect_size'].map(lambda v: round(v, 2) if pd.notna(v) else v)

        print("\n")
        print(self.boxplot_p_values_df)

        after = time.time()
        duration = after - before
        print("\ncalculate_metadata_for_box_plot took " + str(duration) + " seconds.\n")

    def update_feature_in_boxplot(self, event):
        self.calculate_metadata_for_box_plot()

        self.box_plot = pn.bind(pm.create_box_plot, avg_df=self.genomes_subset_selected_size_APSS_df,
                                sorted_genomes_list=self.sorted_selected_genomes_subset,
                                pvalues_df=self.boxplot_p_values_df, color=self.box_plot_color,
                                use_metadata=self.use_metadata_box_plot, feature=self.box_plot_feature_select.value,
                                same_color=self.box_plot_same_color, different_color=self.box_plot_different_color)
        self.box_plot_pane.object = self.box_plot

    def download_box_plot(self, event):
        fformat = self.box_plot_image_format.value

        # Set the directory for saving
        if self.save_box_plot_file_path.value == "":
            box_plot_file_path = self.downloads_dir_path + self.save_box_plot_file_path.placeholder + "." + fformat
        else:
            box_plot_file_path = self.save_box_plot_file_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, box_plot_file_path, re.IGNORECASE):
                box_plot_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(box_plot_file_path):
                box_plot_file_path = self.downloads_dir_path + box_plot_file_path

        self.box_plot().savefig(box_plot_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + box_plot_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_box_plot_column.pop(4)
        self.download_box_plot_column.append(download_floatpanel)

    def download_boxplot_table(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_boxplot_table_path.value == "":
            boxplot_table_path = self.downloads_dir_path + self.save_boxplot_table_path.placeholder + "." + fformat
        else:
            boxplot_table_path = self.save_boxplot_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, boxplot_table_path, re.IGNORECASE):
                boxplot_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(boxplot_table_path):
                boxplot_table_path = self.downloads_dir_path + boxplot_table_path

        self.genomes_subset_selected_size_APSS_df.to_csv(boxplot_table_path, index=False, sep='\t')

        download_message = "The table is successfully saved under:\n" + boxplot_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_boxplot_table_column.pop(2)
        self.download_boxplot_table_column.append(download_floatpanel)

    def download_pvalues_table(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_pvalues_table_path.value == "":
            pvalues_table_path = self.downloads_dir_path + self.save_pvalues_table_path.placeholder + "." + fformat
        else:
            pvalues_table_path = self.save_pvalues_table_path.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, pvalues_table_path, re.IGNORECASE):
                pvalues_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(pvalues_table_path):
                pvalues_table_path = self.downloads_dir_path + pvalues_table_path

        self.boxplot_p_values_df.to_csv(pvalues_table_path, index=False, sep='\t')

        download_message = "The table is successfully saved under:\n" + pvalues_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_pvalues_table_column.pop(2)
        self.download_pvalues_table_column.append(download_floatpanel)

    def create_box_plot_multi_pane_ani(self):
        styling_title = "Plot styling options:"
        metadata_colors_row = pn.Row(self.box_plot_same_color_ani, pn.Spacer(width=3), self.box_plot_different_color_ani)
        metadata_col = pn.Column(self.box_plot_feature_select_ani, metadata_colors_row, styles={'padding': "10x"})
        self.metadata_box_plot_card_ani.append(metadata_col)
        styling_col = pn.Column(pn.pane.Markdown(styling_title, styles={'font-size': "15px", 'font-weight': "bold",
                                                                        'color': config.title_blue_color,
                                                                        'margin': "0"}),
                                self.box_plot_color_ani,
                                pn.Spacer(height=5),
                                self.use_metadata_box_plot_ani,
                                self.metadata_box_plot_card_ani
                                )

        save_file_title = "Plot download options:"
        download_button = pn.widgets.Button(name='Download high-resolution image', button_type='primary')
        download_button.on_click(self.download_box_plot_ani)

        if self.species_num_in_subset_ani == self.number_of_genomes:
            num = "all"
        else:
            num = str(self.species_num_in_subset_ani)
        box_plot_file = "ANI_boxplot_" + num + "_species"
        boxplot_table = "Data_for_ANI_boxplot_" + num + "_species"
        pvalues_table = "P-values_ANI_boxplot_" + num + "_species"

        self.save_box_plot_file_path_ani.placeholder = box_plot_file
        self.save_boxplot_table_path_ani.placeholder = boxplot_table
        self.save_pvalues_table_path_ani.placeholder = pvalues_table

        self.download_box_plot_column_ani = pn.Column(pn.pane.Markdown(save_file_title,
                                                                       styles={'font-size': "15px",
                                                                               'font-weight': "bold",
                                                                               'color': config.title_blue_color,
                                                                               'margin': "0"}),
                                                      self.box_plot_image_format_ani, self.save_box_plot_file_path_ani,
                                                      download_button, pn.pane.Markdown())

        download_table_button = pn.widgets.Button(name='Download data table in tsv format', button_type='primary')
        download_table_button.on_click(self.download_boxplot_table_ani)

        self.download_boxplot_table_column_ani = pn.Column(self.save_boxplot_table_path_ani, download_table_button,
                                                           pn.pane.Markdown())

        self.download_pvalues_table_column_ani = pn.Column(self.save_pvalues_table_path_ani,
                                                           self.download_pvalues_button_ani,
                                                           pn.pane.Markdown())

        self.download_multi_col_ani = pn.Column(styling_col, pn.Spacer(height=30), self.download_box_plot_column_ani,
                                                self.download_boxplot_table_column_ani)

        # There is metadata
        if self.is_metadata:

            self.box_plot_feature_select_ani.options = self.metadata_features_list
            self.box_plot_feature_select_ani.value = self.metadata_features_list[0]

            self.calculate_metadata_for_box_plot_ani()

            # Set watcher for the feature-select widget
            self.feature_select_ani_watcher = self.box_plot_feature_select_ani.param.watch(
                self.update_feature_in_boxplot_ani, 'value', onlychanged=True)

            # Add the p-values download column
            self.download_multi_col_ani.append(self.download_pvalues_table_column_ani)

        # No metadata
        else:
            self.use_metadata_box_plot_ani.disabled = True

        self.box_plot_ani = pn.bind(pm.create_box_plot_ani, ani_df=self.ani_scores_genomes_subset_df,
                                    sorted_genomes_list=self.sorted_selected_genomes_subset_ani,
                                    pvalues_df=self.boxplot_p_values_df_ani, color=self.box_plot_color_ani,
                                    use_metadata=self.use_metadata_box_plot_ani,
                                    feature=self.box_plot_feature_select_ani.value,
                                    same_color=self.box_plot_same_color_ani,
                                    different_color=self.box_plot_different_color_ani)

        self.box_plot_pane_ani = pn.pane.Matplotlib(self.box_plot_ani, width=720, dpi=300, tight=True, format='png')

        box_plot_row = pn.Row(self.download_multi_col_ani, pn.Spacer(width=30), self.box_plot_pane_ani,
                              styles={'padding': "15px"})
        self.box_plot_card_ani.append(box_plot_row)

    def calculate_metadata_for_box_plot_ani(self):

        before = time.time()

        #presented_genomes_list = self.ani_scores_genomes_subset_df['Ref_genome'].unique()
        feature = self.box_plot_feature_select_ani.value

        print("\ncalculate_metadata_for_box_plot_ani:")
        print("Number of species to present: " + str(self.species_num_in_subset_ani))
        print("Compared feature: " + feature)

        if self.species_num_in_subset_ani == self.number_of_genomes:
            num = "all"
        else:
            num = str(self.species_num_in_subset_ani)

        box_plot_file = "ANI_boxplot_" + num + "_species"
        boxplot_table = "Data_for_ANI_boxplot_" + num + "_species"
        pvalues_table = "P-values_ANI_boxplot_" + num + "_species"

        self.save_box_plot_file_path_ani.placeholder = box_plot_file + "_feature_" + feature
        self.save_boxplot_table_path.placeholder = boxplot_table + "_feature_" + feature
        self.save_pvalues_table_path.placeholder = pvalues_table + "_feature_" + feature

        self.ani_scores_genomes_subset_df['Category'] = self.ani_scores_genomes_subset_df.apply(
            lambda row: category_by_feature(row, feature, self.metadata_dict), axis=1)
        #print("\nDF for box_plot with category:")
        #print(self.genomes_subset_selected_size_APSS_df)

        same_feature = 'Same ' + feature
        diff_feature = 'Different ' + feature

        # Calculate P-value and effect size for each genome
        valid_pval_list = []
        genome_pval_dict = {}
        genome_effect_size_dict = {}
        pval_corrected = []
        for genome in self.sorted_selected_genomes_subset_ani:
            #print("\nGenome: " + genome + ", Feature: " + feature)
            same_array = self.ani_scores_genomes_subset_df[
                (self.ani_scores_genomes_subset_df['Ref_genome'] == genome) &
                (self.ani_scores_genomes_subset_df['Category'] == same_feature)]['ANI']
            diff_array = self.ani_scores_genomes_subset_df[
                (self.ani_scores_genomes_subset_df['Ref_genome'] == genome) &
                (self.ani_scores_genomes_subset_df['Category'] == diff_feature)]['ANI']

            # Sample size is enough for P-value calculation
            if len(same_array) >= 1 and len(diff_array) >= 1:
                before = time.time()
                u, p_val = return_p_value_mannwhitneyu(same_array, diff_array)
                after = time.time()
                duration = after - before
                #print("\nP-value = " + str(p_val))
                #print("P-value calculation took " + str(duration) + " seconds.\n")
                if str(p_val) != "NaN":
                    valid_pval_list.append(p_val)
                    genome_pval_dict[genome] = p_val
                    effect_size = abs(1 - (2 * u) / (len(same_array) * len(diff_array)))
                    genome_effect_size_dict[genome] = effect_size
                    #print("Effect size = " + str(effect_size))
            # Sample size is not enough for P-value calculation
            else:
                genome_pval_dict[genome] = np.nan
                genome_effect_size_dict[genome] = np.nan
                #print("\nCannot calculate P-value for feature " + feature + ": Sample size is not enough.")

        #print("\nOriginal p-values:")
        #print(valid_pval_list)

        # Need to apply multiple testing correction
        if len(valid_pval_list) >= 2:
            reject, pval_corrected, _, q_values = multipletests(valid_pval_list, method='fdr_bh')
        else:
            pval_corrected = valid_pval_list
        #print("\nCorrected p-values:")
        #print(pval_corrected)

        valid_counter = 0
        if len(pval_corrected) > 0:
            for genome in self.sorted_selected_genomes_subset_ani:
                if not np.isnan(genome_pval_dict[genome]):
                    genome_pval_dict[genome] = pval_corrected[valid_counter]
                    valid_counter += 1

        updated_pval_list = genome_pval_dict.values()
        effect_size_list = genome_effect_size_dict.values()
        genomes_pvalues_dict = {'Ref_genome': self.sorted_selected_genomes_subset_ani,
                                'P_value': updated_pval_list, 'Effect_size': effect_size_list}
        self.boxplot_p_values_df_ani = pd.DataFrame(genomes_pvalues_dict)

        self.boxplot_p_values_df_ani['Significance'] = self.boxplot_p_values_df_ani.apply(
            lambda row: return_significance(row), axis=1)

        self.boxplot_p_values_df_ani['P_value'] = \
            self.boxplot_p_values_df_ani['P_value'].map(lambda v: f"{v: .2e}" if pd.notna(v) else v)
        self.boxplot_p_values_df_ani['Effect_size'] = \
            self.boxplot_p_values_df_ani['Effect_size'].map(lambda v: round(v, 2) if pd.notna(v) else v)

        print(self.boxplot_p_values_df_ani)

        after = time.time()
        duration = after - before
        print("\ncalculate_metadata_for_box_plot_ani took " + str(duration) + " seconds.\n")

    def update_feature_in_boxplot_ani(self, event):
        self.calculate_metadata_for_box_plot_ani()

        self.box_plot_ani = pn.bind(pm.create_box_plot_ani, ani_df=self.ani_scores_genomes_subset_df,
                                    sorted_genomes_list=self.sorted_selected_genomes_subset_ani,
                                    pvalues_df=self.boxplot_p_values_df_ani, color=self.box_plot_color_ani,
                                    use_metadata=self.use_metadata_box_plot_ani,
                                    feature=self.box_plot_feature_select_ani.value,
                                    same_color=self.box_plot_same_color_ani,
                                    different_color=self.box_plot_different_color_ani)
        self.box_plot_pane_ani.object = self.box_plot_ani

    def download_box_plot_ani(self, event):
        fformat = self.box_plot_image_format_ani.value

        # Set the directory for saving
        if self.save_box_plot_file_path_ani.value == "":
            box_plot_file_path = self.downloads_dir_path + self.save_box_plot_file_path_ani.placeholder + "." + fformat
        else:
            box_plot_file_path = self.save_box_plot_file_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, box_plot_file_path, re.IGNORECASE):
                box_plot_file_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(box_plot_file_path):
                box_plot_file_path = self.downloads_dir_path + box_plot_file_path

        self.box_plot_ani().savefig(box_plot_file_path, format=fformat, dpi=300.0, bbox_inches='tight')

        download_message = "The image is successfully saved under:\n" + box_plot_file_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_box_plot_column_ani.pop(4)
        self.download_box_plot_column_ani.append(download_floatpanel)

    def download_boxplot_table_ani(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_boxplot_table_path_ani.value == "":
            boxplot_table_path = self.downloads_dir_path + self.save_boxplot_table_path_ani.placeholder + "." + fformat
        else:
            boxplot_table_path = self.save_boxplot_table_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, boxplot_table_path, re.IGNORECASE):
                boxplot_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(boxplot_table_path):
                boxplot_table_path = self.downloads_dir_path + boxplot_table_path

        self.ani_scores_genomes_subset_df.to_csv(boxplot_table_path, index=False, sep='\t')

        download_message = "The table is successfully saved under:\n" + boxplot_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_boxplot_table_column_ani.pop(2)
        self.download_boxplot_table_column_ani.append(download_floatpanel)

    def download_pvalues_table_ani(self, event):
        fformat = "tsv"

        # Set the directory for saving
        if self.save_pvalues_table_path_ani.value == "":
            pvalues_table_path = self.downloads_dir_path + self.save_pvalues_table_path_ani.placeholder + "." + fformat
        else:
            pvalues_table_path = self.save_pvalues_table_path_ani.value

            # Add a file-format suffix if there is none
            regex = r"^\S+\." + re.escape(fformat) + r"$"
            if not re.search(regex, pvalues_table_path, re.IGNORECASE):
                pvalues_table_path += "." + fformat

            # If path is not absolute - save file basename under the downloads dir
            if not os.path.isabs(pvalues_table_path):
                pvalues_table_path = self.downloads_dir_path + pvalues_table_path

        self.boxplot_p_values_df_ani.to_csv(pvalues_table_path, index=False, sep='\t')

        download_message = "The table is successfully saved under:\n" + pvalues_table_path
        markdown = pn.pane.Markdown(download_message, styles={'font-size': "12px", 'color': config.title_red_color})
        download_floatpanel = pn.layout.FloatPanel(markdown, name='Download message', margin=10)
        self.download_pvalues_table_column_ani.pop(2)
        self.download_pvalues_table_column_ani.append(download_floatpanel)


