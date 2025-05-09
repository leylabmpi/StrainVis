import bokeh.palettes as bp
import colorcet as cc

col_set = ['Ref_genome', 'Sample1', 'Sample2', 'Region', 'Synteny_score']
sampling_sizes = ['All', '40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
sampling_sizes_wo_all = ['40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
genomes_sorting_options = ['Number of compared pairs', 'Genome name']
contig_sorting_options = ['Contig length', 'Contig name']
catplot_types = ['Scatter (jitter) plot', 'Boxplot']
min_pairs_for_all_regions = 100
max_clustermap_cols = 120
max_network_nodes = 300
network_iterations_options = ['50', '100', '150', '200', '250', '300', '350', '400', '450', '500']
network_thresholds_options = ['Mean APSS', 'Mean APSS+1 STD', 'Mean APSS+2 STD', 'Define another threshold']
APSS_connections_threshold_default = 0.9
region_length = 5000
top_percentile = 0.9
bottom_percentile = 0.1

file_upload_timeout = 20
downloads_dir = "/Downloads/"
manual_file = "/SynTrackerVis_app/manual.md"

## CSS Styles ##
header_color = "#0072b5"
normal_bar_color = "#B048B5"
highlight_bar_color = "#43BFC7"
title_red_color = "#800517"
title_purple_color = "#800080"
title_blue_color = "#002060"
same_color = "#F22C5D"
diff_color = "#47A3E1"
nodes_default_color = 'gray'

header_container_style = {
    'margin': '0 auto',
    'padding': '0',
    'width': "1350px",
}

menu_row_style = {
    'margin': '0',
    'padding': '0',
}

menu_tabs_style = {
    'font-size': "20px",
}

main_area_style = {
    'margin': '0 auto',
    'width': "1200px",
}

single_multi_tabs_style = {
    'font-size': "24px",
    'background': "#f9f9f9",
}

single_tabs_style = {
    'width': "1160px",
    'font-size': "20px",
    'background': "#f3f3f3",
    'padding': "0px",
}

main_column_style = {
    'width': "1200px",
    'background': "#f9f9f9",
    'padding': "20px",
}

plot_card_style = {
    'background': "#ffffff",
    'width': "1150px",
    'font-size': "20px"
}

secondary_button = {
    'background': 'rgba(0, 128, 255, 0.5)',
    'color': 'white'
}

# Export file formats
matplotlib_file_formats = ['png', 'pdf', 'svg', 'eps']
bokeh_file_formats = ['png', 'svg']

# Colormaps
clustermap_colormaps_list = ['Blues', 'Purples', 'Greens', 'Oranges', 'Reds', 'Greys',
                             'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'BuGn', 'YlGn',
                             'YlGnBu', 'PuBuGn', 'YlOrRd',
                             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia']

categorical_colormap_dict = {
    'Glasbey': cc.glasbey,
    'Glasbey_light': cc.glasbey_light,
    'Glasbey_category10': bp.Category10[10],
    'Set1': bp.Set1[9],
    'Set3': bp.Set3[12],
    'Define custom colormap': ['#000000']
}

continuous_colormap_dict = {
    'Rainbow4': cc.m_rainbow4,
    'Isolum': cc.isolum,
    'Plasma': bp.Plasma256,
    'Viridis': bp.Viridis256,
    'Blues': bp.Blues256,
    'Reds': bp.Reds256,
    'Greens': bp.Greens256,
}





