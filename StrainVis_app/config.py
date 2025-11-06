import bokeh.palettes as bp
import colorcet as cc

col_set = ['Ref_genome', 'Sample1', 'Sample2', 'Region', 'Synteny_score']
ANI_col_names = ['Ref_genome', 'Sample1', 'Sample2', 'ANI']
sampling_sizes = ['All', '40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
sampling_sizes_wo_all = ['40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
genomes_sorting_options = ['Number of compared pairs', 'Species name']
contig_sorting_options = ['Contig length', 'Contig name']
catplot_types = ['Scatter (jitter) plot', 'Boxplot']
clustering_methods = ['euclidean', 'correlation', 'cityblock', 'cosine']
min_pairs_for_all_regions = 100
max_clustermap_cols = 200
region_length = 5000
top_percentile = 0.9
bottom_percentile = 0.1

## Network-related parameters
network_iterations_options = ['50', '100', '150', '200', '250', '300', '350', '400', '450', '500']
network_thresholds_options = ['Mean APSS', 'Mean APSS+1 STD', 'Mean APSS+2 STD', 'Define another threshold']
max_network_nodes = 300
APSS_connections_threshold_default = 0.9
hvnx_nodes_size = 20
hvnx_highlighted_nodes_size = 30
nx_nodes_size = 140
nx_highlighted_nodes_size = 280
outline_width = 0.5
highlighted_outline_width = 3
highlighted_outline_width_matplotlib = 2
outline_color = 'black'
highlighted_outline_color = 'cyan'
max_groups_for_legend = 15

file_upload_timeout = 20
downloads_dir = "/Downloads/"
#manual_file = "/StrainVis_app/manual.md"
manual_file = "/StrainVis_app/manual.html"

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
conserved_color = '#E66B77'
variable_color = '#00ffff'

header_container_style = {
    'margin': '0 auto',
    'padding': '0',
    'width': "1300px",
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

main_column_style = {
    'width': "1200px",
    'background': "#f6f6f6",
    'padding': "20px",
    'margin': "0",
    'border-top': "2px solid #0072b5",
}

ani_or_multi_style = {
    'width': "1160px",
    'background': "#eaeaea",
    'margin': "0",
    'padding': "20px",
}

synteny_only_single_style = {
    'width': "1160px",
    'background': "#eaeaea",
    'margin': "0",
}

both_mode_SynTracker_single_style = {
    'width': "1160px",
    'background': "#eaeaea",
    'padding-top': "7px",
    'margin': "0",
    'border-top': "2px solid #0072b5",
}

both_mode_other_style = {
    'width': "1160px",
    'background': "#eaeaea",
    'margin': "0",
    'padding': "20px",
    'border-top': "2px solid #0072b5",
}

single_tabs_style = {
    'width': "1160px",
    'padding': "20px",
    'margin': "0",
    'border-top': "2px solid #800080",
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
    'cet_glasbey': cc.glasbey,
    'cet_glasbey_light': cc.glasbey_light,
    'cet_glasbey_category10': bp.Category10[10],
    'Set1': bp.Set1[9],
    'Set3': bp.Set3[12],
    'Define custom colormap': ['#000000']
}

continuous_colormap_dict = {
    'plasma_r': bp.Plasma256[::-1],
    'cet_rainbow4_r': cc.rainbow4[::-1],
    'cet_isolum_r': cc.isolum[::-1],
    'viridis_r': bp.Viridis256[::-1],
    'Blues_r': bp.Blues256[::-1],
    'Reds_r': bp.Reds256[::-1],
    'Greens_r': bp.Greens256[::-1],
    'Oranges_r': bp.Oranges256[::-1],
    'Purples_r': bp.Purples256[::-1],
}





