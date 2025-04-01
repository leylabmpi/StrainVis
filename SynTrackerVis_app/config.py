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

file_upload_timeout = 20
downloads_dir = "/Downloads/"

## CSS Styles ##
normal_bar_color = "#B048B5"
#highlight_bar_color = "#ba2649"
highlight_bar_color = "#43BFC7"
title_red_color = "#800517"
title_purple_color = "#800080"
title_blue_color = "#002060"
same_color = "#F22C5D"
diff_color = "#47A3E1"
nodes_default_color = 'gray'


main_area_style = {
    'width': "1200px",
    'padding': "20px",
}

single_multi_tabs_style = {
    'width': "1200px",
    'font-size': "24px",
    'background': "#f9f9f9",
    'padding': "0px",
}

single_tabs_style = {
    'width': "1200px",
    'font-size': "20px",
    'background': "#f3f3f3",
    'padding': "0px",
}

main_column_style = {
    'background': "#f9f9f9",
    'padding': "20px",
}

plot_card_style = {
    'background': "#ffffff",
    'width': "1150px",
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
    'Set3': bp.Set3[12]
}

continuous_colormap_dict = {
    'cet_rainbow4': cc.m_rainbow4,
    'cet_isolum': cc.isolum,
    'plasma': bp.Plasma256,
    'viridis': bp.Viridis256,
    'Blues': bp.Blues256,
    'Reds': bp.Reds256,
    'Greens': bp.Greens256,
}





