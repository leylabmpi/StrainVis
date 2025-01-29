import bokeh.palettes as bp

sampling_sizes = ['All_regions', '40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
sampling_sizes_wo_all = ['40', '60', '80', '100', '125', '150', '175', '200', '250', '300', '350', '400']
min_pairs_for_all_regions = 100
max_clustermap_cols = 120
network_iterations_options = ['50', '100', '150', '200', '250', '300', '350', '400', '450', '500']
network_thresholds_options = ['Mean APSS', 'Mean APSS+1 STD', 'Mean APSS+2 STD', 'Define another threshold']
APSS_connections_threshold_default = 0.9
region_length = 5000

file_upload_timeout = 20
downloads_dir = "/Downloads/"

## CSS Styles ##
normal_bar_color = "#B048B5"
highlight_bar_color = "#ba2649"
title_red_color = "#800517"
title_purple_color = "#800080"
title_blue_color = "#002060"

main_area_style = {
    'width': "1200px",
    'padding': "20px",
    #'background': "#b0e0e6",
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
    #'background': "#f9f9f9",
    'padding': "20px",
    #'padding': "0px 20px 10px 20px"
}

plot_card_style = {
    'background': "#ffffff",
    'width': "1150px",
}

## matplotlib patrameters
clustermap_colormaps_list = ['Blues', 'Purples', 'Greens', 'Oranges', 'Reds', 'Greys',
                             'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'BuGn', 'YlGn',
                             'YlGnBu', 'PuBuGn', 'YlOrRd',
                             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia']
matplotlib_file_formats = ['png', 'pdf', 'svg', 'eps']

## Bokeh patrameters
Bokeh_categorical_colormap_dict = {
    'Category20': bp.Category20[20],
    'Category10': bp.Category10[10],
    'Set1': bp.Set1[9],
    'Set3': bp.Set3[12],
    'Spectral': bp.Spectral[11],
    'Bokeh': bp.Bokeh[8],
}
Bokeh_categorical_colormap_list = [bp.Category10[10], bp.Category20[20], bp.Pastel1[9], bp.Set1[9], bp.Set3[12],
                                   bp.Spectral[11], bp.Bokeh[8], bp.Turbo256]
Bokeh_continuous_colormap_dict = {
    'Turbo256': bp.Turbo256,
    'Plasma': bp.Plasma256,
    'Viridis': bp.Viridis256,
    'Blues': bp.Blues256,
    'Reds': bp.Reds256,
    'Greens': bp.Greens256,
}
bokeh_file_formats = ['png', 'svg']

normal_bar_color = "#B048B5"
highlight_bar_color = "#43BFC7"




