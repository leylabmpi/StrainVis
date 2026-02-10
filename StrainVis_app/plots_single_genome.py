import time
import re
import numpy as np
import pandas as pd
import networkx as nx
import hvplot.pandas  # Enable interactive
import holoviews as hv
import hvplot.networkx as hvnx
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import seaborn as sns
from bokeh.plotting import figure
from bokeh.transform import jitter
from bokeh.models import HoverTool
import StrainVis_app.config as config
hv.extension('bokeh')


def plot_pairs_vs_sampling_size_bar(df, sampling_size, is_all_regions):
    normal_bar_color = "#B048B5"
    highlight_bar_color = "#43BFC7"

    value = str(sampling_size)

    # If 'All regions' shouldn't be presented, remove this row from the df
    if is_all_regions == 0:
        df = df.drop(0)

    #print("plot_pairs_vs_sampling_size_bar:")
    #print(df)

    df['color'] = np.where(df['Subsampled_regions'] == value, highlight_bar_color, normal_bar_color)

    tooltips = [
        ('Number of pairs', '@Number_of_pairs'),
    ]
    hover = HoverTool(tooltips=tooltips)

    bar_plot = df.hvplot.bar(x='Subsampled_regions', y='Number_of_pairs', color='color', xlabel='Subsampled regions',
                             ylabel='Number of pairs').opts(shared_axes=False, tools=[hover])

    return bar_plot


def plot_samples_vs_sampling_size_bar(df, sampling_size, is_all_regions):
    value = str(sampling_size)

    # If 'All regions' shouldn't be presented, remove this row from the df
    if is_all_regions == 0:
        df = df.drop(0)

    df['color'] = np.where(df['Subsampled_regions'] == value, config.highlight_bar_color, config.normal_bar_color)

    tooltips = [
        ('Number of samples', '@Number_of_samples'),
    ]
    hover = HoverTool(tooltips=tooltips)

    bar_plot = df.hvplot.bar(x='Subsampled_regions', y='Number_of_samples', color='color', xlabel='Subsampled regions',
                             ylabel='Number of samples').opts(shared_axes=False, tools=[hover])

    return bar_plot


def create_jitter_plot_bokeh(avg_df, color):
    df_for_jitter = avg_df[['APSS']]
    df_for_jitter.insert(1, 'Category', 'Comparisons')
    #print("\nDF for jitter plot:")
    #print(df_for_jitter)

    plot = figure(width=300, height=600, x_range=['Comparisons'])
    plot.scatter(x=jitter('Category', width=0.5, range=plot.x_range), y='APSS', size=9, color=color,
                 alpha=0.4, source=df_for_jitter)
    plot.xgrid.grid_line_color = None
    plot.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    plot.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
    plot.yaxis.axis_label = "Average Synteny Score"
    plot.xaxis.axis_label_text_font_size = "18pt"

    print("\ncreate_jitter_plot:")
    print(plot)

    return plot


def create_clustermap(matrix, type, cmap, method, is_metadata, feature, is_continuous, cmap_metadata, custom_cmap,
                      metadata_dict):

    # The number of columns doesn't exceed the defined maximum - continue creating the clustermap plot
    col_num = len(matrix.columns)

    # Create a masking matrix to display the NA values in grey color
    mask_array = np.isnan(matrix)

    # Fill the NA cells with the mean score
    matrix = matrix.fillna(matrix.mean().mean())

    #print("\ncreate_clustermap: selected cmap: " + cmap)
    colmap = plt.get_cmap(cmap)
    colmap.set_bad("lightgrey")

    if is_metadata:

        # In case of numeric continuous feature:
        if is_continuous:
            print("\ncreate_clustermap: Continuous feature")

            # Extract non-missing values
            non_missing_values = [metadata_dict[feature][sample] for sample in matrix.iloc[:, 0].index
                                  if not np.isnan(metadata_dict[feature][sample])]

            # Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            print("Min value: " + str(min_value))
            print("Max value: " + str(max_value))

            feature_array = np.array([metadata_dict[feature][sample] for sample in matrix.iloc[:, 0].index])
            normalized_values = (feature_array - min_value) / (max_value - min_value)

            # Get the colormap for the metadata feature
            cmap_metadata_mpl = plt.get_cmap(cmap_metadata)

            # Map the feature values to colors using the colormap
            feature_colors = []
            for i, sample in enumerate(matrix.iloc[:, 0].index):
                value = metadata_dict[feature][sample]
                if np.isnan(value):  # If the value is NaN, use the default color
                    feature_colors.append(config.nodes_default_color)
                else:
                    feature_colors.append(cmap_metadata_mpl(normalized_values[i]))  # Apply the colormap to the value
            feature_colors_df = pd.DataFrame({feature: feature_colors}, index=matrix.index)

            clustermap = sns.clustermap(matrix, metric=method, cmap=colmap, row_cluster=True, linewidths=.5,
                                        cbar_pos=(0.04, 0.82, 0.02, 0.15), xticklabels=1, yticklabels=1,
                                        dendrogram_ratio=(0.2, 0.2), mask=mask_array, row_colors=feature_colors_df)

            # Create a ScalarMappable to associate the colormap with the normalized values
            norm = Normalize(vmin=min_value, vmax=max_value)
            sm = plt.cm.ScalarMappable(cmap=cmap_metadata_mpl, norm=norm)
            sm.set_array([])  # Empty array as we don't need to pass actual data to the ScalarMappable

            # Add the colorbar next to the heatmap
            # [left, bottom, width, height] in figure coordinates
            cbar_ax = clustermap.fig.add_axes([0.95, 0.82, 0.02, 0.15])

            cbar = clustermap.fig.colorbar(
                sm, cax=cbar_ax, orientation='vertical'
            )

            # Customize the colorbar label
            cbar.set_label(label=feature, labelpad=8)

        # Categorical feature
        else:

            # Get the unique groups of the selected feature
            unique_groups = sorted(list(set([str(metadata_dict[feature][sample]) for sample in matrix.iloc[:, 0].index])))
            groups_num = len(unique_groups)

            # Move the 'nan' group (if any) to the end of the list
            if 'nan' in unique_groups:
                unique_groups.remove('nan')
                unique_groups.append('nan')
            #print(unique_groups)

            # If the user defined a custom cmap - process it and turn it into a cmap
            if cmap_metadata == 'Define custom colormap':
                if custom_cmap != "":
                    custom_colors_list = custom_cmap
                    cmap_metadata_mpl = re.split(r'\s*,\s*', custom_colors_list)
                    cmap_length = len(cmap_metadata_mpl)
                    print("\nCustom cmap:")
                    print(cmap_metadata_mpl)

                    # Assign each group with a color, according to the colormap order
                    group_to_color = {group: cmap_metadata_mpl[i % cmap_length] for i, group in enumerate(unique_groups)}

                # Custom colormap is not defined yet
                else:
                    print("\nCustom cmap is not defined yet")
                    # Assign each group with black, according to the colormap order
                    group_to_color = {group: 'black' for i, group in enumerate(unique_groups)}

            # Standard colormap
            else:
                cmap_metadata_mpl = plt.get_cmap(cmap_metadata)
                cmap_length = cmap_metadata_mpl.N

                # Assign each group with a color, according to the colormap order
                group_to_color = {group: cmap_metadata_mpl(i % cmap_length) for i, group in enumerate(unique_groups)}

            #print("cmap length: " + str(cmap_length))

            # Create a colors dataframe, with sample names as indices
            colors = [group_to_color[str(metadata_dict[feature][sample])] for sample in matrix.iloc[:, 0].index]
            colors_df = pd.DataFrame({feature: colors}, index=matrix.index)
            #print(colors_df)

            clustermap = sns.clustermap(matrix, metric=method, cmap=colmap, row_cluster=True, linewidths=.5,
                                        cbar_pos=(0.04, 0.82, 0.02, 0.15), xticklabels=1, yticklabels=1,
                                        dendrogram_ratio=(0.2, 0.2), mask=mask_array, row_colors=colors_df)

            # If number of groups <= 10, add legend
            if groups_num <= 10:
                legend_handles = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=group_to_color[group], markersize=8,
                           markeredgecolor='black', markeredgewidth=0.1, label=group) for group in unique_groups]

                # Add the legend to the plot
                clustermap.ax_heatmap.legend(handles=legend_handles, title=feature, loc="upper left",
                                             bbox_to_anchor=(1.06, 1.25), fontsize=10)

    # No metadata
    else:
        clustermap = sns.clustermap(matrix, metric=method, cmap=colmap, row_cluster=True, linewidths=.5,
                                    cbar_pos=(0.04, 0.82, 0.02, 0.15), xticklabels=1, yticklabels=1,
                                    dendrogram_ratio=(0.2, 0.2), mask=mask_array)

    if col_num <= 10:
        font_size = 14
    elif 10 < col_num <= 20:
        font_size = 12
    elif 20 < col_num <= 30:
        font_size = 11
    elif 30 < col_num <= 40:
        font_size = 10
    elif 40 < col_num <= 50:
        font_size = 9
    elif 50 < col_num <= 60:
        font_size = 8
    elif 60 < col_num <= 70:
        font_size = 7
    elif 70 < col_num <= 80:
        font_size = 6
    elif 80 < col_num <= 100:
        font_size = 5
    elif 100 < col_num <= 150:
        font_size = 4
    else:
        font_size = 3

    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=font_size,
                                          rotation='vertical')
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=font_size,
                                          rotation='horizontal')

    # Access the colorbar and add a label
    clustermap.cax.set_ylabel(type, labelpad=8)
    # Turn on the colorbar border
    clustermap.cax.spines['top'].set_visible(True)
    clustermap.cax.spines['right'].set_visible(True)
    clustermap.cax.spines['bottom'].set_visible(True)
    clustermap.cax.spines['left'].set_visible(True)

    # Set the font size of the feature label to be the same as the xticklabels
    if is_metadata:
        for ax in clustermap.fig.axes:
            if len(ax.get_xticklabels()) > 0:
                ax.set_xticklabels(ax.get_xticklabels(), fontsize=font_size)

    plt.close(clustermap.figure)

    return clustermap.figure


def cretae_network_plot(network, is_metadata, nodes_feature, is_continuous, cmap, custom_cmap, node_color, edge_color,
                        is_highlight_group, highlight_feature, highlight_group,
                        is_edge_colorby, edges_feature, within_edge_color, between_edge_color, iterations, pos_dict,
                        show_labels, all_or_highlighted, is_highlight_samples, samples_to_highlight, metadata_dict):
    iter_num = int(iterations)
    #print("\nIn cretae_network_plot.\nIterations number = " + str(iter_num) + "\nFeature: " + nodes_feature)

    pos = nx.layout.fruchterman_reingold_layout(network, iterations=iter_num, pos=pos_dict, k=2)

    is_legend = False

    nodes_to_label = []
    for node in network.nodes:
        # Set the default size and outline first
        network.nodes[node]['node_size'] = config.hvnx_nodes_size
        network.nodes[node]['outline_color'] = config.outline_color
        network.nodes[node]['outline_width'] = config.outline_width

        # If the 'highlight samples' checkbox is checked - set highlighted nodes differently
        if is_highlight_samples:
            highlighted_samples_list = re.split(r'\s*,\s*', samples_to_highlight)
            # Highlighted node
            if node in highlighted_samples_list:
                network.nodes[node]['node_size'] = config.hvnx_highlighted_nodes_size
                network.nodes[node]['outline_color'] = config.highlighted_outline_color
                network.nodes[node]['outline_width'] = config.highlighted_outline_width
                nodes_to_label.append(node)

        # In case the 'highlight nodes by feature' checkbox is checked,
        # highlight the nodes that belong to the selected group
        if is_metadata and is_highlight_group:
            if str(network.nodes[node][highlight_feature]) == highlight_group:
                network.nodes[node]['node_size'] = config.hvnx_highlighted_nodes_size
                network.nodes[node]['outline_color'] = config.highlighted_outline_color
                network.nodes[node]['outline_width'] = config.highlighted_outline_width
                nodes_to_label.append(node)

    if is_metadata:
        tooltips = [
            ('SampleID', '@SampleID'),
            (nodes_feature, '@' + nodes_feature)
        ]

        if is_edge_colorby and edges_feature != nodes_feature:
            tooltips.append((edges_feature, '@' + edges_feature))

        # Set the edges colors for the requested color-by feature
        edges = network.edges()
        for u, v in edges:
            # color the same/different edges according to the set colors
            if is_edge_colorby:
                if metadata_dict[edges_feature][u] == metadata_dict[edges_feature][v]:
                    network.edges[u, v]['edge_color'] = within_edge_color
                else:
                    network.edges[u, v]['edge_color'] = between_edge_color
            # All edges have the same color
            else:
                network.edges[u, v]['edge_color'] = edge_color
            network.edges[u, v]['edge_width'] = network.edges[u, v]['width']

        # In case of numeric continuous feature:
        if is_continuous:
            #print("Continuous feature")

            # Step 1: Extract non-missing values (only if the values are not strings)
            non_missing_values = [network.nodes[node][nodes_feature] for node in network.nodes
                                  if not np.isnan(network.nodes[node][nodes_feature])]

            # Step 2: Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            #print("Min value: " + str(min_value))
            #print("Max value: " + str(max_value))

            # Step 3: Create a colormap for non-missing values and normalize
            cmap = cmap
            norm = Normalize(vmin=min_value, vmax=max_value)

            network_plot = hvnx.draw(network, pos, node_color=nodes_feature, cmap=cmap, norm=norm,
                                     vmin=min_value, vmax=max_value)

            network_plot.opts(colorbar=True)

        # Feature is categorical
        else:
            #print("Categorical feature")

            # Prepare the colors mapping for the legend
            unique_groups = sorted(list(set([str(network.nodes[node][nodes_feature]) for node in network.nodes()])))
            groups_num = len(unique_groups)
            # Move the 'nan' group (if any) to the end of the list
            if 'nan' in unique_groups:
                unique_groups.remove('nan')
                unique_groups.append('nan')
            #print(unique_groups)

            cmap_length = len(cmap)

            # If the user defined a custom cmap - process it and turn it into a cmap
            if cmap_length == 1:
                custom_colors_list = custom_cmap
                cmap = re.split(r'\s*,\s*', custom_colors_list)
                cmap_length = len(cmap)
                #print("Custom cmap:")
                #print(cmap)

            group_to_color = {group: cmap[i % cmap_length] for i, group in enumerate(unique_groups)}
            colors = [group_to_color[str(network.nodes[node][nodes_feature])] for node in network.nodes()]

            network_plot = hvnx.draw(network, pos, node_color=colors)

            # Add a legend if there are up to max. groups allowed (value is defined in the config file)
            if groups_num <= config.max_groups_for_legend:
                legend_items = []
                for i, group in enumerate(unique_groups):
                    legend_items.append(hv.Scatter((0, 1 - i * 0.1), label=str(group)).opts(
                        show_legend=True, color=group_to_color[group], size=0))
                # Combine the plot with the generated legend
                legend = hv.Overlay(legend_items)
                is_legend = True

    # No metadata
    else:
        tooltips = [
            ('SampleID', '@SampleID')
        ]

        # Reset the edges colors
        edges = network.edges()
        for u, v in edges:
            network.edges[u, v]['edge_color'] = edge_color
            network.edges[u, v]['edge_width'] = network.edges[u, v]['width']

        network_plot = hvnx.draw(network, pos, node_color=node_color, node_alpha=0.95, edge_color=edge_color)

    hover = HoverTool(tooltips=tooltips)
    network_plot.opts(tools=[hover],
                      node_size='node_size',
                      node_line_color='outline_color',
                      node_line_width='outline_width',
                      node_alpha=0.95,
                      edge_color='edge_color',
                      edge_line_width='edge_width'
                      )

    if show_labels:
        # Show labels for all nodes
        if all_or_highlighted == 'All':
            labels = hv.Labels(network_plot.nodes, ['x', 'y'], 'index').opts(text_font_size='8pt', text_color='black')

            # Display legend and sample names
            if is_legend:
                hv_layout = network_plot * legend * labels
            else:
                hv_layout = network_plot * labels

        # Show labels for highlighted nodes only (if any)
        elif len(nodes_to_label) > 0:
            nodes_df = network_plot.nodes.data  # this is a pandas DataFrame
            subset_df = nodes_df[nodes_df['index'].isin(nodes_to_label)]
            labels = hv.Labels(subset_df, ['x', 'y'], 'index').opts(text_font_size='8pt', text_color='black')

            # Display legend and sample names
            if is_legend:
                hv_layout = network_plot * legend * labels
            else:
                hv_layout = network_plot * labels

        else:
            if is_legend:
                hv_layout = network_plot * legend
            else:
                hv_layout = network_plot

    else:
        if is_legend:
            hv_layout = network_plot * legend
        else:
            hv_layout = network_plot

    return hv_layout


def cretae_network_plot_matplotlib(network, is_metadata, nodes_feature, is_continuous, cmap, custom_cmap, node_color,
                                   edge_color, is_highlight_group, highlight_feature, highlight_group,
                                   is_edge_colorby, edges_feature, within_edge_color, between_edge_color,
                                   iterations, pos_dict, show_labels, all_or_highlighted, is_highlight_samples,
                                   samples_to_highlight, metadata_dict):
    before = time.time()
    iter_num = int(iterations)
    #print("\nIn cretae_network_plot_matplotlib. Iterations number = " + str(iter_num))
    #print("cmap: " + str(cmap))

    # Preparing network drawing also with matplotlib for better image production
    fig, ax1 = plt.subplots(figsize=(8, 7))
    widths = [network[u][v]['width'] for u, v in network.edges()]

    pos = nx.layout.fruchterman_reingold_layout(network, iterations=iter_num, pos=pos_dict, k=2)

    # Get the colormap
    if cmap != 'Define custom colormap':
        cmap_mpl = plt.get_cmap(cmap)

    outline_width_array = []
    outline_color_array = []
    nodes_size_array = []

    nodes_to_label = []
    for i, node in enumerate(network.nodes):
        # Set the default size and outline first
        outline_width_array.append(config.outline_width)
        outline_color_array.append(config.outline_color)
        nodes_size_array.append(config.nx_nodes_size)

        # If the 'highlight samples' checkbox is checked - set highlighted nodes differently
        if is_highlight_samples:
            highlighted_samples_list = re.split(r'\s*,\s*', samples_to_highlight)
            # Highlighted node
            if node in highlighted_samples_list:
                outline_width_array[i] = config.highlighted_outline_width_matplotlib
                outline_color_array[i] = config.highlighted_outline_color
                nodes_size_array[i] = config.nx_highlighted_nodes_size
                nodes_to_label.append(node)

        # In case the 'highlight nodes by feature' checkbox is checked,
        # highlight the nodes that belong to the selected group
        if is_metadata and is_highlight_group:
            if str(network.nodes[node][highlight_feature]) == highlight_group:
                outline_width_array[i] = config.highlighted_outline_width_matplotlib
                outline_color_array[i] = config.highlighted_outline_color
                nodes_size_array[i] = config.nx_highlighted_nodes_size
                nodes_to_label.append(node)

    if is_metadata:
        # Set the colors of the edges
        edges = network.edges()
        for u, v in edges:
            # Set the edges colors for the requested color-by feature
            if is_edge_colorby:
                if metadata_dict[edges_feature][u] == metadata_dict[edges_feature][v]:
                    network.edges[u, v]['edge_color'] = within_edge_color
                else:
                    network.edges[u, v]['edge_color'] = between_edge_color
            # All edges have the same color
            else:
                network.edges[u, v]['edge_color'] = edge_color

        edge_colors = nx.get_edge_attributes(network, 'edge_color').values()

        # In case of numeric continuous feature:
        if is_continuous:
            #print("Continuous feature")

            # Extract non-missing values
            non_missing_values = [network.nodes[node][nodes_feature] for node in network.nodes
                                  if not np.isnan(network.nodes[node][nodes_feature])]

            # Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            #print("Min value: " + str(min_value))
            #print("Max value: " + str(max_value))

            nodes_feature_array = np.array([network.nodes[node][nodes_feature] for node in network.nodes()])
            normalized_values = (nodes_feature_array - min_value) / (max_value - min_value)

            # Map the node values to colors using the colormap
            node_colors = []
            for i, node in enumerate(network.nodes):
                value = network.nodes[node][nodes_feature]
                if np.isnan(value):  # If the value is NaN, use the default color
                    node_colors.append(config.nodes_default_color)
                else:
                    node_colors.append(cmap_mpl(normalized_values[i]))  # Apply the colormap to the value

            # Create a ScalarMappable to associate the colormap with the normalized values
            norm = Normalize(vmin=min_value, vmax=max_value)
            sm = plt.cm.ScalarMappable(cmap=cmap_mpl, norm=norm)
            sm.set_array([])  # Empty array as we don't need to pass actual data to the ScalarMappable

            nx.draw(network, pos, node_size=nodes_size_array, node_color=node_colors, alpha=0.95,
                    edge_color=list(edge_colors), width=list(widths), with_labels=False,
                    vmin=min_value, vmax=max_value, linewidths=outline_width_array,
                    edgecolors=outline_color_array)

            plt.colorbar(ax=ax1, mappable=sm, label=nodes_feature)

        # Feature is categorical
        else:
            # Get a list of the unique groups of the selected feature
            unique_groups = sorted(list(set([str(network.nodes[node][nodes_feature]) for node in network.nodes()])))
            groups_num = len(unique_groups)
            # Move the 'nan' group (if any) to the end of the list
            if 'nan' in unique_groups:
                unique_groups.remove('nan')
                unique_groups.append('nan')
            #print(unique_groups)

            # If the user defined a custom cmap - process it and turn it into a cmap
            if cmap == 'Define custom colormap':
                custom_colors_list = custom_cmap
                cmap_mpl = re.split(r'\s*,\s*', custom_colors_list)
                cmap_length = len(cmap_mpl)
                #print("Custom cmap:")
                #print(cmap_mpl)

                # Assign each group with a color, according to the colormap order
                group_to_color = {group: cmap_mpl[i % cmap_length] for i, group in enumerate(unique_groups)}

            # Standard colormap
            else:
                cmap_length = cmap_mpl.N
                # Assign each group with a color, according to the colormap order
                group_to_color = {group: cmap_mpl(i % cmap_length) for i, group in enumerate(unique_groups)}

            #print("cmap length: " + str(cmap_length))

            # Get a list of colors for all the nodes
            colors = [group_to_color[str(network.nodes[node][nodes_feature])] for node in network.nodes()]

            nx.draw(network, pos, node_size=nodes_size_array, node_color=colors, alpha=0.95,
                    edge_color=list(edge_colors), width=list(widths), with_labels=False,
                    linewidths=outline_width_array, edgecolors=outline_color_array)

            # Add a legend if there are up to max. groups allowed (value is defined in the config file)
            if groups_num <= config.max_groups_for_legend:
                legend_handles = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=group_to_color[group], markersize=8,
                           markeredgecolor='black', markeredgewidth=0.1, label=group)
                    for group in unique_groups]

                # Add the legend to the plot
                plt.legend(handles=legend_handles, title=nodes_feature, loc="best")

    # No metadata coloring
    else:
        nx.draw(network, pos, node_size=nodes_size_array, node_color=node_color, alpha=0.95,
                edge_color=edge_color, width=list(widths), with_labels=False, linewidths=outline_width_array,
                edgecolors=outline_color_array)

    # Add the labels to the desired nodes
    if show_labels:
        # Show labels for all nodes
        if all_or_highlighted == 'All':
            labels_dict = {n: n for n in network.nodes}
        # Show labels for highlighted nodes only (if any)
        elif len(nodes_to_label) > 0:
            labels_dict = {n: n for n in nodes_to_label}
        else:
            labels_dict = {}

        nx.draw_networkx_labels(
            network, pos,
            labels=labels_dict,
            font_size=5, font_weight='bold'
        )

    plt.close(fig)

    after = time.time()
    duration = after - before
    print("\ncretae_network_plot_matplotlib: saving the network took " + str(duration) + " seconds.\n")

    return fig



