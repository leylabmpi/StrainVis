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
import SynTrackerVis_app.config as config
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


def create_clustermap(matrix, cmap, is_metadata, feature, cmap_metadata, custom_cmap, metadata_dict):

    # The number of columns doesn't exceed the defined maximum - continue creating the clustermap plot
    col_num = len(matrix.columns)
    mask_array = np.full((col_num, col_num), np.where(matrix == 100, True, False))

    print("\ncreate_clustermap: selected cmap: " + cmap)
    colmap = plt.get_cmap(cmap)
    colmap.set_bad("lightgrey")

    if is_metadata:

        # Get the unique groups of the selected feature
        unique_groups = sorted(list(set([str(metadata_dict[feature][sample]) for sample in matrix.iloc[:, 0].index])))
        groups_num = len(unique_groups)

        # Move the 'nan' group (if any) to the end of the list
        if 'nan' in unique_groups:
            unique_groups.remove('nan')
            unique_groups.append('nan')
        print(unique_groups)

        # If the user defined a custom cmap - process it and turn it into a cmap
        if cmap_metadata == 'Define custom colormap':
            custom_colors_list = custom_cmap
            cmap_metadata_mpl = re.split(r'\s*,\s*', custom_colors_list)
            cmap_length = len(cmap_metadata_mpl)
            print("Custom cmap:")
            print(cmap_metadata_mpl)

            # Assign each group with a color, according to the colormap order
            group_to_color = {group: cmap_metadata_mpl[i % cmap_length] for i, group in enumerate(unique_groups)}

        # Standard colormap
        else:
            cmap_metadata_mpl = plt.get_cmap(cmap_metadata)
            cmap_length = cmap_metadata_mpl.N

            # Assign each group with a color, according to the colormap order
            group_to_color = {group: cmap_metadata_mpl(i % cmap_length) for i, group in enumerate(unique_groups)}

        print("cmap length: " + str(cmap_length))

        # Create a colors dataframe, with sample names as indices
        colors = [group_to_color[str(metadata_dict[feature][sample])] for sample in matrix.iloc[:, 0].index]
        colors_df = pd.DataFrame({feature: colors}, index=matrix.index)
        print(colors_df)

        clustermap = sns.clustermap(matrix, cmap=colmap, row_cluster=True, linewidths=.5,
                                    cbar_pos=(0.04, 0.82, 0.02, 0.15), xticklabels=1, yticklabels=1,
                                    dendrogram_ratio=(0.2, 0.2), mask=mask_array, row_colors=colors_df)

    # No metadata
    else:
        clustermap = sns.clustermap(matrix, cmap=colmap, row_cluster=True, linewidths=.5,
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
    else:
        font_size = 4

    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=font_size,
                                          rotation='vertical')
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=font_size,
                                          rotation='horizontal')

    # Set the font size of the feature label to be the same as the xticklabels
    if is_metadata:
        for ax in clustermap.fig.axes:
            if len(ax.get_xticklabels()) > 0:
                ax.set_xticklabels(ax.get_xticklabels(), fontsize=font_size)

        # If number of groups <= 10, add legend
        if groups_num <= 10:
            legend_handles = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor=group_to_color[group], markersize=8,
                       markeredgecolor='black', markeredgewidth=0.1, label=group) for group in unique_groups]

            # Add the legend to the plot
            clustermap.ax_heatmap.legend(handles=legend_handles, title=feature, loc="upper left",
                                         bbox_to_anchor=(1.06, 1.25), fontsize=10)

    plt.close(clustermap.figure)

    return clustermap.figure


def cretae_network_plot(network, is_metadata, nodes_feature, is_continuous, cmap, custom_cmap, node_color, edge_color,
                        is_edge_colorby, edges_feature, within_edge_color, between_edge_color, iterations, pos_dict,
                        show_labels, metadata_dict):
    iter_num = int(iterations)
    print("\nIn cretae_network_plot.\nIterations number = " + str(iter_num) + "\nFeature: " + nodes_feature)

    pos = nx.layout.fruchterman_reingold_layout(network, iterations=iter_num, pos=pos_dict, k=2)

    is_legend = False

    if is_metadata:
        tooltips = [
            ('SampleID', '@SampleID'),
            (nodes_feature, '@' + nodes_feature)
        ]

        # Set the edges colors for the requested color-by feature
        if is_edge_colorby:
            edges = network.edges()
            for u, v in edges:
                if metadata_dict[edges_feature][u] == metadata_dict[edges_feature][v]:
                    network.edges[u, v]['edge_color'] = within_edge_color
                else:
                    network.edges[u, v]['edge_color'] = between_edge_color

            tooltips.append((edges_feature, '@' + edges_feature))

        # In case of numeric continuous feature:
        if is_continuous:
            print("Continuous feature")

            # Step 1: Extract non-missing values (only if the values are not strings)
            non_missing_values = [network.nodes[node][nodes_feature] for node in network.nodes
                                  if not np.isnan(network.nodes[node][nodes_feature])]

            # Step 2: Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            print("Min value: " + str(min_value))
            print("Max value: " + str(max_value))

            # Step 3: Create a colormap for non-missing values and normalize
            cmap = cmap
            norm = Normalize(vmin=min_value, vmax=max_value)

            if is_edge_colorby:
                # Plot using holoviews
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, norm=norm,
                                         node_alpha=0.95, edge_color=hv.dim('edge_color'),
                                         edge_width=hv.dim('weight') / 5, vmin=min_value, vmax=max_value)

            else:
                # Plot using holoviews
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, norm=norm,
                                         node_alpha=0.95, edge_color=edge_color, edge_width=hv.dim('weight')/5,
                                         vmin=min_value, vmax=max_value)

            network_plot.opts(colorbar=True)

        # Feature is categorical
        else:
            print("Categorical feature")

            # Prepare the colors mapping for the legend
            unique_groups = sorted(list(set([str(network.nodes[node][nodes_feature]) for node in network.nodes()])))
            groups_num = len(unique_groups)
            # Move the 'nan' group (if any) to the end of the list
            if 'nan' in unique_groups:
                unique_groups.remove('nan')
                unique_groups.append('nan')
            print(unique_groups)

            cmap_length = len(cmap)
            print("Cmap length = " + str(cmap_length))

            # If the user defined a custom cmap - process it and turn it into a cmap
            if cmap_length == 1:
                custom_colors_list = custom_cmap
                cmap = re.split(r'\s*,\s*', custom_colors_list)
                cmap_length = len(cmap)
                print("Custom cmap:")
                print(cmap)
                print("custom cmap length: " + str(cmap_length))

            group_to_color = {group: cmap[i % cmap_length] for i, group in enumerate(unique_groups)}
            colors = [group_to_color[str(network.nodes[node][nodes_feature])] for node in network.nodes()]

            if is_edge_colorby:
                # Plot using holoviews
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=colors, cmap=cmap,
                                         node_alpha=0.95,
                                         edge_color=hv.dim('edge_color'), edge_width=hv.dim('weight') / 5)

            else:
                # Plot using holoviews
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=colors, node_alpha=0.95,
                                         edge_color=edge_color, edge_width=hv.dim('weight')/5)

            # Add a legend if there are up to 10 groups
            if groups_num <= 10:
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

        # Plot using holoviews
        network_plot = hvnx.draw(network, pos, node_size=300, node_color=node_color, node_alpha=0.95,
                                 edge_color=edge_color, edge_width=hv.dim('weight')/5)

    hover = HoverTool(tooltips=tooltips)
    network_plot.opts(tools=[hover])

    if show_labels:
        labels = hv.Labels(network_plot.nodes, ['x', 'y'], 'index').opts(text_font_size='8pt', text_color='black')

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

    return hv_layout


def cretae_network_plot_matplotlib(network, is_metadata, nodes_feature, is_continuous, cmap, custom_cmap, node_color,
                                   edge_color, is_edge_colorby, edges_feature, within_edge_color, between_edge_color,
                                   iterations, pos_dict, show_labels, metadata_dict):
    iter_num = int(iterations)
    print("\nIn cretae_network_plot_matplotlib. Iterations number = " + str(iter_num))
    print("cmap: " + str(cmap))

    # Preparing network drawing also with matplotlib for better image production
    fig, ax1 = plt.subplots(figsize=(8, 7))
    widths = [network[u][v]['weight']/10 for u, v in network.edges()]

    pos = nx.layout.fruchterman_reingold_layout(network, iterations=iter_num, pos=pos_dict, k=2)

    # Get the colormap
    if cmap != 'Define custom colormap':
        cmap_mpl = plt.get_cmap(cmap)

    if is_metadata:
        # Set the edges colors for the requested color-by feature
        if is_edge_colorby:
            edges = network.edges()
            for u, v in edges:
                if metadata_dict[edges_feature][u] == metadata_dict[edges_feature][v]:
                    network.edges[u, v]['edge_color'] = within_edge_color
                else:
                    network.edges[u, v]['edge_color'] = between_edge_color

            edge_colors = nx.get_edge_attributes(network, 'edge_color').values()

        # In case of numeric continuous feature:
        if is_continuous:
            print("Continuous feature")

            # Extract non-missing values
            non_missing_values = [network.nodes[node][nodes_feature] for node in network.nodes
                                  if not np.isnan(network.nodes[node][nodes_feature])]

            # Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            print("Min value: " + str(min_value))
            print("Max value: " + str(max_value))

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

            if is_edge_colorby:
                if show_labels:
                    nx.draw(network, pos, node_size=80, node_color=node_colors, alpha=0.95,
                            edge_color=list(edge_colors), width=list(widths), vmin=min_value, vmax=max_value,
                            with_labels=True, linewidths=0.5, edgecolors='black', font_size=6)
                else:
                    nx.draw(network, pos, node_size=80, node_color=node_colors, alpha=0.95,
                            edge_color=list(edge_colors), width=list(widths), with_labels=False,
                            vmin=min_value, vmax=max_value, linewidths=0.5, edgecolors='black')

            else:
                if show_labels:
                    nx.draw(network, pos, node_size=80, node_color=node_colors, alpha=0.95,
                            edge_color=edge_color, width=list(widths), vmin=min_value, vmax=max_value, with_labels=True,
                            linewidths=0.5, edgecolors='black', font_size=6)
                else:
                    nx.draw(network, pos, node_size=80, node_color=node_colors, alpha=0.95,
                            edge_color=edge_color, width=list(widths), with_labels=False,
                            vmin=min_value, vmax=max_value, linewidths=0.5, edgecolors='black')

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
            print(unique_groups)

            # If the user defined a custom cmap - process it and turn it into a cmap
            if cmap == 'Define custom colormap':
                custom_colors_list = custom_cmap
                cmap_mpl = re.split(r'\s*,\s*', custom_colors_list)
                cmap_length = len(cmap_mpl)
                print("Custom cmap:")
                print(cmap_mpl)

                # Assign each group with a color, according to the colormap order
                group_to_color = {group: cmap_mpl[i % cmap_length] for i, group in enumerate(unique_groups)}

            # Standard colormap
            else:
                cmap_length = cmap_mpl.N
                # Assign each group with a color, according to the colormap order
                group_to_color = {group: cmap_mpl(i % cmap_length) for i, group in enumerate(unique_groups)}

            print("cmap length: " + str(cmap_length))

            # Get a list of colors for all the nodes
            colors = [group_to_color[str(network.nodes[node][nodes_feature])] for node in network.nodes()]

            if is_edge_colorby:
                if show_labels:
                    nx.draw(network, pos, node_size=80, node_color=colors, alpha=0.95,
                            edge_color=list(edge_colors), width=list(widths), with_labels=True, linewidths=0.5,
                            edgecolors='black', font_size=6)
                else:
                    nx.draw(network, pos, node_size=80, node_color=colors, alpha=0.95,
                            edge_color=list(edge_colors), width=list(widths), with_labels=False, linewidths=0.5,
                            edgecolors='black')

            else:
                if show_labels:
                    nx.draw(network, pos, node_size=80, node_color=colors, alpha=0.95,
                            edge_color=edge_color, width=list(widths), with_labels=True, linewidths=0.5,
                            edgecolors='black', font_size=6)
                else:
                    nx.draw(network, pos, node_size=80, node_color=colors, alpha=0.95,
                            edge_color=edge_color, width=list(widths), with_labels=False, linewidths=0.5,
                            edgecolors='black')

            # If number of groups <= 10, add legend
            if groups_num <= 10:
                legend_handles = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=group_to_color[group], markersize=8,
                           markeredgecolor='black', markeredgewidth=0.1, label=group)
                    for group in unique_groups]

                # Add the legend to the plot
                plt.legend(handles=legend_handles, title=nodes_feature, loc="best")

    else:
        if show_labels:
            nx.draw(network, pos, node_size=100, node_color=node_color, alpha=0.95, edge_color=edge_color,
                    width=list(widths), with_labels=True, linewidths=0.5, edgecolors='black', font_size=6)
        else:
            nx.draw(network, pos, node_size=100, node_color=node_color, alpha=0.95, edge_color=edge_color,
                    width=list(widths), with_labels=False, linewidths=0.5, edgecolors='black')

    plt.close(fig)

    return fig


def create_coverage_plot(contig_name, score_per_pos_contig, avg_score_per_pos_contig,
                         start_pos, end_pos,
                         show_avg, avg_color, show_scores, scores_color,
                         show_hyper_var, hyper_var_color, hyper_var_alpha,
                         show_hyper_cons, hyper_cons_color, hyper_cons_alpha, bottom_val):

    before = time.time()
    print("\ncreate_coverage_plot:\nContig name: " + contig_name)
    print("Start position: " + start_pos)
    print("End position: " + end_pos)

    # Consider only the positions within the requested range
    score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] >= int(start_pos)]
    score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] < int(end_pos)]
    avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] >= int(start_pos)]
    avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

    #print("\nScore per position table:")
    #print(score_per_pos_contig)
    #print("\nAVG score per position table:")
    #print(avg_score_per_pos_contig)

    after = time.time()
    duration = after - before
    print("Setting the correct range took " + str(duration) + " seconds")

    # Prepare data for plotting the avg scores as lines
    avg_score_per_pos_contig_end_pos = avg_score_per_pos_contig.copy()
    avg_score_per_pos_contig_end_pos['Position'] = avg_score_per_pos_contig['Position'] + config.region_length - 100
    #print(avg_score_per_pos_contig_end_pos)
    #print("\nAfter creating the second df:")
    #print(avg_score_per_pos_contig)

    avg_score_per_pos_contig_for_line_plot = pd.concat([avg_score_per_pos_contig, avg_score_per_pos_contig_end_pos],
                                                       ignore_index=True).sort_values(by='Position')

    #print("\nAverage data for line plot after concatenating:")
    #print(avg_score_per_pos_contig_for_line_plot)

    pos_array = np.full((2, len(score_per_pos_contig.index)), 0)
    pos_array[1, :] = config.region_length - 100

    avg_pos_array = np.full((2, len(avg_score_per_pos_contig.index)), 0)
    avg_pos_array[1, :] = config.region_length

    fig, ax1 = plt.subplots(figsize=(11, 5))

    if show_scores:
        coverage_plot = plt.errorbar(score_per_pos_contig['Position'], score_per_pos_contig['Synteny_score'],
                                     xerr=pos_array, color=scores_color, fmt='none', elinewidth=0.7, zorder=1,
                                     label='Synteny scores')

    if show_avg:
        line_avg = plt.plot(avg_score_per_pos_contig_for_line_plot['Position'],
                            avg_score_per_pos_contig_for_line_plot['Avg_synteny_score'],
                            color=avg_color, zorder=2, label='Average synteny scores')

    plt.xticks(avg_score_per_pos_contig['Position'], fontsize=6, rotation=90)
    ax1.locator_params(axis='x', tight=True, nbins=40)
    plt.xlabel("Position in reference genome/contig", labelpad=8)

    if show_scores:
        plt.ylabel("Synteny Score")
    else:
        plt.ylabel("Average synteny Score")

    if show_hyper_var:
        hypervar_bars = plt.bar(avg_score_per_pos_contig['Position'],
                        avg_score_per_pos_contig['Hypervariable'], align='edge',
                        width=config.region_length, bottom=bottom_val, color=hyper_var_color, alpha=hyper_var_alpha,
                        label='Hypervariable regions')

    if show_hyper_cons:
        hypercons_bars = plt.bar(avg_score_per_pos_contig['Position'],
                         avg_score_per_pos_contig['Hyperconserved'], align='edge',
                         width=config.region_length, bottom=bottom_val, color=hyper_cons_color, alpha=hyper_cons_alpha,
                         label='Hyperconserved regions')

    plt.legend(fontsize='small', loc=(0, 1.02))

    plt.close(fig)

    after = time.time()
    duration = after - before
    print("Create/update the coverage plot took " + str(duration) + " seconds")

    return fig



