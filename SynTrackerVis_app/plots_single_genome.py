import numpy as np
import pandas as pd
import networkx as nx
import hvplot.pandas  # Enable interactive
import holoviews as hv
import hvplot.networkx as hvnx
import matplotlib as mpl
import matplotlib.pyplot as plt
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

    # DEbug
    #export_png(plot, filename="/Users/ipaz/ownCloud/Projects/SynTracker_Vis/Downloads/Jitter_try.png")
    ###

    return plot


def jitter_dots(dots):
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.1, 0.1, offsets.shape[0])
    dots.set_offsets(jittered_offsets)


def category_by_feature(row, feature, metadata_dict):
    if metadata_dict[feature][row['Sample1']] == metadata_dict[feature][row['Sample2']]:
        return 'Same ' + feature
    else:
        return 'Different ' + feature


def create_jitter_plot(avg_df, color, use_metadata, metadata_dict, feature, same_color, different_color):
    df_for_jitter = avg_df.loc[:, ['Sample1', 'Sample2', 'APSS']].copy()

    # Use metadata to separate plot to same/different feature
    if use_metadata:
        df_for_jitter['Category'] = df_for_jitter.apply(lambda row: category_by_feature(row, feature, metadata_dict),
                                                        axis=1)
        same_feature = 'Same ' + feature
        diff_feature = 'Different ' + feature
        jitter_plot = sns.catplot(data=df_for_jitter, x="Category", y="APSS", order=[same_feature, diff_feature],
                                  hue="Category", hue_order=[same_feature, diff_feature],
                                  palette=[same_color, different_color], edgecolor="gray", linewidth=0.1)

    # Do not use metadata in plot - show all the comparisons together
    else:
        df_for_jitter['Category'] = 'All Comparisons'
        jitter_plot = sns.catplot(data=df_for_jitter, x="Category", y="APSS", color=color, edgecolor="gray",
                                  linewidth=0.1)

    print("\nDF for jitter plot:")
    print(df_for_jitter)

    return jitter_plot.figure


def create_clustermap(matrix, cmap):

    np.fill_diagonal(matrix.values, 1.0)
    matrix = matrix.fillna(100.0)

    # The number of columns doesn't exceed the defined maximum - continue creating the clustermap plot
    col_num = len(matrix.columns)
    mask_array = np.full((col_num, col_num), np.where(matrix == 100, True, False))

    print("\ncreate_clustermap: selected cmap: " + cmap)
    colmap = mpl.colormaps.get_cmap(cmap)
    colmap.set_bad("lightgrey")

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

    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xmajorticklabels(), fontsize=font_size,
                                          rotation='vertical')
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_ymajorticklabels(), fontsize=font_size,
                                          rotation='horizontal')

    return clustermap.figure


def cretae_network_plot(network, is_metadata, nodes_feature, is_continuous, cmap, node_color, edge_color,
                        is_edge_colorby, edges_feature, within_edge_color, between_edge_color, iterations, pos_dict,
                        show_labels, metadata_dict):
    iter_num = int(iterations)

    pos = nx.layout.fruchterman_reingold_layout(network, iterations=iter_num, pos=pos_dict, k=2)

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

            # Step 1: Extract non-missing values
            non_missing_values = [network.nodes[node][nodes_feature] for node in network.nodes
                                  if network.nodes[node][nodes_feature] is not np.nan]

            # Step 2: Define min and max values based only on non-missing values
            min_value = min(non_missing_values)
            max_value = max(non_missing_values)
            print("Min value: " + str(min_value))
            print("Max value: " + str(max_value))

            # Step 3: Create a colormap for non-missing values and normalize
            cmap = cmap
            norm = Normalize(vmin=min_value, vmax=max_value)

            if is_edge_colorby:
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, norm=norm,
                                         node_alpha=0.95, edge_color=hv.dim('edge_color'),
                                         edge_width=hv.dim('weight') / 5, vmin=min_value, vmax=max_value)

            else:
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, norm=norm,
                                         node_alpha=0.95, edge_color=edge_color, edge_width=hv.dim('weight')/5,
                                         vmin=min_value, vmax=max_value)

            network_plot.opts(colorbar=True)

        # Feature is categorical
        else:
            if is_edge_colorby:
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, node_alpha=0.95,
                                         edge_color=hv.dim('edge_color'), edge_width=hv.dim('weight') / 5)
            else:
                network_plot = hvnx.draw(network, pos, node_size=300, node_color=nodes_feature, cmap=cmap, node_alpha=0.95,
                                         edge_color=edge_color, edge_width=hv.dim('weight')/5)

        hover = HoverTool(tooltips=tooltips)
        network_plot.opts(tools=[hover])

    else:
        network_plot = hvnx.draw(network, pos, node_size=300, node_color=node_color, node_alpha=0.95,
                                 edge_color=edge_color, edge_width=hv.dim('weight')/5)

    if show_labels:
        labels = hv.Labels(network_plot.nodes, ['x', 'y'], 'index').opts(text_font_size='8pt', text_color='black')
        hv_layout = network_plot * labels
    else:
        hv_layout = network_plot

    return hv_layout


def create_coverage_plot(contig_name, score_per_pos_contig, avg_score_per_pos_contig,
                         start_pos, end_pos,
                         show_avg, avg_color, show_scores, scores_color,
                         show_hyper_var, hyper_var_color, hyper_var_alpha,
                         show_hyper_cons, hyper_cons_color, hyper_cons_alpha, bottom_val):

    print("\ncreate_coverage_plot:\nContig name: " + contig_name)
    print("Start position: " + start_pos)
    print("End position: " + end_pos)

    #print("\nScore per contig df sorted:")
    #print(score_per_pos_contig)

    # Consider only the positions within the requested range
    score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] >= int(start_pos)]
    score_per_pos_contig = score_per_pos_contig[score_per_pos_contig['Position'] < int(end_pos)]
    avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] >= int(start_pos)]
    avg_score_per_pos_contig = avg_score_per_pos_contig[avg_score_per_pos_contig['Position'] < int(end_pos)]

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

    fig, ax1 = plt.subplots(figsize=(10, 5))

    pos_array = np.full((2, len(score_per_pos_contig.index)), 0)
    pos_array[1, :] = config.region_length - 100

    avg_pos_array = np.full((2, len(avg_score_per_pos_contig.index)), 0)
    avg_pos_array[1, :] = config.region_length

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
        plt.bar(avg_score_per_pos_contig['Position'],
                avg_score_per_pos_contig['Hypervariable'], align='edge',
                width=config.region_length, bottom=bottom_val, color=hyper_var_color, alpha=hyper_var_alpha,
                label='Hypervariable regions')

    if show_hyper_cons:
        plt.bar(avg_score_per_pos_contig['Position'],
                avg_score_per_pos_contig['Hyperconserved'], align='edge',
                width=config.region_length, bottom=bottom_val, color=hyper_cons_color, alpha=hyper_cons_alpha,
                label='Hyperconserved regions')

    plt.legend(fontsize='small', loc=(0, 1.02))
    #plt.legend(fontsize='small', loc=(1.01, 0.8))

    plt.close(fig)

    return fig



