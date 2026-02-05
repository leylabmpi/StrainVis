import numpy as np
import pandas as pd
import hvplot.pandas  # Enable interactive
from bokeh.models import HoverTool
import seaborn as sns
import matplotlib.pyplot as plt
import StrainVis_app.config as config


def plot_pairs_vs_sampling_size_bar(df, sampling_size, is_all_regions):
    value = str(sampling_size)

    # If 'All regions' shouldn't be presented, remove this row from the df
    if is_all_regions == 0:
        df = df.drop(0)

    df['color'] = np.where(df['Subsampled_regions'] == value, config.highlight_bar_color, config.normal_bar_color)

    tooltips = [
        ('Number of pairs', '@Number_of_pairs'),
    ]
    hover = HoverTool(tooltips=tooltips)

    bar_plot = df.hvplot.bar(x='Subsampled_regions', y='Number_of_pairs', color='color', xlabel='Subsampled regions',
                             ylabel='Number of pairs').opts(shared_axes=False, tools=[hover])

    return bar_plot


def plot_species_vs_sampling_size_bar(df, sampling_size, is_all_regions):
    value = str(sampling_size)

    # If 'All regions' shouldn't be presented, remove this row from the df
    if is_all_regions == 0:
        df = df.drop(0)

    df['color'] = np.where(df['Subsampled_regions'] == value, config.highlight_bar_color, config.normal_bar_color)

    tooltips = [
        ('Number of species', '@Number_of_species'),
    ]
    hover = HoverTool(tooltips=tooltips)

    bar_plot = df.hvplot.bar(x='Subsampled_regions', y='Number_of_species', color='color', xlabel='Subsampled regions',
                             ylabel='Number of species').opts(shared_axes=False, tools=[hover])

    return bar_plot


def create_box_plot(avg_df, sorted_genomes_list, pvalues_df, color, use_metadata, feature, same_color, different_color):

    #print("\ncreate_box_plot: Feature is " + feature)
    #print("Genomes list:")
    #print(genomes_list)
    #print("\nAPSS dataframe:")
    #print(avg_df)

    presented_genomes_list = avg_df['Ref_genome'].unique()
    genomes_num = len(presented_genomes_list)

    if genomes_num <= 8:
        fig_height = 4
    else:
        fig_height = 0.5 * genomes_num
    plt.figure(figsize=(7.2, fig_height))

    # Use metadata to separate plot to same/different feature
    if use_metadata:

        same_feature = 'Same ' + feature
        diff_feature = 'Different ' + feature

        box_plot = sns.boxplot(data=avg_df, x="APSS", y="Ref_genome", order=sorted_genomes_list,
                               hue="Category", width=0.8, showcaps=False, gap=0.1,
                               hue_order=[same_feature, diff_feature],
                               palette=[same_color, different_color]
        )

        # Get the x,y axes positions (corresponding to the groups)
        y_ticks = box_plot.get_yticks()
        x_ticks = box_plot.get_xticks()
        space_x = x_ticks[1] - x_ticks[0]
        space_y = y_ticks[1] - y_ticks[0]
        x_range = box_plot.get_xlim()
        y_range = box_plot.get_ylim()
        x_upper_limit = x_range[1]
        x_pos_for_stars = x_upper_limit + space_x * 0.2
        x_pos_for_effect_size = x_pos_for_stars + space_x * 0.7
        x_pos_for_title = x_upper_limit - space_x * 0.2
        y_pos_for_title = y_range[1] - space_y * 0.4
        is_significance = 0
        is_effect_size = 0
        for i, row in pvalues_df.iterrows():
            y_pos = y_ticks[i]  # This gives the y-position of the i-th group
            if row['Significance'] != "NS":
                is_significance = 1
                y_pos_for_stars = y_pos + space_y * 0.06
                box_plot.text(x=x_pos_for_stars, y=y_pos_for_stars, s=row['Significance'], fontsize=14)
                if pd.notna(row['Effect_size']):
                    is_effect_size = 1
                    box_plot.text(x=x_pos_for_effect_size, y=y_pos, s=row['Effect_size'], fontsize=11)
        # If there is at least one row with significance / effect_size calculation, print title for them
        if is_significance:
            if is_effect_size:
                box_plot.text(x=x_pos_for_title, y=y_pos_for_title, s='Significance  Effect-size', fontsize=10)
            else:
                box_plot.text(x=x_pos_for_title, y=y_pos_for_title, s='Significance', fontsize=11)

        box_plot.legend(title="", loc='lower left', bbox_to_anchor=(0, 1), fontsize=13)

    # Do not use metadata in plot - show all the comparisons together
    else:
        box_plot = sns.boxplot(data=avg_df, x="APSS", y="Ref_genome", order=sorted_genomes_list, color=color, width=0.5,
                               showcaps=False)

    box_plot.yaxis.grid(True)  # Hide the horizontal gridlines
    box_plot.xaxis.grid(True)  # Show the vertical gridlines
    box_plot.set_xlabel('APSS', fontsize=14, labelpad=7)
    box_plot.set_ylabel('Species', fontsize=14, labelpad=7)

    fig = box_plot.figure
    plt.close()

    return fig


def create_box_plot_ani(ani_df, sorted_genomes_list, pvalues_df, color, use_metadata, feature, same_color,
                        different_color):

    #print("\ncreate_box_plot_ani: Feature is " + feature)
    #print("Genomes list:")
    #print(genomes_list)
    #print("\nANI dataframe:")
    #print(ani_df)

    genomes_num = len(sorted_genomes_list)
    if genomes_num <= 8:
        fig_height = 4
    else:
        fig_height = 0.5 * genomes_num
    plt.figure(figsize=(7.2, fig_height))

    # Use metadata to separate plot to same/different feature
    if use_metadata:

        same_feature = 'Same ' + feature
        diff_feature = 'Different ' + feature

        box_plot = sns.boxplot(data=ani_df, x="ANI", y="Ref_genome", order=sorted_genomes_list,
                               hue="Category", width=0.8, showcaps=False, gap=0.1,
                               hue_order=[same_feature, diff_feature],
                               palette=[same_color, different_color]
        )

        # Get the x,y axes positions (corresponding to the groups)
        y_ticks = box_plot.get_yticks()
        x_ticks = box_plot.get_xticks()
        space_x = x_ticks[1] - x_ticks[0]
        space_y = y_ticks[1] - y_ticks[0]
        x_range = box_plot.get_xlim()
        y_range = box_plot.get_ylim()
        x_upper_limit = x_range[1]
        x_pos_for_stars = x_upper_limit + space_x * 0.2
        x_pos_for_effect_size = x_pos_for_stars + space_x * 0.9
        x_pos_for_title = x_upper_limit - space_x * 0.2
        y_pos_for_title = y_range[1] - space_y * 0.4
        is_significance = 0
        is_effect_size = 0

        for i, row in pvalues_df.iterrows():
            y_pos = y_ticks[i]  # This gives the y-position of the i-th group
            if row['Significance'] != "NS":
                is_significance = 1
                y_pos_for_stars = y_pos + space_y * 0.06
                box_plot.text(x=x_pos_for_stars, y=y_pos_for_stars, s=row['Significance'], fontsize=14)
                if pd.notna(row['Effect_size']):
                    is_effect_size = 1
                    box_plot.text(x=x_pos_for_effect_size, y=y_pos, s=row['Effect_size'], fontsize=11)
        # If there is at least one row with significance / effect_size calculation, print title for them
        if is_significance:
            if is_effect_size:
                box_plot.text(x=x_pos_for_title, y=y_pos_for_title, s='Significance  Effect-size', fontsize=10)
            else:
                box_plot.text(x=x_pos_for_title, y=y_pos_for_title, s='Significance', fontsize=11)

        box_plot.legend(title="", loc='lower left', bbox_to_anchor=(0, 1), fontsize=13)

    # Do not use metadata in plot - show all the comparisons together
    else:
        box_plot = sns.boxplot(data=ani_df, x="ANI", y="Ref_genome", order=sorted_genomes_list, color=color,
                               width=0.5, showcaps=False)

    box_plot.yaxis.grid(True)  # Hide the horizontal gridlines
    box_plot.xaxis.grid(True)  # Show the vertical gridlines
    box_plot.set_xlabel('ANI', fontsize=14, labelpad=7)
    box_plot.set_ylabel('Species', fontsize=14, labelpad=7)

    fig = box_plot.figure
    plt.close()

    return fig
