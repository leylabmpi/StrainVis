import numpy as np
import pandas as pd
import hvplot.pandas  # Enable interactive
from bokeh.models import HoverTool
import SynTrackerVis_app.config as config


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

