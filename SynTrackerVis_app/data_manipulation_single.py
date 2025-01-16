import pandas as pd
import numpy as np
import re
import SynTrackerVis_app.config as config


def return_selected_genome_table(score_per_region_df, selected_genome):
    selected_genome_df = score_per_region_df[score_per_region_df['Ref_genome'] == selected_genome]
    selected_genome_df = selected_genome_df[['Sample1', 'Sample2', 'Region', 'Synteny_score']]

    #print(selected_genome_df)
    return selected_genome_df


def return_selected_genome_avg_table(avg_big_df, selected_genome):
    selected_genome_avg_df = avg_big_df[avg_big_df['Ref_genome'] == selected_genome]
    selected_genome_avg_df = selected_genome_avg_df[['Sample1', 'Sample2', 'Avg_synteny_score', 'Compared_regions']]

    #print(selected_genome_avg_df)
    return selected_genome_avg_df


def calculate_avg_scores_selected_genome_size(score_per_region_selected_genome_df, genome, size):

    # Taking all available regions - no subsampling
    if size == 'All_regions':
        avg_scores_one_size_df = score_per_region_selected_genome_df.groupby(['Sample1', 'Sample2'])['Synteny_score'].\
            mean().reset_index().rename(columns={"Synteny_score": "Avg_synteny_score"})

    # All the other optional sizes -
    # need to include only the sample pairs that have a result for at least the requested number of regions
    else:
        filtered_df = score_per_region_selected_genome_df[['Sample1', 'Sample2', 'Synteny_score']].\
            groupby(['Sample1', 'Sample2']).filter(lambda x: x['Synteny_score'].count() >= int(size))
        #print(filtered_df)

        sampled_regions_df = filtered_df[['Sample1', 'Sample2', 'Synteny_score']].\
            groupby(['Sample1', 'Sample2']).sample(n=int(size), random_state=1).reset_index()
            #print(sampled_regions_df)
            #print("\n")

        avg_scores_one_size_df = sampled_regions_df[['Sample1', 'Sample2', 'Synteny_score']]. \
            groupby(['Sample1', 'Sample2']).mean().reset_index().\
            rename(columns={"Synteny_score": "Avg_synteny_score"})

    if not avg_scores_one_size_df.empty:
        avg_scores_one_size_df['Compared_regions'] = size
        avg_scores_one_size_df['Ref_genome'] = genome

    #print("\nCalculate_avg_scores_selected_genome_size:")
    #print(avg_scores_one_size_df)

    return avg_scores_one_size_df


def create_pairs_num_per_sampling_size(score_per_region_df):

    regions_num_per_pair_df = score_per_region_df[['Sample1', 'Sample2', 'Synteny_score']].\
        groupby(['Sample1', 'Sample2']).count().reset_index(). \
        rename(columns={"Synteny_score": "Num_of_compared_regions"})

    #print(regions_num_per_pair_df)

    # Add a column for each subsampling size (to calculate how many pairs have results for at least this size)
    for size in config.sampling_sizes:
        if size == 'All_regions':
            regions_num_per_pair_df[size] = np.where(regions_num_per_pair_df['Num_of_compared_regions'] >= 1,
                                                     1, 0)
        else:
            regions_num_per_pair_df[size] = np.where(regions_num_per_pair_df['Num_of_compared_regions'] >= int(size),
                                                     1, 0)

    #print(regions_num_per_pair_df)

    pairs_num_per_sampling_size_df = regions_num_per_pair_df[['All_regions', '40', '60', '80', '100',
                                                              '125', '150', '175', '200', '250', '300', '350',
                                                              '400']].sum().reset_index()
    pairs_num_per_sampling_size_df.columns.values[0] = "Subsampled regions"
    pairs_num_per_sampling_size_df.columns.values[1] = "Number of pairs"
    pairs_num_per_sampling_size_df['Pairs lost percent'] = (1 - pairs_num_per_sampling_size_df['Number of pairs'] / \
                                                            pairs_num_per_sampling_size_df.at[0, 'Number of pairs']) * \
                                                            100
    pairs_num_per_sampling_size_df['Pairs lost percent'] = \
        pairs_num_per_sampling_size_df['Pairs lost percent'].apply(lambda x: round(x, 2))

    #print(pairs_num_per_sampling_size_df)
    return pairs_num_per_sampling_size_df


def create_score_per_region_sorted_contigs_table(score_per_region_df):
    contigs_dict = dict()
    contig_length_dict = dict()

    # Split the 'Region' column into
    score_per_region_df[['Contig_name', 'Position']] = score_per_region_df['Region'].str.extract(r'(\S+)_(\d+)_\d+')
    score_per_region_df['Position'] = score_per_region_df['Position'].astype(int)

    # If the contig names contain numbers, sort them numerically
    if re.search(r"^\S+_\d+$", score_per_region_df.iloc[0]['Contig_name']):

        # Create a temporary column 'contigs_sort' to sort the contig names numericlly
        score_per_region_df['Contig_number'] = score_per_region_df['Contig_name'].str.extract(r'\S+_(\d+)')\
            .astype(int)

        contigs_list_by_name = list(score_per_region_df.sort_values('Contig_number').groupby(['Contig_name'],
                                                                                             sort=False).groups)

    else:
        contigs_list_by_name = list(score_per_region_df.groupby(['Contig_name']).groups)

    #print("\ncreate_score_per_region_sorted_contigs_table:")
    #print(score_per_region_df)
    #print("\nContigs list sorted by name:")
    #print(contigs_list_by_name)

    # Create a dictionary for the contigs, sorted by their names.
    for contig in contigs_list_by_name:

        score_per_region_contig_df = score_per_region_df[score_per_region_df['Contig_name'] == contig]
        contigs_dict[contig] = score_per_region_contig_df[['Contig_name', 'Position', 'Synteny_score']]

        # Find contig length by the last position
        score_per_region_contig_df = score_per_region_contig_df.sort_values('Position')
        contig_length = score_per_region_contig_df.iloc[-1]['Position'] + config.region_length
        contig_length_dict[contig] = contig_length

    # Sort the contig lengths dict by the lengths and return a sorted list of names
    #sorted_dict = dict(sorted(contig_length_dict.items(), key=lambda item: item[1], reverse=True))
    sorted_by_length_list = sorted(contig_length_dict, key=contig_length_dict.get, reverse=True)

    return contigs_dict, contigs_list_by_name, sorted_by_length_list


