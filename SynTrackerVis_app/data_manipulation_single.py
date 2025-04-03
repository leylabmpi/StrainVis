import pandas as pd
import numpy as np
import re
import time
import SynTrackerVis_app.config as config


def return_selected_genome_table(score_per_region_df, selected_genome):
    selected_genome_df = score_per_region_df[score_per_region_df['Ref_genome'] == selected_genome]
    selected_genome_df = selected_genome_df[['Sample1', 'Sample2', 'Region', 'Synteny_score']]

    #print(selected_genome_df)
    return selected_genome_df


def return_selected_genome_avg_table(avg_big_df, selected_genome):
    selected_genome_avg_df = avg_big_df[avg_big_df['Ref_genome'] == selected_genome]
    selected_genome_avg_df = selected_genome_avg_df[['Sample1', 'Sample2', 'APSS', 'Compared_regions']]

    #print(selected_genome_avg_df)
    return selected_genome_avg_df


def calculate_avg_scores_selected_genome_size(score_per_region_selected_genome_df, genome, size):

    # Taking all available regions - no subsampling
    if size == 'All':
        avg_scores_one_size_df = score_per_region_selected_genome_df.groupby(['Sample1', 'Sample2'])['Synteny_score'].\
            mean().reset_index().rename(columns={"Synteny_score": "APSS"})

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
            rename(columns={"Synteny_score": "APSS"})

    if not avg_scores_one_size_df.empty:
        avg_scores_one_size_df['Compared_regions'] = size
        avg_scores_one_size_df['Ref_genome'] = genome

    #print("\nCalculate_avg_scores_selected_genome_size:")
    #print(avg_scores_one_size_df)

    return avg_scores_one_size_df


def count_samples_num(row, df):
    condition_df = df[df[row['Subsampled_regions']] > 0]

    unique_samples = pd.concat([condition_df['Sample1'], condition_df['Sample2']]).unique()
    samples_num = len(unique_samples)

    return samples_num


def create_pairs_num_per_sampling_size(score_per_region_df):
    print("\ncreate_pairs_num_per_sampling_size:")

    regions_num_per_pair_df = score_per_region_df[['Sample1', 'Sample2', 'Synteny_score']].\
        groupby(['Sample1', 'Sample2']).count().reset_index(). \
        rename(columns={"Synteny_score": "Num_of_compared_regions"})

    #print(regions_num_per_pair_df)

    # Add a column for each subsampling size (to calculate how many pairs have results for at least this size)
    for size in config.sampling_sizes:
        if size == 'All':
            regions_num_per_pair_df[size] = np.where(regions_num_per_pair_df['Num_of_compared_regions'] >= 1,
                                                     1, 0)
        else:
            regions_num_per_pair_df[size] = np.where(regions_num_per_pair_df['Num_of_compared_regions'] >= int(size),
                                                     1, 0)

    #print(regions_num_per_pair_df)

    pairs_num_per_sampling_size_df = regions_num_per_pair_df[['All', '40', '60', '80', '100',
                                                              '125', '150', '175', '200', '250', '300', '350',
                                                              '400']].sum().reset_index()
    pairs_num_per_sampling_size_df.columns.values[0] = "Subsampled_regions"
    pairs_num_per_sampling_size_df.columns.values[1] = "Number_of_pairs"
    pairs_num_per_sampling_size_df['Pairs_lost_percent'] = (1 - pairs_num_per_sampling_size_df['Number_of_pairs'] / \
                                                            pairs_num_per_sampling_size_df.at[0, 'Number_of_pairs']) * \
                                                            100
    pairs_num_per_sampling_size_df['Pairs_lost_percent'] = \
        pairs_num_per_sampling_size_df['Pairs_lost_percent'].apply(lambda x: round(x, 2))

    # Calculate the number of samples in each sampling size
    pairs_num_per_sampling_size_df['Number_of_samples'] = \
        pairs_num_per_sampling_size_df.apply(lambda row: count_samples_num(row, regions_num_per_pair_df), axis=1)

    print(pairs_num_per_sampling_size_df)

    return pairs_num_per_sampling_size_df


def return_sorted_contigs_lists(score_per_region_df):
    before = time.time()
    # Split the 'Region' column into Contig_name and Position
    pattern = re.compile(r'(\S+)_(\d+)_\d+')
    score_per_region_df[['Contig_name', 'Position']] = score_per_region_df['Region'].str.extract(pattern)
    score_per_region_df['Position'] = score_per_region_df['Position'].astype(int)

    after = time.time()
    duration = after - before
    print("Extract position from region took " + str(duration) + " seconds.\n")

    # Get a list of contigs, sorted by name
    # If the contig names contain numbers, sort them numerically
    #if re.search(r"^\S+_\d+$", score_per_region_df.iloc[0]['Contig_name']):

        # Create a temporary column 'contigs_sort' to sort the contig names numericlly
    #    score_per_region_df['Contig_number'] = score_per_region_df['Contig_name'].str.extract(r'\S+_(\d+)')\
    #        .astype(int)

    #    contigs_list_by_name = list(score_per_region_df.sort_values('Contig_number').groupby(['Contig_name'],
    #                                                                                         sort=False).groups)

    #else:
    #    contigs_list_by_name = list(score_per_region_df.groupby(['Contig_name']).groups)

    before = time.time()
    contigs_list_by_name = list(score_per_region_df.groupby(['Contig_name']).groups)
    after = time.time()
    duration = after - before
    print("Sort by name took " + str(duration) + " seconds.\n")

    # Get a list of contigs, sorted by length
    before = time.time()
    contigs_list_by_length = list(score_per_region_df.sort_values('Position', ascending=False).groupby(['Contig_name'],
                                                                                                       sort=False).
                                  groups)
    after = time.time()
    duration = after - before
    print("Sort by length took " + str(duration) + " seconds.\n")
  
    return contigs_list_by_name, contigs_list_by_length


