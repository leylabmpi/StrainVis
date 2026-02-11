import pandas as pd
import numpy as np
import StrainVis_app.config as config


def complete_metadata(score_per_region_df, metadata_df):

    metadata_dict = dict()
    groups_per_feature_dict = dict()
    error_msg = ""
    sample_ids_column_name = ""

    # Extract the metadata feature names and save them in the list
    metadata_features_list = list(metadata_df.columns)
    print("\nMetadata features:")
    print(metadata_features_list)

    # There are less than two columns - probably a delimiter problem
    if len(metadata_features_list) < 2:
        error_msg = "The provided metadata file contains only one column - probably the delimiter is not correct.  \n" \
                    "Please provide a new metadata file, saved in tab-delimited format."
        return metadata_dict, metadata_features_list, sample_ids_column_name, error_msg

    # Extract the name of the first column in the metadata (the sample_IDs column)
    sample_ids_column_name = metadata_features_list.pop(0)
    #print("\nFirst column name:")
    #print(sample_ids_column_name)

    # Extract a list of unique sample IDs from the metadata
    metadata_orig_sample_list = metadata_df[sample_ids_column_name].to_list()
    print("\nMetadata contains " + str(len(metadata_orig_sample_list)) + " unique samples")

    # Extract a list of unique sample IDs from SynTracker results
    unique_samples_list = pd.concat([score_per_region_df['Sample1'], score_per_region_df['Sample2']],
                                    ignore_index=True). \
        drop_duplicates().to_list()
    print("Data contains " + str(len(unique_samples_list)) + " unique samples")

    # Go over the samples in the data. For those which are missing from the metadata, fill all the fields with 'NaN'
    features_num = len(metadata_features_list)
    for sample in unique_samples_list:
        if sample not in metadata_orig_sample_list:
            #print("Sample " + sample + " is missing from metadata")
            new_row = [sample]
            for i in range(features_num):
                new_row.append(np.nan)
                #new_row.append("NaN")
            metadata_df.loc[len(metadata_df)] = new_row

    #print("\nMetadata after filling missing samples:")
    #print(metadata_df)

    # Go over the features
    for feature in metadata_features_list:
        # 1. Create a dictionary to map the samples to feature values from metadata_df
        metadata_dict[feature] = metadata_df.set_index(sample_ids_column_name)[feature].to_dict()

        # 2. Create another dict, with features -> feature groups
        groups_per_feature_dict[feature] = metadata_df[feature].unique().tolist()

    return metadata_dict, groups_per_feature_dict, metadata_features_list, error_msg


def return_genomes_subset_table(score_per_region_df, genomes_list):
    genomes_subset_df = score_per_region_df[score_per_region_df['Ref_genome'].isin(genomes_list)]
    genomes_subset_df = genomes_subset_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']]

    #print("\nreturn_genomes_subset_table:")
    #print(genomes_subset_df)
    return genomes_subset_df


def count_species_num(row, df):
    condition = df[row['Subsampled_regions']] > 0
    species_num = condition.sum()
    return species_num


def filter_genomes_ani(ani_scores_all_genomes_df):

    pairs_num_df = ani_scores_all_genomes_df[['Ref_genome', 'ANI']].groupby('Ref_genome').count().\
        sort_values('ANI', ascending=False).reset_index()
    pairs_num_df.columns.values[1] = "Number_of_pairs"

    # Leave only genomes, which have at lease 10 pairwise comparisons
    pairs_num_filtered_df = pairs_num_df[pairs_num_df["Number_of_pairs"] >= 10]
    #print("\nfilter_genomes_ani:")
    #print(pairs_num_filtered_df)

    filtered_genomes_list = list(pairs_num_filtered_df['Ref_genome'])

    ani_scores_filtered_genomes_df = ani_scores_all_genomes_df[
        ani_scores_all_genomes_df['Ref_genome'].isin(filtered_genomes_list)]

    return ani_scores_filtered_genomes_df, filtered_genomes_list


def create_sorted_by_pairs_genomes_list_syntracker(score_per_region_all_genomes_df):
    regions_num_per_pair_df = score_per_region_all_genomes_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']]. \
        groupby(['Ref_genome', 'Sample1', 'Sample2']).count().reset_index(). \
        rename(columns={"Synteny_score": "Num_of_compared_regions"})

    regions_num_per_pair_df['40'] = np.where(regions_num_per_pair_df['Num_of_compared_regions'] >= 40, 1, 0)

    pairs_num_at_40_regions_df = regions_num_per_pair_df[['Ref_genome', '40']].groupby('Ref_genome').sum().\
        sort_values('40', ascending=False).reset_index()
    pairs_num_at_40_regions_df.columns.values[1] = "Number_of_pairs"

    #print("\ncreate_sorted_by_pairs_genomes_list_syntracker:")
    #print(pairs_num_at_40_regions_df)

    genomes_list_by_pairs_num = list(pairs_num_at_40_regions_df['Ref_genome'])
    #print("\nGenomes list sorted by pairs number:")
    #print(genomes_list_by_pairs_num)

    return genomes_list_by_pairs_num


def create_sorted_by_pairs_genomes_list_ani(ani_scores_all_genomes_df):

    pairs_num_df = ani_scores_all_genomes_df[['Ref_genome', 'ANI']].groupby('Ref_genome').count().\
        sort_values('ANI', ascending=False).reset_index()
    pairs_num_df.columns.values[1] = "Number_of_pairs"

    print("\ncreate_sorted_by_pairs_genomes_list_ani:")
    print(pairs_num_df)

    genomes_list_by_pairs_num = list(pairs_num_df['Ref_genome'])
    #print("\nGenomes list sorted by pairs number:")
    #print(genomes_list_by_pairs_num)

    return genomes_list_by_pairs_num


def create_pairs_num_per_sampling_size(score_per_region_selected_genomes_df):

    regions_num_per_pair_df = score_per_region_selected_genomes_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']].\
        groupby(['Ref_genome', 'Sample1', 'Sample2']).count().reset_index(). \
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

    pairs_num_per_sampling_size_df = regions_num_per_pair_df[['Ref_genome', 'All', '40', '60', '80', '100',
                                                              '125', '150', '175', '200', '250', '300', '350',
                                                              '400']].groupby('Ref_genome').sum().reset_index()
    #print(pairs_num_per_sampling_size_df)

    for size in config.sampling_sizes:
        pairs_num_per_sampling_size_df[size] = np.where(pairs_num_per_sampling_size_df[size] >= 10,
                                                        pairs_num_per_sampling_size_df[size], 0)
    #print(pairs_num_per_sampling_size_df)

    summary_df = pairs_num_per_sampling_size_df[['All', '40', '60', '80', '100', '125', '150', '175', '200',
                                                 '250', '300', '350', '400']].sum().reset_index()

    summary_df.columns.values[0] = "Subsampled_regions"
    summary_df.columns.values[1] = "Number_of_pairs"

    summary_df['Pairs_lost_percent'] = (1 - summary_df['Number_of_pairs'] / summary_df.at[0, 'Number_of_pairs']) * 100
    summary_df['Pairs_lost_percent'] = summary_df['Pairs_lost_percent'].apply(lambda x: round(x, 2))

    summary_df['Number_of_species'] = \
        summary_df.apply(lambda row: count_species_num(row, pairs_num_per_sampling_size_df), axis=1)

    #print("\ncreate_pairs_num_per_sampling_size:")
    #print(summary_df)

    return summary_df


def calculate_APSS_all_genomes_sampling_size(score_per_region_df, size):

    # Taking all available regions - no subsampling
    if size == 'All':
        avg_scores_one_size_df = score_per_region_df.groupby(['Ref_genome', 'Sample1', 'Sample2'])['Synteny_score'].\
            mean().reset_index().rename(columns={"Synteny_score": "APSS"})

    # All the other optional sizes -
    # need to include only the sample pairs that have a result for at least the requested number of regions
    else:
        filtered_df = score_per_region_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']].\
            groupby(['Ref_genome', 'Sample1', 'Sample2']).filter(lambda x: x['Synteny_score'].count() >= int(size))
        #print(filtered_df)

        sampled_regions_df = filtered_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']].\
            groupby(['Ref_genome', 'Sample1', 'Sample2']).sample(n=int(size), random_state=1).reset_index()
        #print(sampled_regions_df)

        avg_scores_one_size_df = sampled_regions_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']]. \
            groupby(['Ref_genome', 'Sample1', 'Sample2']).mean().reset_index().\
            rename(columns={"Synteny_score": "APSS"})

    if not avg_scores_one_size_df.empty:
        avg_scores_one_size_df['Compared_regions'] = size

    #print("\ncalculate_APSS_all_genomes_sampling_size:")
    #print(avg_scores_one_size_df)

    # Filter out species with less than 10 pairs
    samples_per_genome_df = avg_scores_one_size_df[['Ref_genome', 'APSS']].groupby('Ref_genome').count().reset_index().\
        rename(columns={"APSS": "count"})
    #print(samples_per_genome_df)
    merged_df = avg_scores_one_size_df.merge(samples_per_genome_df[['Ref_genome', 'count']], on='Ref_genome',
                                             how='left')
    avg_scores_one_size_filtered_df = merged_df[merged_df['count'] >= 10].drop(columns='count')
    #print("\ncalculate_APSS_all_genomes_sampling_size after genomes filtering:")
    #print(avg_scores_one_size_filtered_df)

    return avg_scores_one_size_filtered_df


def return_genomes_subset_APSS_selected_size_table(all_genomes_selected_size_APSS_df, genomes_list):
    genomes_subset_selected_size_APSS_df = \
        all_genomes_selected_size_APSS_df[all_genomes_selected_size_APSS_df['Ref_genome'].isin(genomes_list)].\
        reset_index()

    #print("\nreturn_genomes_subset_APSS_selected_size_table:")
    #print(genomes_subset_selected_size_APSS_df)
    return genomes_subset_selected_size_APSS_df


def return_genomes_subset_ani_table(all_genomes_ani_df, genomes_list):
    genomes_subset_ani_df = all_genomes_ani_df[all_genomes_ani_df['Ref_genome'].isin(genomes_list)].reset_index()

    #print("\return_genomes_subset_ani_table:")
    #print(genomes_subset_ani_df)
    return genomes_subset_ani_df


