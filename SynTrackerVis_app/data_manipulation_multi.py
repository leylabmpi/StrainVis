import pandas as pd
import numpy as np
import SynTrackerVis_app.config as config


def create_avg_scores_table(score_per_region_df):
    dfs_list = []
    #regions_num_per_pair_df = score_per_region_df.groupby(['Ref_genome', 'Sample1', 'Sample2'])['Synteny_score'].\
        #count().reset_index().rename(columns={"Synteny_score": "Available_regions"})
    #print(regions_num_per_pair_df)

    for size in config.sampling_sizes:
        # Taking all available regions - no subsampling
        if size == 'All_regions':
            avg_scores_one_size_df = score_per_region_df.groupby(['Ref_genome', 'Sample1', 'Sample2'])['Synteny_score'].\
                mean().reset_index().rename(columns={"Synteny_score": "Avg_synteny_score"})
            avg_scores_one_size_df['Compared_regions'] = size

            if not avg_scores_one_size_df.empty:
                dfs_list.append(avg_scores_one_size_df)

            #print(avg_scores_one_size_df)
            #print("\n")

        else:
            filtered_df = score_per_region_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']].\
                groupby(['Ref_genome', 'Sample1', 'Sample2']).\
                filter(lambda x: x['Synteny_score'].count() >= int(size))
            #print(filtered_df)

            if not filtered_df.empty:

                sampled_regions_df = filtered_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']].\
                    groupby(['Ref_genome', 'Sample1', 'Sample2']).sample(n=int(size), random_state=1).reset_index()
                #print(sampled_regions_df)
                #print("\n")

                avg_scores_one_size_df = sampled_regions_df[['Ref_genome', 'Sample1', 'Sample2', 'Synteny_score']]. \
                    groupby(['Ref_genome', 'Sample1', 'Sample2']).mean().reset_index().\
                    rename(columns={"Synteny_score": "Avg_synteny_score"})
                avg_scores_one_size_df['Compared_regions'] = size
                #print(avg_scores_one_size_df)
                #print("\n")

                if not avg_scores_one_size_df.empty:
                    dfs_list.append(avg_scores_one_size_df)

    # Concatenate the 'all regions' DF to the main dataframe including all the sizes
    avg_scores_all_genomes_all_sizes_df = pd.concat(dfs_list)
    print(avg_scores_all_genomes_all_sizes_df)
    print("\n")

    return avg_scores_all_genomes_all_sizes_df


def complete_metadata(score_per_region_df, metadata_df):

    metadata_dict = dict()

    # Extract the metadata feature names and save them in the list
    metadata_features_list = list(metadata_df.columns)
    print("\nMetadata features:")
    print(metadata_features_list)

    # Extract the name of the first column in the metadata (the sample_IDs column)
    sample_ids_column_name = metadata_features_list.pop(0)
    print("\nFirst column name:")
    print(sample_ids_column_name)

    # Extract a list of unique sample IDs from the metadata
    metadata_orig_sample_list = metadata_df[sample_ids_column_name].to_list()
    print("\nMetadata contains " + str(len(metadata_orig_sample_list)) + " unique samples")

    # Extract a list of unique sample IDs from SynTracker results
    unique_samples_list = pd.concat([score_per_region_df['Sample1'], score_per_region_df['Sample2']],
                                    ignore_index=True). \
        drop_duplicates().to_list()
    print("\nData contains " + str(len(unique_samples_list)) + " unique samples")

    # Go over the samples in the data. For those which are missing from the metadata, fill all the fields with 'NaN'
    features_num = len(metadata_features_list)
    for sample in unique_samples_list:
        if sample not in metadata_orig_sample_list:
            print("Sample " + sample + " is missing from metadata")
            new_row = [sample]
            for i in range(features_num):
                new_row.append(np.nan)
                #new_row.append("NaN")
            metadata_df.loc[len(metadata_df)] = new_row

    print("\nMetadata after filling missing samples:")
    print(metadata_df)

    # Create a dictionary to map the samples to feature values from metadata_df
    for feature in metadata_features_list:
        metadata_dict[feature] = metadata_df.set_index(sample_ids_column_name)[feature].to_dict()

    return metadata_dict, metadata_features_list, sample_ids_column_name
