
def create_pairs_lost_text(df, sampling_size):
    pairs_num = df.loc[df['Subsampled_regions'] == sampling_size, 'Number_of_pairs'].values[0]
    pairs_lost = df.loc[df['Subsampled_regions'] == sampling_size, 'Pairs_lost_percent'].values[0]

    return {
        'object': "Number of compared pairs: " + str(pairs_num) + ".&nbsp;&nbsp;&nbsp;Lost pairs: " +
                  str(pairs_lost) + " %",
        'styles': {'font-size': "16px", 'text-align': "center"}
    }


def create_species_num_text(df, sampling_size):
    species_num = df.loc[df['Subsampled_regions'] == sampling_size, 'Number_of_species'].values[0]
    total_species_num = df.at[0, 'Number_of_species']

    return {
        'object': "Number of species: " + str(species_num) + " out of " + str(total_species_num),
        'styles': {'font-size': "16px", 'text-align': "center"}
    }


def create_samples_num_text(df, sampling_size):
    samples_num = df.loc[df['Subsampled_regions'] == sampling_size, 'Number_of_samples'].values[0]
    total_samples_num = df.at[0, 'Number_of_samples']

    return {
        'object': "Number of samples: " + str(samples_num) + " out of " + str(total_samples_num),
        'styles': {'font-size': "16px", 'text-align': "center"}
    }


def get_contig_title(contig_name):
    contig_title = "Contig name: " + contig_name
    return contig_title
