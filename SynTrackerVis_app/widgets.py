
def create_pairs_lost_text(df, sampling_size):
    pairs_num = df.loc[df['Subsampled regions'] == sampling_size, 'Number of pairs'].values[0]
    pairs_lost = df.loc[df['Subsampled regions'] == sampling_size, 'Pairs lost percent'].values[0]

    return {
        'object': "Number of compared pairs: " + str(pairs_num) + ".&nbsp;&nbsp;&nbsp;Lost pairs: " +
                  str(pairs_lost) + " %",
        'styles': {'font-size': "16px"}
    }


def get_contig_title(contig_name):
    contig_title = "Contig name: " + contig_name
    return contig_title
