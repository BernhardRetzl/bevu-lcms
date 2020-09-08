from Bio.SeqUtils import molecular_weight


water = 2*molecular_weight('G', monoisotopic=True, seq_type='protein') - \
        molecular_weight('GG', monoisotopic=True, seq_type='protein')

small_fragments_1 = ['TGSPTCCSNS',
                     'TGSPTCCS',
                     'TGSPTCC',
                     'GSPTCCSNS',
                     'GSPTCCS',
                     'GSPTCC',
                     'SPTCCSNS',
                     'SPTCCS',
                     'SPTCC'
                     'CCSNS',
                     'CCS',
                     'CC']


def calculate_mass(my_list):
    my_dict = dict()
    for item in my_list:
        my_dict[item] = molecular_weight(item, monoisotopic=True, seq_type='protein')+water-2*1.007825035
    return my_dict


masses_of_interest = calculate_mass(my_list=small_fragments_1)



