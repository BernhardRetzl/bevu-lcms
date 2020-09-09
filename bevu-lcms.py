from Bio.SeqUtils import molecular_weight
import glob
from pyopenms import *

water = 2*molecular_weight('G', monoisotopic=True, seq_type='protein') - \
        molecular_weight('GG', monoisotopic=True, seq_type='protein')


def get_small_fragments(file):
    help_list = list()
    with open(file) as in_file:
        for line in in_file:
            line = line.strip()
            help_list.append(line)
    return help_list


class ItemToFind:
    def __init__(self, sequence, single, double, triple):
        self.sequence = sequence
        self.single=single
        self.double=double
        self.triple=triple
        self.single_found = list()
        self.double_found = list()
        self.triple_found = list()

    def __lt__(self, other):
        return len(self.single_found+self.double_found+self.triple_found)<\
        len(other.single_found+other.double_found+other.triple_found)


class Hit:
    def __init__(self, mass, intensity, retention_tiem):
        self.mass = mass
        self.intensity = intensity
        self.retention_time = retention_tiem


def calculate_mass(my_list):
    help_list = list()
    for item in my_list:
        single = molecular_weight(item, monoisotopic=True, seq_type='protein')+water-1*1.007825035
        double = (molecular_weight(item, monoisotopic=True, seq_type='protein')+water-0*1.007825035)/2
        triple = (molecular_weight(item, monoisotopic=True, seq_type='protein')+water+1*1.007825035)/3
        help_list.append(ItemToFind(sequence=item,
                         single=single,
                         double=double,
                         triple=triple))
    return help_list


def analyse_spectrum(spectrum, mass_list, output_name):
    exp = MSExperiment()
    MzMLFile().load(spectrum, exp)

    for s in exp.getSpectra():
        peaks = s.get_peaks()
        mass_int = list(zip(peaks[0], peaks[1]))
        retention_time = s.getRT()
        for mass in mass_int:
            for item in mass_list:
                if abs(mass[0]-item.single) <= 0.1 and mass[1] > 10000:
                    item.single_found.append(Hit(mass=mass[0], intensity=mass[1], retention_tiem=retention_time))
                if abs(mass[0]-item.double) <= 0.1 and mass[1] > 10000:
                    item.double_found.append(Hit(mass=mass[0], intensity=mass[1], retention_tiem=retention_time))
                if abs(mass[0]-item.triple) <= 0.1 and mass[1] > 10000:
                    item.triple_found.append(Hit(mass=mass[0], intensity=mass[1], retention_tiem=retention_time))

    mass_list = sorted(mass_list)

    with open('/home/b/Desktop/LC-MS/Analysis/'+output_name, 'wt') as out_file:
        for item in mass_list:
            single_maximum = max(item.single_found, key=lambda x: x.intensity, default=Hit(0,0,0))
            double_maximum = max(item.double_found, key=lambda x: x.intensity, default=Hit(0,0,0))
            triple_maximum = max(item.triple_found, key=lambda x: x.intensity, default=Hit(0,0,0))
            out_file.write('Sequence\t'+item.sequence+'\n'+ 'Number of single found\t'+
                           str(len(item.single_found))+'\n'+ 'max intensity, retime, mass single'+'\t'+
                           str(single_maximum.intensity)+'\t'+str(single_maximum.retention_time)+'\t'+str(round(single_maximum.mass,2))+'\n'+

            'Number of double found\t' +
            str(len(item.double_found)) + '\n' + 'max intensity, retime, mass double' + '\t' +
            str(double_maximum.intensity) + '\t' + str(double_maximum.retention_time) + '\t' + str(
                round(double_maximum.mass, 2)) + '\n'+

            'Number of triple found\t' +
            str(len(item.triple_found)) + '\n' + 'max intensity, retime, mass triple' + '\t' +
            str(triple_maximum.intensity) + '\t' + str(triple_maximum.retention_time) + '\t' + str(
                round(triple_maximum.mass, 2)) + '\n')


def start(fragment_file):
    small_fragments = get_small_fragments(fragment_file)
    file_list = glob.glob('/home/b/Desktop/LC-MS/*.mzML')
    file_list.sort()
    for file in file_list:
        my_mass_list = calculate_mass(my_list=small_fragments)
        print(file)
        sample_number = int(file.split('sample')[1][0])
        control = True if '_K_' in file else False
        day_2 = True if '_2d_' in file else False
        if day_2:
            sample_number += 6
        new_name = fragment_file+'/'+'sample'+str(sample_number)
        new_name += '_K' if control else ''
        print(new_name)
        analyse_spectrum(spectrum=file, mass_list=my_mass_list, output_name=new_name)


start(fragment_file='fragments_3')
