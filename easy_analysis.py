import glob

class Hit:
    def __init__(self, sequence, max_int, ret_time, mass):
        self.sequence = sequence
        self.max_int = max_int
        self.ret_time = ret_time
        self.mass = mass


def read_file():
    my_dict = dict()
    with open('/home/b/Desktop/LC-MS/Analysis/fragments_1/') as in_file:
        for line in in_file:
            line = line.strip().split()
            if line[0] == 'Sequence':
                sequence = line[1]
                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
