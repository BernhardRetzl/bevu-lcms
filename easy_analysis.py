import glob


class Hit:
    def __init__(self, sequence, max_int, ret_time, mass, n_found):
        self.sequence = sequence
        self.max_int = max_int
        self.ret_time = ret_time
        self.mass = mass
        self.n_found = n_found


def read_file():
    my_dict = dict()
    with open('/home/b/Desktop/LC-MS/Analysis/fragments_1/sample1') as in_file:
        for line in in_file:
            line = line.strip().split()
            if line[0] == 'Sequence':
                sequence = line[1]
                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                if sequence in my_dict:
                    my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                                 ret_time=ret_time, mass=mass, n_found=n_found))
                else:
                    my_dict[sequence] = [Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found)]

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found))

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found))

    return my_dict


def read_control_file():
    my_dict = dict()
    with open('/home/b/Desktop/LC-MS/Analysis/fragments_1/sample1_K') as in_file:
        for line in in_file:
            line = line.strip().split()
            if line[0] == 'Sequence':
                sequence = line[1]
                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                if sequence in my_dict:
                    my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                                 ret_time=ret_time, mass=mass, n_found=n_found))
                else:
                    my_dict[sequence] = [Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found)]

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found))

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].append(Hit(sequence=sequence, max_int=float(max_int),
                                             ret_time=ret_time, mass=mass, n_found=n_found))

    return my_dict


this_dict = read_file()
this_control_dict = read_control_file()


for item in this_dict:
    if this_dict[item][0].max_int > this_control_dict[item][0].max_int*5:
        print(this_dict[item][0].max_int, this_control_dict[item][0].max_int)