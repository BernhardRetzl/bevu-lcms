import glob


file_list_junk = glob.glob('/home/b/Desktop/LC-MS/Analysis/*')
file_list_junk.sort()
file_list = [i for i in file_list_junk if '_K' not in i]
file_list_junk = glob.glob('/home/b/Desktop/LC-MS/Analysis/*')
file_list_junk.sort()
file_list += [i for i in file_list_junk if '_K' in i]
my_dict = dict()


class Hit:
    def __init__(self):
        self.single_items = list()
        self.double_items = list()
        self.triple_items = list()


for file in file_list:
    with open(file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')
            if line[0] == 'Sequence':
                sequence = line[1]
                if sequence not in my_dict:
                    my_dict[sequence] = Hit()
                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].single_items.append((n_found, max_int, ret_time, mass))

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].double_items.append((n_found, max_int, ret_time, mass))

                n_found = in_file.readline().strip().split('\t')[1]
                junk, max_int, ret_time, mass = in_file.readline().strip().split('\t')
                my_dict[sequence].triple_items.append((n_found, max_int, ret_time, mass))

for item in my_dict:
    with open('/home/b/Desktop/LC-MS/Individual_analysis_1/'+item, 'wt') as out_file:
        out_file.write(item+'\n')
        for i in my_dict[item].single_items:
            out_file.write('\t'.join(i)+'\n')
        for i in my_dict[item].double_items:
            out_file.write('\t'.join(i)+'\n')
        for i in my_dict[item].triple_items:
            out_file.write('\t'.join(i)+'\n')