import glob
import pandas as pd
file_list = glob.glob('/home/b/Desktop/LC-MS/individual_analysis_1/*')


table=list()

for file in file_list:
    a=list()
    b=list()
    c=list()
    d=list()

    with open(file) as in_file:
        for line in in_file:
            line = line.strip().split('\t')
            a.append(line[0])
            b.append(line[0])
            c.append(line[0])
            d.append(line[0])
    table.append(a)
    table.append(b)
    table.append(c)
    table.append(d)

df = pd.DataFrame(table)
df.to_csv('test')