import re, os
import pandas as pd

clade = "Reina"

def get_length(row):
    [start_loc, end_loc] = row['Length'].split("..")
    return int(end_loc) - int(start_loc)

def get_athila_lengths(data):
    data = data[data["Clade"] == clade]
    data['Length'] = data['#TE'].str.extract("(\d+\.\.\d+)", expand=True)
    data["Length"] = data.apply(get_length, axis=1)
    return data 

    #TODO get length of all rows for both files then seperate by clade into files.


directory = 'out/data'

if not os.path.exists('out/clade_length'):
    os.makedirs('out/clade_length')



processed_data = get_athila_lengths(pd.read_csv(directory + "/ltr_gypsy_processed.tsv", sep='\t'))
processed_data.to_csv('out/clade_length/' + clade + '_length.tsv', sep="\t")

