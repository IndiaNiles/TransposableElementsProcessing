import re, os
import pandas as pd


def get_length(row):
    [start_loc, end_loc] = row['Length'].split("..")
    return int(end_loc) - int(start_loc)

def get_athila_lengths(data):
    data = data[data["Clade"] == "Athila"]
    data['Length'] = data['#TE'].str.extract("(\d+\.\.\d+)", expand=True)
    data["Length"] = data.apply(get_length, axis=1)
    return data 


directory = 'out/data'

if not os.path.exists('out/athila_length'):
    os.makedirs('out/athila_length')



processed_data = get_athila_lengths(pd.read_csv(directory + "/ltr_gypsy_processed.tsv", sep='\t'))
processed_data.to_csv('out/athila_length/athila_length.tsv', sep="\t")

