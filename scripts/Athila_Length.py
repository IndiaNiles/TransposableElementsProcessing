import re, os
import pandas as pd


def get_length(row):
    [start_loc, end_loc] = row["Length"].split("..")
    return int(end_loc) - int(start_loc)


def add_lengths(data):
    data["Length"] = data["#TE"].str.extract("(\d+\.\.\d+)", expand=True)
    data["Length"] = data.apply(get_length, axis=1)
    return data


def clades_to_file(clades, processed_data):
    for clade in clades:
        out_data = processed_data[processed_data["Clade"] == clade]
        out_data.to_csv(out_directory + "/" + clade + "_length.tsv", sep="\t")


in_directory = "out/data"
out_directory = "out/clade_lengths"

if not os.path.exists(out_directory):
    os.makedirs(out_directory)

for filename in os.listdir(in_directory):
    data = add_lengths(pd.read_csv(in_directory + "/" + filename, sep="\t"))
    clades_to_file(data["Clade"].unique(), data)
