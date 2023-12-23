import re, os
import pandas as pd


def process_data(name, data):
    # create stats data
    ltr_gypsy_stats = {}
    ltr_copia_stats = {}

    # add genome column
    genome = re.findall("^[a-zA-Z]*", name)[0]
    data['Genome'] = genome
    ltr_gypsy_stats['Genome'] = genome
    ltr_copia_stats['Genome'] = genome

    # filter for intact elements
    if name.find("intact") == -1:
        data = data[data['#TE'].str.find("intact") != -1]
    ltr_gypsy_stats['Intact_Elements'] = len(data.index)
    ltr_copia_stats['Intact_Elements'] = len(data.index)

    # EDTA/TEsorter classification count
    ltr_gypsy = data[
        (data['#TE'].str.find('LTR') != -1) & (data['#TE'].str.find('Gypsy') != -1) & (data.Order == 'LTR') & (
                data.Superfamily == 'Gypsy')]
    ltr_copia = data[
        (data['#TE'].str.find('LTR') != -1) & (data['#TE'].str.find('Copia') != -1) & (data.Order == 'LTR') & (
                data.Superfamily == 'Copia')]
    ltr_gypsy_stats['EDTA/TEsorter_Classification_Count'] = len(ltr_gypsy.index)
    ltr_copia_stats['EDTA/TEsorter_Classification_Count'] = len(ltr_copia.index)

    # Remove incorrect coding domains classifications
    def has_correct_domains(clade, domains):
        split_domains = domains.split(' ')
        for domain in split_domains:
            if clade not in domain:
                return False
        return True

    ltr_gypsy = ltr_gypsy[ltr_gypsy.apply(lambda row: has_correct_domains(row['Clade'], row['Domains']), axis=1)]
    ltr_copia = ltr_copia[ltr_copia.apply(lambda row: has_correct_domains(row['Clade'], row['Domains']), axis=1)]
    ltr_gypsy_stats['Correct_Domains_Count'] = len(ltr_gypsy.index)
    ltr_copia_stats['Correct_Domains_Count'] = len(ltr_copia.index)

    # Remove clade from domains
    ltr_gypsy.Domains = ltr_gypsy.Domains.str.replace("\|[^\s]*", "", regex=True)
    ltr_copia.Domains = ltr_copia.Domains.str.replace("\|[^\s]*", "", regex=True)

    # Removing extra coding domains
    # - aRH to RH
    ltr_gypsy.Domains = ltr_gypsy.Domains.apply(lambda domains: domains.replace('aRH', 'RH'))
    ltr_copia.Domains = ltr_copia.Domains.apply(lambda domains: domains.replace('aRH', 'RH'))

    # - Remove CHD & CHDCR
    ltr_gypsy.Domains = ltr_gypsy.Domains.str.replace("CHDCR", "", regex=True).str.strip()
    ltr_gypsy.Domains = ltr_gypsy.Domains.str.replace("CHD", "", regex=True).str.strip()
    ltr_gypsy.Domains = ltr_gypsy.Domains.str.replace("\s\s", " ", regex=True).str.strip()

    ltr_copia.Domains = ltr_copia.Domains.str.replace("CHDCR", "", regex=True).str.strip()
    ltr_copia.Domains = ltr_copia.Domains.str.replace("CHD", "", regex=True).str.strip()
    ltr_copia.Domains = ltr_copia.Domains.str.replace("\s\s", " ", regex=True).str.strip()

    # Remove elements with singular gene (except GAG)
    def has_multiple_domain(domains):
        split_domains = domains.split(' ')
        return len(split_domains) > 1 or (len(split_domains) == 1 and split_domains[0] == "GAG")

    ltr_gypsy = ltr_gypsy[ltr_gypsy.apply(lambda row: has_multiple_domain(row['Domains']), axis=1)]
    ltr_copia = ltr_copia[ltr_copia.apply(lambda row: has_multiple_domain(row['Domains']), axis=1)]

    # Remove elements with incorrect domain order
    gypsy_domain_order = ["GAG", "PROT", "RT", "RH", "INT"]
    copia_domain_order = ["GAG", "PROT", "INT", "RT", "RH"]

    def has_correct_domain_order(domain_order, domains):
        split_domains = domains.split(' ')
        place_in_sequence = 0
        for domain in split_domains:
            index_of_value = domain_order.index(domain)
            if index_of_value < place_in_sequence:
                return False
            else:
                place_in_sequence = index_of_value

        return True

    ltr_gypsy = ltr_gypsy[
        ltr_gypsy.apply(lambda row: has_correct_domain_order(gypsy_domain_order, row['Domains']), axis=1)]
    ltr_copia = ltr_copia[
        ltr_copia.apply(lambda row: has_correct_domain_order(copia_domain_order, row['Domains']), axis=1)]

    # Remove elements with multiple deletion events
    def does_not_have_multiple_deletion_events(domain_order, domains):
        split_domains = domains.split(' ')
        place_in_sequence = -1
        has_deletion_event = False

        for domain in split_domains:
            index_of_value = domain_order.index(domain)

            if (index_of_value - place_in_sequence) > 1:
                if has_deletion_event:
                    return False
                has_deletion_event = True

            place_in_sequence = index_of_value

        if has_deletion_event & ((domain_order.index(split_domains[-1])) != (len(domain_order) - 1)):
            return False
        return True

    ltr_gypsy = ltr_gypsy[
        ltr_gypsy.apply(lambda row: does_not_have_multiple_deletion_events(gypsy_domain_order, row['Domains']), axis=1)]
    ltr_copia = ltr_copia[
        ltr_copia.apply(lambda row: does_not_have_multiple_deletion_events(copia_domain_order, row['Domains']), axis=1)]
    ltr_gypsy_stats['Processed_Domains_Count'] = len(ltr_gypsy.index)
    ltr_copia_stats['Processed_Domains_Count'] = len(ltr_copia.index)

    return [ltr_gypsy, ltr_copia, pd.DataFrame.from_dict(ltr_gypsy_stats, orient='index').T,
            pd.DataFrame.from_dict(ltr_copia_stats, orient='index').T]


directory = 'all-data'
stat_out = [pd.DataFrame(), pd.DataFrame()]
data_out = [pd.DataFrame(), pd.DataFrame()]

for filename in os.listdir(directory):
    print(filename)
    print("\n")
    processed_data = process_data(filename, pd.read_csv(directory + "/" + filename, sep='\t'))
    data_out = [pd.concat([data_out[0], processed_data[0]], ignore_index=True),
                pd.concat([data_out[1], processed_data[1]], ignore_index=True)]
    stat_out = [pd.concat([stat_out[0], processed_data[2]], ignore_index=True),
                pd.concat([stat_out[1], processed_data[3]], ignore_index=True)]

if not os.path.exists('out/data'):
    os.makedirs('out/data')

if not os.path.exists('out/stats'):
    os.makedirs('out/stats')

data_out[0].to_csv('out/data/ltr_gypsy_processed.tsv', sep="\t")
data_out[1].to_csv('out/data/ltr_copia_processed.tsv', sep="\t")
stat_out[0].to_csv('out/stats/ltr_gypsy_stats.tsv', sep="\t")
stat_out[1].to_csv('out/stats/ltr_copia_stats.tsv', sep="\t")
