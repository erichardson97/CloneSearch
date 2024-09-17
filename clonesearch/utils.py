#!/usr/bin/env python3
from Levenshtein import ratio
from collections import defaultdict
from typing import DefaultDict, Callable
import pandas as pd
import os
import requests
def ascs_ogrdb(ogrdb_file_path: (str | type(None)) = None):
    if ogrdb_file_path != None:
        ogrdb_file = pd.read_json(ogrdb_file_path)
    else:
        download = download_ogrdb()
        ogrdb_file = pd.read_json(download)
    alleles = ogrdb_file.iloc[0]['GermlineSet']['allele_descriptions']
    asc_mapping = {}
    for k in alleles:
        asc_mapping[k['label']] = k['label'] if k['allele_similarity_cluster_designation'] == None else k['allele_similarity_cluster_designation']
    return asc_mapping

def ascs_github(tsv_file_path: (str | type(None)) = f'asc_cluster_cut.tsv', key = 'ASC_Cluster_0.25'):
    if (tsv_file_path != None) & (os.path.exists(tsv_file_path)):
        asc_file = pd.read_csv(tsv_file_path, sep = '\t')
    else:
        output = requests.get(f'https://raw.githubusercontent.com/yaarilab/asc_archive/main/asc_cluster_cut.tsv').content.decode()
        with open(f'asc_cluster_cut.tsv', 'w') as k:
            k.writelines(output)
        asc_file = pd.read_csv(f'asc_cluster_cut.tsv', sep = '\t')
    asc_mapping = dict(asc_file[['allele', 'ASC_Cluster_0.05']].values)
    return asc_mapping

def download_ogrdb():
    return None


ascs_dict = ascs_github('asc_cluster_cut.tsv')
ascs_family = ascs_github('asc_cluster_cut.tsv', key='ASC_Cluster_0.25')

def unpack_genes(v_field: str):
    return ','.join(set([p.split('*')[0] for p in v_field.split(',')]))

def make_hash(dataframe: pd.DataFrame, v_field: str = 'v_call', cdr3_field: str = 'cdr3_aa',
              allele: bool = True, sequence_id: str = 'sequence_id', use_v: bool = True,
              use_j: bool = False, j_field: str = 'j_call', use_asc: bool = False) -> DefaultDict:
    '''

    :param dataframe:
    :param v_field:
    :param cdr3_field:
    :param allele:
    :param sequence_id:
    :param use_j:
    :param j_field:
    :return: Dictionary with keys equal to V, (J), CDR3 Length,
    '''
    group = ['cdr3_length']
    dataframe = dataframe[(dataframe[cdr3_field].notna())&(dataframe[v_field].notna())]
    dataframe['cdr3_length'] = dataframe[cdr3_field].fillna('').map(len)
    if use_v:
        group.append('v')
        if allele:
            dataframe['v'] = dataframe[v_field]
            if use_asc == 'cluster':
                dataframe['v'] = dataframe['v'].map(lambda x:','.join([str(ascs_dict[k]) if k in ascs_dict else k for k in x.split(',')]))
            elif use_asc == 'family':
                dataframe['v'] = dataframe['v'].map(
                    lambda x: ','.join([str(ascs_family[k]) if k in ascs_family else k for k in x.split(',')]))
        else:
            dataframe['v'] = dataframe[v_field].map(unpack_genes)
            if use_asc == 'cluster':
                dataframe['v'] = dataframe['v'].map(
                    lambda x: ','.join([str(ascs_dict[k+'*01']) if k+'*01' in ascs_dict else k for k in x.split(',')]))
            elif use_asc == 'family':
                dataframe['v'] = dataframe['v'].map(
                    lambda x: ','.join([str(ascs_family[k+'*01']) if k+'*01' in ascs_family else k for k in x.split(',')]))
        array = dataframe[[sequence_id, 'cdr3_length', cdr3_field, 'v']].values
    else:
        array = dataframe[[sequence_id, 'cdr3_length', cdr3_field]].values
    if use_j:
        group.append('j')
        if allele:
            dataframe['j'] = dataframe[j_field]
        else:
            dataframe['j'] = dataframe[j_field].map(unpack_genes)
        array = dataframe[[sequence_id, 'cdr3_length', cdr3_field, 'v', 'j']].values
    output = defaultdict(lambda: defaultdict(set))
    for row in array:
        if use_v:
            for v in set(row[3].split(',')):
                if use_j:
                    for j in set(row[4].split(',')):
                        output[(v, j, row[1])][row[2]].add(row[0])
                else:
                    output[(v, row[1])][row[2]].add(row[0])
        else:
            output[(row[1])][row[2]].add(row[0])
    return output

def search_two_dbs(hash1: dict, hash2: dict, threshold: float = .7, score_func: Callable = ratio,
                   use_threshold: bool = True) -> DefaultDict:
    overlap = set(hash1.keys()).intersection(set(hash2.keys()))
    if use_threshold:
        matches = defaultdict(set)
    else:
        matches = defaultdict(lambda: defaultdict(int))
    for key in overlap:
        h3s_1 = hash1[key]
        h3s_2 = hash2[key]
        for h3_1 in h3s_1:
            for h3_2 in h3s_2:
                match = False
                if use_threshold:
                    if h3_1 == h3_2:
                        match = True
                    elif score_func(h3_1, h3_2) >= threshold:
                        match = True
                    if match:
                        for entry1 in hash1[key][h3_1]:
                            matches[entry1].update(hash2[key][h3_2])
                else:
                    score = score_func(h3_1, h3_2)
                    for entry1 in hash1[key][h3_1]:
                        for entry2 in hash2[key][h3_2]:
                            matches[entry1][entry2] = score
    return matches

def annotate_og_file(dataframe: pd.DataFrame, matches: DefaultDict, sequence_id: str = 'sequence_id', column_name: str = 'query',
                     annot_func: Callable = lambda x:','.join(x)) -> pd.DataFrame:
    dataframe[column_name] = dataframe[sequence_id].map(matches).map(annot_func)
    return dataframe

def handle_error(func: Callable) -> Callable:
    def handle(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            raise Exception
    return handle



