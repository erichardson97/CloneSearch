#!/usr/bin/env python3
from Levenshtein import ratio
from collections import defaultdict
from typing import DefaultDict, Callable
import pandas as pd
import logging
import requests
from functools import wraps

def download_covabdab(outf: str, filter_human: True) -> str:
    f = requests.get('https://opig.stats.ox.ac.uk/webapps/covabdab/static/downloads/CoV-AbDab_130623.csv')
    with open(outf, 'w') as k:
        k.write(f.content.decode())
    db = pd.read_csv(outf)
    if filter_human:
        db = db[(db['Heavy V Gene'].str.contains('Human'))&(db["Origin"].str.contains('Convalescent|Vaccinee|Human'))]
    db['v_call'] = db['Heavy V Gene'].map(lambda x:x.split('(')[0])
    db['j_call'] = db['Heavy J Gene'].map(lambda x:x.split('(')[0])
    hash = make_hash(db, cdr3_field = 'CDRH3', allele = False, sequence_id = 'Name')
    return hash

def unpack_genes(v_field: str) -> str:
    return ','.join(set([p.split('*')[0] for p in v_field.split(',')]))

def make_hash(dataframe: pd.DataFrame, v_field: str = 'v_call', cdr3_field: str = 'cdr3_aa',
              allele: bool = True, sequence_id: str = 'sequence_id', use_v: bool = True,
              use_j: bool = False, j_field: str = 'j_call') -> DefaultDict:
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
        else:
            dataframe['v'] = dataframe[v_field].map(unpack_genes)
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
    @wraps(func)
    def handle(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in: {func.__name__}: {e}")
            raise
    return handle
