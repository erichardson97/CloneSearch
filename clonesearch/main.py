#!/usr/bin/env python3
from .utils import *
from typing import Callable
# from data_file import iedb, ebola

class AirrFile():
    def __init__(self, path: str, use_v: bool = True, use_j: bool = False, use_allele: bool = False):
        self.file = pd.read_csv(path, sep = '\t')
        self.use_v = use_v
        self.use_j = use_j
        self.use_allele = use_allele
        self.dbs = {}
        self.__init__ = handle_error(self.__init__)
        self.process = handle_error(self.process)
        self.add_db = handle_error(self.add_db)
        self.query_db = handle_error(self.query_db)
    def rename(self, rename: dict):
        self.file.rename(columns = rename, inplace=True)
    def process(self):
        self.master = make_hash(self.file, use_v = self.use_v, use_j = self.use_j, allele = self.use_allele)
    def add_db(self, path: str, name: str):
        self.dbs[name] = make_hash(pd.read_csv(path, sep = '\t'), use_v = self.use_v, use_j = self.use_j, allele = self.use_allele)
    def query_db(self, name: str, column_name: str, threshold: float = 0.7, function: Callable = ratio,
                 use_threshold: bool = True, annot_func: Callable = lambda x:','.join(x)):
        if name not in self.dbs:
            return logging.error(f"The database name you provided for query func doesn't exist, run add_db before.")
        else:
            hits = search_two_dbs(self.master, self.dbs[name], threshold = threshold,
                                  score_func = function, use_threshold = use_threshold)
            self.file = annotate_og_file(self.file, hits, column_name = column_name, annot_func=annot_func)

class ObjectSet():
    def __init__(self, path: str, db: list, names: list, use_v: bool = True, use_j: bool = False):
        airr_f = AirrFile(path, use_v = use_v, use_j = use_j)
        airr_f.process()
        for db_path, db_name in zip(db, names):
            if os.path.exists(db_path) is False:
                return 'No database found.'
            airr_f.add_db(db_path, db_name)
        self.airr_f = airr_f
    def query_db(self, name: str, column_name: str, threshold: float = 0.7, function: Callable = ratio):
        self.airr_f.query_db(name, column_name, threshold, function)
    def save(self, path: str, delim: str = '\t'):
        self.airr_f.file.to_csv(path, sep = delim)
