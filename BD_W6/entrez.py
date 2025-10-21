#!/usr/bin/python3

#modified from https://entrezpy.readthedocs.io/en/master/functions/efetch_func.html

import entrezpy.esearch.esearcher
import entrezpy.log.logger
import entrezpy.efetch.efetcher

entrezpy.log.logger.set_level('WARN')

e = entrezpy.esearch.esearcher.Esearcher("entrezpy",
                                         "s2906787@bioinfmsc9.bio.ed.ac.uk",
                                         apikey=None,
                                         apikey_var=None,
                                         threads=None,
                                         qid=None)
analyzer_result = e.inquire({'db' : 'pubmed',
                      'term' : 'Mouse [orgn] and Nanog',
                      'retmax' : '20',
                      'rettype' : 'uilist'})
print(analyzer_result.result.count, analyzer_result.result.uids)



#just fecth a couple of the hits
e = entrezpy.efetch.efetcher.Efetcher("entrezpy",
                                      "s2906787@bioinfmsc9.bio.ed.ac.uk",
                                      apikey=None,
                                      apikey_var=None,
                                      threads=None,
                                      qid=None)
analyzer = e.inquire({'db' : 'pubmed',
                      'id' : [36243240,36213683],
                      'retmode' : 'text',
                      'rettype' : 'abstract'})
print(analyzer.get_result())

