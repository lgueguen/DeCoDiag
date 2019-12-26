#!/usr/bin/python3
#coding:utf8

import sys
import re
import networkx as nx

import iO
import functools

#  Categorization of cycle, for a given set of families
#  Inherits from dict as (species,lists of genes)

class Cycle(dict):
  def __init__(self, fam, graph):
    self.fam=fam
    self.graph=graph
    self.compute_degree()
    
  def compute_degree(self):
    """ computes degree of genes for all species."""
    self.deg={sp:[self.graph[sp].degree(gene) for gene in genes] for sp,genes in self.items()}

  
    
def check_singleton(graph, gene):
  """
  Identification if is singleton
  Return list of lists of nodes involved
  """
  lok=[]
  neighbors = graph.neighbors(gene)
  for neigh in neighbors:
    if graph.degree(neigh) == 2:
      for k in graph.neighbors(neigh):
        if k != gene and graph.degree(k) == 1:
          lok.append([neigh, k])
          
  return lok

def check_OEN(graph, gene):
  """
  Identification if is OEN
  Return list of lists of nodes involved
  """
  lok=[]
  neighbors = graph.neighbors(gene)
  for neigh in neighbors:
    if graph.degree(neigh) == 1:
      lok.append([neigh])

  return lok

class Ancestral:
  """
  All about ancestral adjacencies
  """
  
  def __init__(self, gene_file, species_tree, adj_file):

    self.species_tree=species_tree    

    # dict of genes
    self.__dic_gene = self.__get_dic_gene(gene_file)

    # dict of adj graphs
    self.__dic_graph = {node.number:nx.Graph() for node in species_tree.traverse("postorder")}
    self.__read_genes(gene_file)
    self.__read_adj(adj_file)

    # look for conflicts
    self.__conf=self.list_conflicts()
    self.OEN = self.class_conflict(check_OEN)
    self.sing = self.class_conflict(check_singleton)
    self.__cycles ={}
    self.find_cycles()

    
    
  #############################
  # Elements
  #############################

  def dic_graph(self):
    return self.__dic_graph

  def dic_gene(self):
    return self.__dic_gene

  def get_graph(self, sp):
    return self.dic_graph()[sp]

  def __read_genes(self, gene_file):
    f = open(gene_file, "r")
    for i in f.readlines():
      line = i.split()
      self.__dic_graph[int(line[0])].add_node((line[1],'start'))
      self.__dic_graph[int(line[0])].add_node((line[1],'stop'))
      self.__dic_graph[int(line[0])].add_edge((line[1],'start'),(line[1],'stop'),eq_class="gene", edge_type="gene", weight = 1)
    f.close()
    if len(self.__dic_graph)==0:
      raise Exception("No genes in " + gene_file)
    
  def __get_dic_gene(self, fgenes):
    """
    Dictionnary of genes as {leaf label : gene family number}
    """
    gene=open(fgenes,"r")
    dic_gene={}
    for line in gene.readlines():
      ll=line.split()
      if len(ll)>2:
        fam=ll[1].split("|")[0]
        for x in ll[2:]:
          if '@' in x:
            dic_gene[x]=fam
      
    gene.close()
    return dic_gene

  ###############################
  # Identification of conflicts
  ##############################

  def nb_conflict(self): 
    """
    Return number of genes with conflicts.
    """
    graph=self.__dic_graph
    nbconflit=0
    for species in graph.keys():
      deg = list(graph[species].degree().values())
      nbconflit+=functools.reduce(lambda n,d: n+(d>2), deg, 0)

    return nbconflit

  def list_conflicts(self): 
    """
    Return dictionary {species: list of genes with conflicts}.
    """
    graph=self.__dic_graph
    nbconflit=0
    dconf={}
    for species in graph.keys():
      deg = graph[species].degree()
      dconf[species]=[d for d,v in deg.items() if v>2]
      
    return dconf

  def family_numbers(self):
    """ Return the ordered list of family numbers."""
    return list(set(map(int,self.dic_gene.values())))
  
  ###############################
  # Conflicts categorization
  ##############################

  def class_conflict(self,function):   
    dicsp={}
    dicg={}
    for species in self.__dic_graph:
      ds={}
      graph= self.__dic_graph[species]
      conf = self.__conf[species]

      lseen={}
      for node in conf:
        if node in lseen:
          continue
        
        sing=function(graph,node)
      
        if len(sing)!=0:
          family=iO.find_family(node,self.__dic_gene)
          if family in dicsp.keys():
            if not species in dicsp[family]:
              dicsp[family].append(species)
              dicg[family].append(tuple([node]+sing[0]))
          else:           
            dicsp[family]=[species]
            dicg[family]=[tuple([node]+sing[0])]
          for ls in sing:
            for gene in ls:
              lseen[gene]=""

    return (dicsp,dicg)
    
  def find_cycles(self):
    """
    Find cycles in graphs, sort them by length, and by sets of gene families
    """
    cycles={}
    graph=self.__dic_graph
    for sp,gr in graph.items():
      cycles[sp]=nx.cycle_basis(gr)
      
    for sp in graph.keys():
      if len(cycles[sp])==0:
        del cycles[sp]
    
    ## cycle lengths
    lsp={sp:[len(x) for x in cysp] for sp,cysp in cycles.items()}
    lengths=set(functools.reduce(lambda x,y:x+y, lsp.values()))
    
    self.__cycles={}
    for leng in lengths:
      dtri0={sp:[x for x in cysp if len(x)==leng] for sp,cysp in cycles.items()}
      dtri={sp:l for sp,l in dtri0.items() if len(l)>0}
      dtrifam={sp:[iO.find_families(tri,self.dic_gene()) for tri in ltri] for sp,ltri in dtri.items()}
      fr = functools.reduce(lambda x,y:x+y, dtrifam.values())
      dfam=set([x[0] for x in fr])
      
      self.__cycles[leng]={}
      cycleng=self.__cycles[leng]
      
      for fam in dfam:
        cycleng[fam]=Cycle(fam, self.dic_graph())
      
        for sp,vsp in dtrifam.items():
          for v in vsp:
            if v[0]==fam:
              cycleng[fam][sp]=v[1]
              break
  

  def cycle_lengths(self):
    return self.__cycles.keys()

  def get_cycles(self, length):
    return self.__cycles[length]

  def __str__(self):
    dres=[]
    dres.append(["Conflicts",self.nb_conflict()])
      
    dres.append(["Unique singletons",len(self.sing[0])])
    dres.append(["Singletons",sum(len(k) for k in self.sing[0].values())])
    
    dres.append(["Unique OEN",len(self.OEN[0])])
    dres.append(["OEN",sum(len(k) for k in self.OEN[0].values())])

    for lc in sorted(self.cycle_lengths()):
      dres.append(["cycle %d"%lc,len(self.get_cycles(lc))])
    
    return("\n".join([k+" \t"+str(v) for [k,v] in dres]))

  #Parsing of DeCoSTAR files

  def __read_adj(self,adj_file):
    """
    Parsing of ancestral adjacencies file
    """
    
    f=open(adj_file,'r') 
    for i in f.readlines():
      line = i.split()
      ## add edge if proba > 0
      if float(line[6]) == 0:
        continue

      deb = (line[1],['stop','start'][line[3] == "-"])
      fin = (line[2],['stop','start'][line[4] == "+"])

      self.__dic_graph[int(line[0])].add_edge(deb,fin,weight = float(line[6]))
      
    f.close()

