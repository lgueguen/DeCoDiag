#!/usr/bin/python3
#coding:utf8

import re
import networkx as nx
import time
import os, shutil
import sys
from ete3 import Tree  

from tkinter import Tk
from tkinter.filedialog import askopenfilename

import ancestral 
import iO

import cycles


## Main class     
class Fix_bubbles:
  def __init__(self, diagnostic):

    ## Diagnostic that has all information
    self.diag=diagnostic

    ## list of tree files
    self.__lfiles_gene_tree = self.diag.get_files_gene_tree()

    ## dictionnary of family trees shared with diagnostic
    self.fam_trees = self.diag.fam_trees

    ## species trees shared with diagnostic
    self.species_tree=self.diag.species_tree

    
  ##################################################
  #### Functions for zip
  
  ## Zip (one step) all duplicated families in cycles.
  ##
  ## ln list of cycle sizes
  
  def zip_cycles_dup(self, ln):
    cndup=cycles.filter(self.diag.anc, ln, 2, 1)

    cdup={}
    for ndup in cndup.values():
      cdup.update(ndup)

    ## associate duplicated family and associated species 
    csp={}
    for n in ln:
      cyN=self.diag.anc.get_cycles(n)
      for lfam,ldup in cdup.items():
        if lfam in cyN:
          lsp=list(cyN[lfam].keys())
          for fam in ldup:
            csp[fam]=lsp
        
    zipfam=cycles.zip_dup(csp, self.fam_trees, self.__lfiles_gene_tree, self.species_tree)
    return(zipfam)

  
  def output_gene_trees(self, gene_trees, new_directory, new_config_file):
    """Output gene_trees.
    
    gene_trees: dict of gene_trees, with same keys as in dfile "gene.distribution.file"
    new_directory: directory with new gene tree files
    new_config_file:  file similar to content of dfile["gene.distribution.file"], but updated for new gene trees.
    """

    if not os.path.exists(new_directory):
      os.mkdir(new_directory)
    
    gene_tree_file=self.diag.dfile["gene.distribution.file"]
    ftree=open(gene_tree_file,"r")
    itf=list(map(lambda x:x.strip(),ftree.readlines()))
    ftree.close()
    
    new_distrib_lines=[]
    keys=sorted(gene_trees.keys())
    ik=0
    
    for numtree in range(len(itf)):
      if ik<len(keys) and numtree==keys[ik]:
        new_file=os.path.join(new_directory,itf[numtree][itf[numtree].rfind(os.sep)+1:])        
        gene_trees[numtree].write(format=9,outfile=new_file)
        ik+=1
      else:
        new_file=itf[numtree]
        #shutil.copyfile(os.getcwd()+os.sep+itf[numtree],new_file)
      
      new_distrib_lines.append(new_file)

    f=open(new_config_file,"w")
    f.write("\n".join(new_distrib_lines)+"\n")
    f.close()


  def increment_suffix_in_param_file(self, fname):
    # increment suffix of fname in param file from dfile    
    pos=self.diag.dfile[fname].rfind("_")
    if (pos==-1):
      return(self.diag.dfile[fname]+"_1")
    else:
      suff=self.diag.dfile[fname][pos+1:]
      if suff.isdigit():
        return(self.diag.dfile[fname][:pos]+"_"+str(int(suff)+1))
      else:
        return(self.diag.dfile[fname]+"_1")

      
  def new_param_file(self):
    dfile=self.diag.dfile
    for df in self.diag.dfile:
      if df in ["param.file","output.dir", 'gene.distribution.file']:
        dfile[df]=self.increment_suffix_in_param_file(df)

    paramf=open(dfile["param.file"],"w")
    paramf.write("\n".join([k+"="+v for k,v in self.diag.dfile.items()])+"\n")
    paramf.close()


