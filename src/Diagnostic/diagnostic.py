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
import cycles, parsimony


## Main class     
class Diagnostic:
  def __init__(self, param_file):

    # dict of param files
    self.dfile=iO.read_param(param_file)
    self.dfile["param.file"]=param_file.split(os.sep)[-1]
    
    # Get species_tree
    directory=self.dfile["output.dir"]
    prefix=self.dfile["output.prefix"]

    f=open(os.path.join(directory,prefix+".speciesTree.newick"),"r")
    l=f.readline()
    f.close()
    self.species_tree=Tree(l)

    ## Get gene tree files
    ftree=open(self.dfile["gene.distribution.file"],"r")
    self.gtree_files=sorted(map(lambda x:x.strip(),ftree.readlines()))
    ftree.close()

    # add feature number on nodes
    for l in self.species_tree.traverse("postorder"):
      name=l.name
      if name.find("@")!=-1: #leaf
        n2=name.split("@")
        l.name=n2[0]
        l.number=int(n2[1])
      else:
        l.number=int(l.support)

    ## Dictionnary of family trees, read when needed
    self.fam_trees={}

  ## Build ancestral graphs
  def build_ancestral(self):
    directory=self.dfile["output.dir"]
    prefix=self.dfile["output.prefix"]
    self.anc= ancestral.Ancestral(os.path.join(directory,prefix + ".genes.txt"), self.species_tree, os.path.join(directory, prefix + ".adjacencies.txt"))

# #def main():
# if len(sys.argv)>=3:
#   param_file=os.path.join(sys.argv[1],sys.argv[2])
# elif len(sys.argv)>=2:
#   param_file=sys.argv[1]
# else:
#   Tk().withdraw() 
#   param_file=askopenfilename()  

# if len(os.path.dirname(param_file))!=0:
#   os.chdir(os.path.dirname(param_file))

# diag=Diagnostic(os.path.basename(param_file))

  def compute_parsimony(self, fam1, fam2):
    """Compute the parsimony score of switch in gene presence between
    fam1 & fam2, along the species tree.

    param fam1, fam2: numbers of studied families

    """ 
    ## Look for gene trees for these families
    for f in [fam1,fam2]:
      if not f in self.fam_trees:
        self.fam_trees[f]=Tree(self.gtree_files[f])
    
    ## Look for extant species in each tree
    lv={}
    for f in [fam1,fam2]:
      lv[f]=[l[:l.find("@")] for l in self.fam_trees[f].get_leaf_names()]
    
    return parsimony.parsimony(self.species_tree,lv[fam1],lv[fam2])
  
  def compute_venn(self, fam1, fam2):
    """Compute the Venn scores of gene presence between fam1 & fam2,.

    param fam1, fam2: numbers of studied families

    """ 
    ## Look for gene trees for these families
    for f in [fam1,fam2]:
      if not f in self.fam_trees:
        self.fam_trees[f]=Tree(self.gtree_files[f])
    
    ## Look for extant species in each tree
    lv={}
    for f in [fam1,fam2]:
      lv[f]=[l[:l.find("@")] for l in self.fam_trees[f].get_leaf_names()]
    
    return parsimony.venn(lv[fam1],lv[fam2])
  

  def output_cycles(self):

    dup_file= input(" Name of the file for families involved on 4-cycles? ").strip()

    if (len(dup_file)!=0):
      d2b=cycles.filter(self.anc,[4],1,1)
      fdup=open(dup_file,"w")
      fdup.write(iO.output_genes_sp(self.anc.get_cycles(4), d2b[4]))
      fdup.close()

    ######################
    ## get cycles with duplicated genes

    dup_file=input(' Name of the file for duplicated families ? : ').strip()
    if (len(dup_file)!=0):
      ddup=cycles.filter(self.anc,self.anc.cycle_lengths(),2,1)
      fdup=open(dup_file,"w")
      for k,v in sorted(ddup.items()):
        if k>10:
          break
      fdup.write(iO.output_genes_sp(self.anc.get_cycles(k), v))
      fdup.close()

    ######################
    ## get parallel families to try clustering

    br_file=input(' Name of the file for parallel families ? : ').strip()
    # while os.path.exists(br_file):
    #   print("This file already exists.")
    #   br_file=input('Name of the file for broken families ? : ')

    if (len(br_file)!=0):
      dpar=cycles.filter(self.anc,self.anc.cycle_lengths(),1,2)
      fbr=open(br_file,"w")
      for k,v in sorted(dpar.items()):
        if k>10:
          break
        fbr.write(iO.output_genes_sp(self.anc.get_cycles(k), v))
      fbr.close()

    ######################
    ## get parsimony and venn scores for parallel families

    pars_file=input(' Name of the file for parsimony scores ? : ').strip()
    # while os.path.exists(br_file):
    #   print("This file already exists.")
    #   br_file=input('Name of the file for broken families ? : ')

    if (len(pars_file)!=0):
      if not "dpar" in locals():
        dpar=cycles.filter(self.anc,self.anc.cycle_lengths(),1,2)

      dsc={}
      for k,lv in sorted(dpar.items()):
        if k>10:
          break

        for v in lv:
          for f1 in v:
            for f2 in v:
              if f2<=f1 or (f1,f2) in dsc:
                continue
              
              p=self.compute_parsimony(f1, f2)
              v=self.compute_venn(f1, f2)
              dsc[(f1,f2)]=[p]+v

              
      fbr=open(pars_file,"w")
      fbr.write("\t".join(["Fam1","Fam2","Pars","N1.2","N2.1","N12"])+"\n")
      for k,pk in sorted(dsc.items()):
        fbr.write(iO.joinstr("\t",list(k)+pk)+"\n")

      fbr.close()

  

  
# if __name__=="__main__":
#   anc=main()

