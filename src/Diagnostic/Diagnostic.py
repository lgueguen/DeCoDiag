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
import inputs 

import cycles6

## Main class     
class Diagnostic:
  def __init__(self, param_file):

    # dict of param files
    self.dfile=inputs.read_param(param_file)
    self.dfile["param.file"]=param_file.split(os.sep)[-1]
    
    # Get species_tree
    directory=self.dfile["output.dir"]
    if directory[-1]!=os.sep:
      directory+=os.sep

    prefix=self.dfile["output.prefix"]

    f=open(directory+prefix+".speciesTree.newick","r")
    l=f.readline()
    f.close()
    self.species_tree=Tree(l)

    # add feature number on nodes
    for l in self.species_tree.traverse("postorder"):
      name=l.name
      if name.find("@")!=-1: #leaf
        n2=name.split("@")
        l.name=n2[0]
        l.number=int(n2[1])
      else:
        l.number=int(l.support)

  ## Build ancestral graphs
  def build_ancestral(self):
    directory=self.dfile["output.dir"]
    if directory[-1]!=os.sep:
      directory+=os.sep
    prefix=self.dfile["output.prefix"]
    self.anc= ancestral.Ancestral(directory + prefix + ".genes.txt", self.species_tree, directory + prefix + ".adjacencies.txt")


  ## get 6-cycles without duplications to try clustering
  def output_broken_families_in_cycles6(self):
    c6br=cycles6.get_broken_families(self.anc)
    return cycles6.output_broken_families(c6br)

  ## Zip (one step) all cycle6 with duplications
  def zip_cycles6_dup(self):
    c6dup=cycles6.get_dup_families(self.anc)
    zipfam=cycles6.zip_dup(c6dup, self.dfile["gene.distribution.file"], self.species_tree)
    return(zipfam)

  ## output gene_trees
  ## input: dict of gene_trees, with same keys as in dfile "gene.distribution.file"
  ## build: same file names as in dfile, in a given directory
  ##        new file similar to content of dfile["gene.distribution.file"], but updated for new gene trees

  
  def output_gene_trees(self, gene_trees, new_directory, new_config_file):
    if not os.path.exists(new_directory):
      os.mkdir(new_directory)
    
    if new_directory[-1]!=os.sep:
      new_directory+=os.sep
    
    gene_tree_file=self.dfile["gene.distribution.file"]
    ftree=open(gene_tree_file,"r")
    itf=sorted(map(lambda x:x.strip(),ftree.readlines()))
    ftree.close()
    
    new_distrib_lines=[]
    keys=sorted(gene_trees.keys())
    ik=0
    
    for numtree in range(len(itf)):
      if ik<len(keys) and numtree==keys[ik]:
        new_file=new_directory+itf[numtree][itf[numtree].rfind(os.sep)+1:]
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
    
    pos=self.dfile[fname].rfind("_")
    if (pos==-1):
      return(self.dfile[fname]+"_1")
    else:
      suff=self.dfile[fname][pos+1:]
      if suff.isdigit():
        return(self.dfile[fname][:pos]+"_"+str(int(suff)+1))
      else:
        return(self.dfile[fname]+"_1")

      
  def new_param_file(self):
    self.dfile["param.file"]=self.increment_suffix_in_param_file("param.file")
    self.dfile["output.dir"]=self.increment_suffix_in_param_file("output.dir")
    paramf=open(self.dfile["param.file"],"w")
    paramf.write("\n".join([k+"="+v for k,v in self.dfile.items()])+"\n")
    paramf.close()



#def main():
if len(sys.argv)>=2:
  param_file=sys.argv[1]
else:
  Tk().withdraw() 
  param_file=askopenfilename()  

os.chdir(os.path.dirname(param_file))

diag=Diagnostic(param_file)
diag.build_ancestral()

print(diag.anc)


######################
## zip 6-cycles with duplications

choice=input('Build new trees ? (y/n) : ')
  
if choice=='y':                   
  new_directory=input('Name of the new directory for trees ? : ')
  # while os.path.exists(new_directory):
  #   print("This directory already exists.")
  #   new_directory=input('Name of the new directory for trees ? : ')
    
  zipfam=diag.zip_cycles6_dup()
  
  new_config_file=diag.increment_suffix_in_param_file("gene.distribution.file")
  diag.output_gene_trees(zipfam, new_directory, new_config_file)
  diag.dfile["gene.distribution.file"]=new_config_file
  diag.new_param_file()
  

######################
## get 6-cycles without duplications to try clustering

choice=input('Get parallel families in 6-cycles ? (y/n) : ')
  
if choice=='y':                   
  br_file=input('Name of the file for broken families ? : ')
  # while os.path.exists(br_file):
  #   print("This file already exists.")
  #   br_file=input('Name of the file for broken families ? : ')

  fbr=open(br_file,"w")
  fbr.write(diag.output_broken_families_in_cycles6())
  fbr.close()

  
# if __name__=="__main__":
#   anc=main()

