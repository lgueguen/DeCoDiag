#!/usr/bin/python3
#coding:utf8

import ete3
import glob
import sys, functools


def filter(anc, lsize, nbgenes, nbfam):
  """Filter in the size-cycles with at least nbfam families present
  complete nbgenes times.

  anc: Ancestral object 
  lsize: List of sizes of the considered cycles
  nbgenes: Minimal number of complete genes present for any family
  nbfam: Minimal number of families with given condition on the number
         of genes.

  Output: {size of cycles: {cycle description: list of families fulfilling the conditions}}
  """

  dd={}
  for size in lsize:
    c=anc.get_cycles(size)
    famnum={lfam:{x:lfam.count(x) for x in set(lfam)} for lfam in c.keys()}
    d2b0={lfam:[g for g in lg if lg[g]>=2*nbgenes] for lfam,lg in famnum.items()}
    d2b1={lfam:lg for lfam,lg in d2b0.items() if len(lg)>=nbfam}
    dd[size]=d2b1
    
  return(dd)


# def simple_cycles(d6_2_deg):
#   """ Filter in cycles with sorted pattern [2,2,2,2,3,3]."""
#   c6=anc.get_cycles(6)
#   d6_simp0={k:[sp for sp,ld in v.deg.items() if  sorted(ld)==[2,2,2,2,3,3]] for k,v in c6.items()}
#   d6_simp={lfam:ld for lfam,ld in d6_simp0.items() if len(ld)!=0}
#   return(d6_simp)

#########################################
#### MGMT OF DUPLICATIONS


## check if node in gene_tree is a duplication node
def isDupNode(dupNode, spec_tree):
  lmrca=[spec_tree.get_common_ancestor([l[:l.find("@")] for l in child.get_leaf_names()]) for child in dupNode.children]
  for i in range(len(lmrca)):
    for j in range(i+1,len(lmrca)):
      if lmrca[i]==lmrca[j] or lmrca[i] in lmrca[j]:
        return True
  return False


## Zip one duplication (ie one step towards the leaves)
def zip1dup(geneTree, dup_nodes, spec_tree):
  alln=geneTree.get_leaf_names()
  Tn=max(dup_nodes)
  for node in spec_tree.traverse():
    if node.number==Tn:
      specNode=node
      break
  
  dupNode=geneTree.get_common_ancestor([l for l in alln if l[:l.find("@")] in specNode.get_leaf_names()])
  
  #look for just above duplication node
  while not isDupNode(dupNode, spec_tree) and dupNode.up:    
    dupNode=dupNode.up
  
  ln=dupNode.get_leaf_names()
  # node in species tree at same location
  specNode=spec_tree.get_common_ancestor([l[:l.find("@")] for l in ln])
  
  ## build return tree
  headTree=geneTree.copy()
    
  # get dupnode in result tree
  rootdup=headTree.get_common_ancestor(ln)
  if rootdup.up:  #rootdup not at the root
    rootdup.add_child(name="Dup")
    dup=headTree.search_nodes(name="Dup")[0]
    headTree.prune([l for l in alln if not l in ln]+["Dup"])
  else:
    headTree=ete3.Tree()
    dup=headTree
    
  for child in specNode.children:
    lf=[l for l in ln if l[:l.find("@")] in child.get_leaf_names()]
    if (lf!=[]):
        copdupTree=dupNode.copy()
        copdupTree.prune(lf)
        if len(lf)==1: # To avoid node with 1 son
          copdupTree=copdupTree.get_leaves()[0].copy()
    
    dup.add_child(copdupTree)
        
  return(headTree)
  
# tree_rep="../instance_test_Docker/DeCoSTAR_Anopheles_Xtopo_gene_trees"
# spec_tree_nf="../instance_test_Docker/Anopheles_species_tree_X_topology.nwk"

## get repertory of gene trees (to avoid loading everything)
## and an output of get_dup_fam_species {fam number: ids of species node to zip}
##
## return dict {family number:zipped ete3 gene Tree}

## dict_dup_fam: lfam in which a gene is duplicated
def zip_dup(dict_dup_fam, itf, spec_tree):

  ## get list of duplicated genes
  famnum=set(functools.reduce(lambda a,b:a+b, dict_dup_fam.values()))
  dupTrees={}
  dupNames={}
  compt=0
  nbi=0
  for tf in itf:
    if compt==famnum[nbi]:
      dupTrees[compt]=ete3.Tree(tf)
      dupNames[compt]=tf
      nbi+=1
      if nbi==len(famnum):
        break
    compt+=1

  dcorTree={}
  for dTk in dupTrees:
    try:
      dcorTree[dTk]=zip1dup(dupTrees[dTk], dict_dup_fam[dTk], spec_tree)
    except:
      print(dTk)
      print(dupTrees[dTk])
      print(dict_dup_fam[dTk])
      print(lnodes)
      
  return(dcorTree)
  

