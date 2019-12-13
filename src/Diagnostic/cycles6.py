#!/usr/bin/python3
#coding:utf8

import ete3
import glob
import sys

def get_par_fam(anc):  
  """ Get 6-cycles with 2 families present at least 2 times."""
  c6=anc.get_cycles(6)

  # 6-cycles with nb of occurences for each fam
  dlen6={lfam:{fam:len([v for v in lfam if v==fam]) for fam in lfam} for lfam in c6}
  # get 6-cycles with 2 families present at least 2 times
  d6_2_deg={lfam:min(v.values()) for lfam,v in d6.items() if len([g for g in dlen6[lfam].values() if g>=2])==2}  
  return d6_2_deg

def get_dup_fam(anc):  
  """ Get 6-cycles exactly one family present at least 4 times."""
  c6=anc.get_cycles(6)

  # 6-cycles with nb of occurences for each fam
  dlen6={lfam:{fam:len([v for v in lfam if v==fam]) for fam in lfam} for lfam in c6}
  # get 6-cycles exactly one family present at least 4 times
  d6_2_deg={lfam:min(v.values()) for lfam,v in d6.items() if len([g for g in dlen6[lfam].values() if g>=4])==1}  
  return d6_2_deg


def get_dup_fam_species(anc):  
  """Get families in 6-cycles with two occurences of a same gene (ie at
least 4 ends), associated with list of species where there is the
6-cycle."""
  
  c6=anc.get_cycles(6)

  # 6-cycles with {species:degree of gene} for each fam
  
  d6={fam:{sp:[anc.get_graph(sp).degree(gene) for gene in genes] for sp,genes in dicfam.items()} for fam,dicfam in c6.items()}

  # 6-cycles with nb of occurences for each fam
  dlen6={lfam:{fam:len([v for v in lfam if v==fam]) for fam in lfam} for lfam in d6}

  # 6-cycles with family with max nb of occurences
  dmaxlen6={lfam:max(flen,key=lambda x: flen[x]) for lfam,flen in dlen6.items()}
  dmax6={lfam:(m,dlen6[lfam][m]) for lfam, m in dmaxlen6.items()}

  # 6-cycles with two occurences of a same gene (ie at least 4 ends)
  d6_4_sp={dmax6[lfam][0]:list(v.keys())  for lfam,v in d6.items() if dmax6[lfam][1]>=4}

  return d6_4_sp

def simple_cycles(d6_2_deg):
  """ filter in cycles with pattern [3,2,2,3,2,2]."""
  d6_simp={k:v for k,v in d6_2_deg.items() if sorted(v)==[2,2,2,2,3,3]}
  
###########
### Output

def output(d6_2_deg):
  d6_2_count={lfam:{fam:len([w for w in lfam if w==fam]) for fam in lfam} for lfam in d6_2_deg}
  str6_2="\n".join(["\t".join(map(str,[w for w in count if count[w]>=2]))+"\t:\t"+"\t".join(map(str,[w for w in count if count[w]==1])) for count in d6_2_count.values()])
  return(str6_2)

# def output_cycles6(d6_2_deg):
#   d6_2_count={lfam:{fam:len([w for w in lfam if w==fam]) for fam in lfam} for lfam in d6_2_deg}
#   str6_2="\n".join(["\t".join(map(str,[w for w in count if count[w]>=2]))+"\t:\t"+"\t".join(map(str,[w for w in count if count[w]==1])) for count in d6_2_count.values()])
#   return(str6_2)

def output_species_dup_fam(d6_4_sp):
  str6_dup="\n".join([str(fam)+"\t"+"\t".join(map(str,sp)) for fam,sp in d6_4_sp.items()])
  return str6_dup


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

def zip_dup(dict_dup_fam, gene_tree_file, spec_tree):
  ftree=open(gene_tree_file,"r")
  itf=sorted(map(lambda x:x.strip(),ftree.readlines()))
  ftree.close()

  famnum=list(dict_dup_fam.keys())
  famnum.sort()
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
  

