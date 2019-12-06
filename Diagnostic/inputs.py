#!/usr/bin/python3
#coding:utf8

import re
import os                      
import ete3
import glob

def find_family(gene, dic_gene):
  """
  from (gene,end) name to (family,end)
  """
  g0=gene[0]
  if '|'  in g0:                          
    return int(g0.split('|')[0])
  return int(dic_gene[g0])

def find_families(lgene, dic_gene):
  """
  from list of (gene,end) name to list of (family,end)
  """
  g=[[find_family(gene, dic_gene),gene] for gene in lgene]
  lg=len(g)
  mg=min(x[0] for x in g)
  indm=[i for i in range(lg) if g[i][0]==mg] # minimum values
                                                # positions
  indv=set([i+1 for i in indm if i<lg-1] + [i-1 for i in indm if i>0] + [0]*(lg-1 in indm) + [lg-1]*(0 in indm)) #neighbour positions
  Mv=max(g[i][0] for i in indv)
  liMv=[i for i in range(lg) if g[i][0]==Mv] # index of maximum neighbour
  
  for iMv in liMv:
    if (g[iMv-1][0]==mg): # g must be reversed    
      vind=[iMv-1-i for i in range(lg)]
      break
    elif (iMv<lg-1 and g[iMv+1][0]==mg) or (iMv==lg-1 and g[0][0]==mg): #iMv+1 is the starting min, and g in correct orientation
      vind=[i-(iMv+1) for i in range(lg)]
      break
    
  gs=[0]*lg
  try:
    for i in range(lg):  # min element first
      gs[vind[i]]=g[i]
  except:
    print(g)

  return (tuple([s[0] for s in gs]),tuple([s[1] for s in gs]))
  
      
def read_param (param_file) :
  """ Parsing of parameter file of DeCoSTAR. """

  param={}
  
  f=open(param_file,"r")
  for l in f.readlines():
    if l.find("=")!=-1:
      line = l.split("=")
      param[line[0].strip()]=line[1].strip()

  f.close()

  required=["species.file","output.dir","output.prefix"]

  for req in required:
    if not req in param.keys():
      raise Exception("Can't find " + req + " in " + param_file)

  return(param)

  
