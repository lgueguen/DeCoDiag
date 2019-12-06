#coding:utf8

def study_poly3(anc):
  c3=anc.get_cycles(3)
  graph=anc.dic_graph()
  d3={fam:{sp:[graph[sp].degree(gene) for gene in genes] for sp,genes in dicfam.items()} for fam,dicfam in c3.items()}

  # split according to nb of genes
  lc3=[filter(lambda x:len(set(x))==i,c3) for i in range(1,4)]
  ld3=[[d3[k] for k in c3] for c3 in lc3]
  return(ld3)

