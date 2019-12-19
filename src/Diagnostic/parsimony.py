from ete3 import Tree
import random


def random_leaf(t):
  for leaf in t:
    leaf.add_features(state=[[0,1000],[1000,0]][random.sample([0,1],1)[0]])


def forward(node):
  if len(node.children)==0:
    return

  nbs=len(node.children[0].state)
  st=[0]*nbs
  for c in node.children:
    stc=c.state
    for i in range(nbs):
      st[i]+=min([stc[j]+[1,0][i==j] for j in range(nbs)])
  node.add_features(state=st)

  
def forward_tree(t):
  for n in t.traverse("postorder"):
    forward(n)


def parsimony(t, spfam1, spfam2):
  """Compute parsimony scory.
  t: species tree
  spfam1: list of species having genes in 1st family
  spfam2: list of species having genes in 2nd family
  """
  ## First assign states on leaves
  for leaf in t:
    leaf.add_features(state=[[1000,0][leaf.name in spfam1],[1000,0][leaf.name in spfam2]])
    if leaf.state==[1000,1000]:
      leaf.state=[0,0]

  forward_tree(t)
  return min(t.state)

def venn(spfam1,spfam2):
  """Compute [#(Fam1\Fam2), #(Fam2\Fam1), #(Fam1 & Fam2)]
  spfam1: list of species having genes in 1st family
  spfam2: list of species having genes in 2nd family
  """
  return list(map(len,[[x for x in spfam1 if not x in spfam2],[x for x in spfam2 if not x in spfam1],[x for x in spfam1 if x in spfam2]]))

#print(t.get_ascii(attributes=["name","state"]))
  

  
