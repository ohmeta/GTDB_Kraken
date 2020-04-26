import os,sys,re,logging
from Bio import SeqIO

class taxon_node:
    def __init__(self,taxid,sciName,rank=None,parent=None,children=[]):
        taxon_node._taxid=taxid
        taxon_node._sciName=sciName
        taxon_node._rank=rank
        taxon_node._parent=parent
        taxon_node._children=children
    
    def getName(self):
        return self._sciName
    
    def getChildren(self):
        return self._children
    
    def getParent(self):
        return self._parent


class taxon_tree:
    def __init__(self,rootnode)
        self._standard_ranks = ['species','genus','family','order','class','phylum','superkingdom']
        self._dic={1:[]}
        self._name=[]
        self._root=rootnode
        self._lastTaxid=1
    
    def IsEmpty(self):
        if self._root==None: return True
    
    def root(self):
        return self._root
    
    def setRoot(self,rootnode):
        self._root=rootnode
    def getRoot(self):
        return self._root
    
    def getNode(self,sciname):

    def lastTaxid(self):return self._lastTaxid

    def Isintree(self,sciname):
        if sciname in self._name: return True
        else: return False
    
    def addNode(self,TreeNode):



def lineage(lineage_str,t):
    taxon_line=lineage_str.split(';')
    parent_node=t.getRoot()

    for name in taxon_line:
        tmp=name.split('__')
        sciName=tmp[1]
        rank=tmp[0]
        if not t.Isintree(sciName):
            taxid=tax_count+1
            A=taxon_node(taxid,sciName,rank,parent_node)
        else: 
        parent_node=A
        
    


        


def addKrakenTXid(fasta,ktaxid,out):
    tmp=[]
    print('Transforming '+fasta+' to kraken format library fna\n',file=sys.stderr)
    for seq_rec in SeqIO.parse(fasta,'fasta'):
        seq_rec.id=seq_rec.id+'|kraken:taxid|'+ktaxid
        tmp.append(seq_rec)
    
    SeqIO.write(tmp,out,"fasta")
    tmp=[]

def 
