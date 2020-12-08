import sys,os,re,argparse,logging,time
from Bio import SeqIO
import multiprocessing

## Build taxonomy from the beginning.
## GTDBtk result based
def parse_gtdb_result(outdir,gtdb_result):
    fasta2id={}
    taxdir=os.path.join(outdir,"taxonomy")
    if not os.path.exists(taxdir):
        os.mkdir(taxdir)
    node_out=open(os.path.join(taxdir,"nodes.dmp"),'w')
    name_out=open(os.path.join(taxdir,"names.dmp"),'w')
    rank={'d':'domain','p':'phylum','c':"class",'o':'order','f':'family','g':'genus','s':'species'}
    taxon_dict={'root':1}
    node.write('\t|\t'.join([str(1),str(1),'no rank','',8,0,1,0,0,0,1,0,1,'']+'\t|\n'))
    taxon_count=1
    with open(gtdb_result,'r') as handle:
        next(handle)
        for line in handle.readlines():
            line=line.strip()
            arr=line.split('\t')
            fastaid=arr[0]
            if re.search(r'GCF|GCA',fastaid):
                fastaid=re.sub(r'\.\d_ASM.*','',fastaid)
            taxonomy=arr[1].split(';')
            for tax in taxonomy:
                if len(tax)==3:
                    if re.match(r'^s',tax):
                        name_parent=re.sub('\w__','',parent)
                        spid=re.sub(r'GCA_|GCF_','',fastaid)
                        tax+=name_parent+' sp.'+spid
                    else:
                        tax+=fastaid
                rans,sciname=tax.split('__')
                if rans == 'd': parent='root'
                #print(tax)
                if tax not in taxon_dict:
                    taxon_count+=1
                    taxid=taxon_count
                    tax_node=taxid
                    node_out.write('\t|\t'.join([str(taxid),str(taxon_dict[parent]),rank[rans],'','0', '1', '11', '1', '0', '1', '1', tax])+'\t|\n')
                    name_out.write('\t|\t'.join([str(taxid),sciname,'','scientific name'])+'\t|\n')
                    taxon_dict[tax]=taxid
                parent=tax
            fasta2id[fastaid]=taxon_dict[tax]
    node_out.close()
    name_out.close()
    return fasta2id

def parsingGTDBtax(outdir,gtdbtax_file): #Building a nodes.dmp and names.dmp based on the GTDB taxonomy File
    taxdir=os.path.join(outdir,"taxonomy")
    if not os.path.exists(taxdir):
        os.mkdir(taxdir)
    rank={'d':'domain','p':'phylum','c':"class",'o':'order','f':'family','g':'genus','s':'species'}
    node_out=open(os.path.join(taxdir,"nodes.dmp"),'w')
    name_out=open(os.path.join(taxdir,"names.dmp"),'w')
    taxon_dict={'root':1}
    node_out.write('\t|\t'.join([str(1),str(1),'no rank','',8,0,1,0,0,0,1,0,1,'']+'\t|\n'))
    asm2id={}
    taxon_count=1
    with open(gtdbtax_file,'r') as handle:
        for line in handle.readlines():
            line=line.strip()
            arr=line.split('\t')
            asm=arr[0]
            asm=re.sub(r'(RS_|GB_)','',asm)
            tax_arr=arr[1].split(';')
            for tax in tax_arr:
                rans,sciname=tax_arr[i].split('__')
                if rans == 'd': parent='root'
                if tax not in taxon_dict:
                    taxon_count+=1
                    taxid=taxon_count
                    tax_node=taxid
                    node_out.write('\t|\t'.join([str(taxid),str(taxon_dict[parent]),rank[rans],'','0', '1', '11', '1', '0', '1', '1', tax]),'\t|\n')
                    name_out.write('\t|\t'.join([str(taxid),sciname,'','scientific name']))
                    taxon_dict[tax]=taxid
                parent=tax
            asm2id[asm]=taxon_dict[tax]
    return asm2id

def addKraken(outdir,fna_file,ass2id,mode): #Transforming the input fna in a Kraken2-readable fna
    #all_file=os.path.join(outdir,'library.fna')
    file_name=os.path.basename(fna_file)
    if re.search(r'gz',file_name): #Support Gzipped file
        handle=gzip.open(fna_file,'rt')
    else: 
        handle=open(fna_file,'r')
    
    if mode=='gtdbtax': #GTDB Database itself
        asmid=re.sub(r'_genomic\.fna.*','',file_name)
    elif mode='gtdbres': #GTDBtk result
        asmid=res.sub(r'\.\d_ASM.*|\.fa','',file_name)
    
    if asmid in ass2id:
        logging.info('Parsing '+fna_file+'...')
        arr=[]
        for seq_rec in SeqIO.parse(handle,"fasta"):
            seq_rec.id=seq_rec.id+'|kraken:taxid|'+taxid #This is still able to be recognized by the scan_fasta_file.pl
            arr.append(seq_rec)
        SeqIO.write(arr,os.path.join(outdir,file_name),"fasta")
        #SeqIO.write(arr,all_file,'fasta') #To output to a general library.fna
        logging.info('Finished\n')
    else:
        logging.error(asmid+' is not in taxon list\n')
    handle.close()

#def addKraken(outdir,fna_file,ass2id): #Transforming the input fna in a Kraken2-readable fna
#    file_name=os.path.basename(fna_file)
#    asmid=re.sub(r'\.\d_ASM.*|.fa','',file_name)
#    handle=open(fna_file,'r')
#    if asmid in ass2id:
#        logging.info('Parsing '+fna_file+'...')
#        arr=[]
#        taxid=ass2id[asmid]
#        for seq_rec in SeqIO.parse(handle,"fasta"):
#            seq_rec.id=seq_rec.id+'|kraken:taxid|'+taxid 
#            arr.append(seq_rec)
#        SeqIO.write(arr,os.path.join(outdir,file_name),"fasta") #To output to a general library.fna
#        logging.info('Finished\n')
#    else:
#        logging.error(asmid+' is not in taxon list\n')
#    handle.close()
            
