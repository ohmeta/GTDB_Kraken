import sys,os,re,argparse,logging
import multiprocessing
from Bio import SeqIO

def parseAsm2ID(asm2id): #Parse Assembly ID to Taxid List
    ass2id={}
    with open(asm2id,'r') as handle:
        for line in handle.readlines():
            line=line.strip()
            arr=line.split('\t')
            asmid=arr[0]
            asmid=re.sub(r'(RS_|GB_)','',asmid)
            taxid=arr[1]
            ass2id[asmid]=taxid
    return ass2id

def parsingGTDBtax(outdir,gtdbtax_file):
    taxdir=os.path.join(outdir,"taxonomy")
    if not os.path.exists(taxdir):
        os.mkdir(taxdir)
    rank={'d':'domain','p':'phylum','c':"class",'o':'order','f':'family','g':'genus','s':'species'}
    node_out=open(os.path.join(taxdir,"nodes.dmp"),'w')
    name_out=open(os.path.join(taxdir,"names.dmp"),'w')
    taxon_dict={'root':1}
    asm2id={}
    taxon_count=1
    with open(gtdbtax_file,'r') as handle:
        for line in handle.readlines():
            line=line.strip()
            arr=line.split('\s+')
            asm=arr[0]
            tax_arr=arr[1].split(';')
            for tax in tax_arr:
                rans,sciname=tax_arr[i].split('__')
                if rans == 'd': parent='root'
                if tax not in taxon_dict:
                    taxon_count+=1
                    taxid=taxon_count
                    tax_node=taxid
                    node_out.write('\t|\t'.join([str(taxid),str(taxon_dict[parent]),rank[rans],'','0', '1', '11', '1', '0', '1', '1', '0']),'\t|\n')
                    name_out.write('\t|\t'.join([str(taxid),sciname,'','scientific name']))
                    taxon_dict[tax]=taxid
                parent=tax
            asm2id[asm]=taxon_dict[tax]
    return asm2id

def addKraken(outdir,fna_file,ass2id):
    file_name=os.path.basename(fna_file)
    asmid=re.sub(r'_genomic.fna','',file_name)
    if asmid in ass2id:
        logging.info('Parsing '+fna_file+'...')
        arr=[]
        for seq_rec in SeqIO.parse(fna_file,"fasta"):
            seq_rec.id=seq_rec.id+'|kraken:taxid|'+taxid
            arr.append(seq_rec)
        SeqIO.write(arr,os.path.join(outdir,file_name),"fasta")
        logging.info('Finished\n')
    else:
        logging.error(asmid+' is not in taxon list\n')

    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--fna_dir",help='Directory to store')
    parser.add_argument('--taxid',help='Assembly ID to Taxid file')
    parser.add_argument('--outdir',help='Output Directory')
    args=parser.parse_args()
    
    logging.basicConfig(filename='./genome2kraken.log', level=logging.INFO)
    #Parse Assembly ID to Taxid List
    logging.info('Starting to Parse Assembly2Taxid File\n')
    ass2id=parseAsm2ID(args.taxid)
    logging.info('Finished\n')

    #Add taxid to kraken lib fna
    t_start=time.time()
    pool=multiprocessing.Pool(10)

    for fna_file in os.listdir(args.fna_dir):
        fna_abs=os.path.abspath(os.path.join(args.fna_dir,fna_file))
        pool.apply_async(run_ak,(args.outdir,fna_abs,ass2id))
    
    pool.close()
    pool.join()
    t_end=time.time()
    t=t_end-t_start
    logging.info('the programe time is %s') %t            
    logging.info('Done\n')

main()