import sys,os,re,argparse,logging,gzip
import multiprocessing
from Bio import SeqIO
import gtdb_result_transfer

    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fa_dir','-f',help='Directory to fasta files')
    parser.add_argument('--suffix','-x',default='fna',help='Suffix of the fasta files')
    #parser.add_argument('--taxid',help='Assembly ID to Taxid file')
    parser.add_argument('--outdir','-o',default='./kraken_library',help='Output kraken Directory',required=True)
    parser.add_argument('--gtdbtax',help='GTDB Taxonomy File')
    parser.add_argument('--gtdbres',help='GTDBtk Result. All summary.tsv has to be combined')
    args=parser.parse_args()
    try:
        if args.gtdbres is None and args.gtdbtax is None:
            parser.print_help()
            sys.exit(1)
        elif args.gtdbres and args.gtdbtax:
            parser.print_help()
            sys.exit(1)

    logging.basicConfig(filename=os.path.join(outdir,'genome2kraken.log'), level=logging.INFO)

    #Parse Assembly ID to Taxid List
    logging.info('Starting to Parse Taxonomy File\n')
    if args.gtdbres:
        ass2id=parse_gtdb_result(args.taxid)
        mode='gtdbres'
    elif args.gtdbtax:
        ass2id=parsingGTDBtax(args.outdir,args.gtdbtax)
        mode='gtdbtax'
    logging.info('Finished\n')
    
    if args.fna_dir: 
    #Add taxid to kraken lib fna
        t_start=time.time()
        pool=multiprocessing.Pool(10)
        for fna_file in os.listdir(args.fna_dir):
            fna_abs=os.path.abspath(os.path.join(args.fna_dir,fna_file))
            if os.path.splitext(fna_abs)[1]!='.'+args.suffix: continue #Checking file suffix
            pool.apply_async(addKraken,(args.outdir,fna_abs,ass2id))
    
        pool.close()
        pool.join()
        t_end=time.time()
        t=t_end-t_start
        logging.info('the programe time is %f' % t) 
    logging.info('Done\n')

main()
