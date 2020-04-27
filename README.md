# GTDB_Kraken
Kraken(2) database based on the GTDB project

This page will contain links, notes and scripts for kraken databases we produced based on the GTDB database (http://gtdb.ecogenomic.org/). The kraken database contains 1 genome per taxon recognized in the GTDB database. Given that the the taxonomy of the GTDB database is based on genomic data, we expect this database to result in better performance for taxonomy-based classifier algorithms.    

## Version 1
Notes:
- Based on GTDB Release 03-RS86 (19th August 2018)
- Not the complete GTDB database, taxa without a genus designation have been omitted (10,994 out of 11042 taxa in this version)
- Reference sequences have been selected ad hoc, without regard for the quality of the sequence. This will be dealt with in future versions  

### Kraken2 db (65 Gb, tarball 48 Gb):

`gtdbk2_bacterial_v1.tar.gz md5sum: eda7855cb38a14b4222381ac5b27fe4b`

https://drive.google.com/file/d/18E0W_ezNLAhxxZjjelQYwoLA_0oBm4Lo/view?usp=sharing

Use `sh scripts/google_download.sh` to download the database from your command line. 

## Building from source

To download all fasta files and compile from source, you can run the provided script with the input file from the GTDB project.  Next, use the `--add-to-library` and `--build` functions in Kraken to format the database.  Example commands are below.

### Requirements

* Edirect - https://www.ncbi.nlm.nih.gov/books/NBK179288/
* Perl
  * BioPerl - https://bioperl.org/INSTALL.html
* Kraken or Kraken2

### Example commands for building from source

    # Creates the directory library/gtdb and adds fasta files.
    # Also creates a taxonomy folder, compatible with Kraken.
    perl scripts/gtdbToTaxonomy.pl --infile data/gtdb.2018-12-10.tsv
    # Format the database.
    db=GTDB_Kraken
    for i in library/gtdb/*.fna; do 
      kraken-build --add-to-library $i --db $db
    done
    mv -v taxonomy $db
    kraken-build --build --db $db --threads 8

### Cleanup (optional)

You can optionally remove all intermediate fasta files and also run the clean utility in Kraken.

    rm -rvf library
    kraken-build --clean $db

## Known issues

* Error 429 too many requests. NCBI is receiving too many general requests, but you can carve out a special place for yourself by getting an NCBI API key.  Log into your NCBI profile and copy your key.  Then, add it to your environment like so: `export NCBI_API_KEY=1fe2...` 

## Authors

* Henk den Bakker @hcdenbakker
* Lee Katz @lskatz


## Modifier

* Jiahui Zhu @magcurly

## Following Rewrite information

### Low efficiency in the original script

I found it's quite slow using the original scripts since the low efficiency by processing genomes one by one. Besides, I found some bugs in the original scripts. I fixed some of the bugs in the perl script while I rewrote the perl scripts to a python one. This python script (Under Python3) is fit for those who already download the gtdb release 89 data pack (especially for those who use GTDBtk). I testified it with 24706 representative genomes selected by GTDB. You can find your taxfile in gtdb taxonomy/ dicrectory (unpack the gtdb tar.gz first) and your fna files in fastani/database/. Since I have to test it with the original perl script, I gunzipped all fna.gz.

You have to noticed that I use multiprocessing in the python script in order to improve the process speed.

### Requirement package for the python script

 - biopython
 - multiprocessing
 - gzip
