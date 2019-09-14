## Script to produce nucleotide distribution and quality plots from a collection
## of fastq files.  Fastq files should be in a 'fastq' sub-folder of the current
## directory.  Script expects files ending in '.fastq.gz'.  If the fastq files
## are not gzipped, then this code will generate an error
##
## -----------------------------------------------------------------------------


## Make directory to save stats files
mkdir ../Stats

## Produce the stats files
for i in ../fastq/*.fastq.gz; do
    base=${i##*/}    
    newfile=${base%.fastq.gz}.stats
    [ -e "$i" ] && gunzip -c "$i" | /Users/puzey/homebrew-20170717/local/Cellar/fastx_toolkit/0.0.14/bin/fastx_quality_stats -o "../Stats/$newfile" || echo "Somthing went wrong while making the stats files"
done

## Make directories to save plots
mkdir ../Qual
mkdir ../NucDist

