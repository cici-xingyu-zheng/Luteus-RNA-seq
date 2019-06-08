## Script to produce nucleotide distribution and quality plots from a collection
## of fastq files.  Fastq files should be in a 'fastq' sub-folder of the current
## directory.  Script expects files ending in '.fastq.gz'.  If the fastq files
## are not gzipped, then this code will generate an error
##
##
## Ron Smith - 2/13/2017
## rdsmith@email.wm.edu
##
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

## Produce the plots
#for i in ../Stats/*.stats; do
#    base=${i##*/}
#    newfile=${base%.stats}_QualBox.png
#    [ -e "$i" ] && /Users/puzey/homebrew-20170717/local/Cellar/fastx_toolkit/0.0.14/bin/fastq_quality_boxplot_graph.sh -i "$i" -o "../Qual/$newfile" || echo "Somthing went wrong while making the qual plots"
    
#    newfile=${base%.stats}_NucDist.png
#    [ -e "$i" ] && /Users/puzey/homebrew-20170717/local/Cellar/fastx_toolkit/0.0.14/bin/fastx_nucleotide_distribution_graph.sh -i "$i" -o "../NucDist/$newfile" || echo "Somthing went wrong while making the nuc dists"
#done
