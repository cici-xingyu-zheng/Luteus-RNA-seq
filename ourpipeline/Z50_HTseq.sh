## Count reads from a collection of BAM files -> single gff file.
##
## Ron Smith - 2/13/2017
## rdsmith@email.wm.edu
##
## --------------------------------------------------------------

## ----- Settings ---------------------------------
## Path to main project directory
MAIN=/Users/puzey/Desktop/Cici

## Path to gff file
GFF=$MAIN/Lut_Genome/MimLut_lcl.gff

## Path to the BAM files
BAM=$MAIN/SAMBAM

## Path to output directory
OUT=$MAIN/Coverage

## Name of a log file
LOG=$MAIN/Log_Z50_HTSeq.txt
## ----- End of settings --------------------------


## Ensure output directory exists
mkdir $OUT

## Create log file (overwrite existing)
echo 'Starting : ' $(date) > $LOG

## Loop through all the BAM files
for file in $BAM/*.bam; do
    ## Make new file name
    base=${file##*/}
    outfile=${base/.bam/.counts.txt}
    outfile=$OUT/$outfile

    ## Print file names to screen
    echo -e '\n'
    echo -e '\tFile:\t\t' $file
    echo -e '\tOutput:\t\t' $outfile
    
    ## Run HTseq
    [ -e "$file" ] && htseq-count -f bam -i Parent -t exon -s no "$file" $GFF > "$outfile" || echo 'Problem with file ' $file >> $LOG
done

echo 'Process complete : ' $(date) >> $LOG
