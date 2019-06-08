## Script to sort SAM files and output as BAM. Original SAM files are backed up
## to a different folder.
##
## Ron Smith - 2/13/2017
## rdsmith@email.wm.edu
##
## ----------------------------------------------------------------------------

## ----- Settings --------------------------------------------
## Path to main project directory
MAIN=/Users/puzey/Desktop/Cici

## Path to the SAM files
SAM=$MAIN/SAMBAM

## Path to backup the original files
BAK=$MAIN/BAK

## Name of a log file
LOG=$MAIN/Log_Z30_Sam2Bam.txt
## ----- End of settings -------------------------------------


## Ensure backup directory exists
mkdir $BAK

## Start a new log file (overwrite existing)
echo 'Starting : ' $(date) > $LOG

## Loop through the SAM files
for file in $SAM/*.sam; do
    
    ## Make the output file name.  Print file names to screen
    outfile=${file/.sam/.sorted.bam}
    echo -e '\n'
    echo -e '\tFile:\t\t' $file
    echo -e '\tOutput:\t\t' $outfile
    
    ## Run samtools
    [ -e "$file" ] && samtools sort -O bam -T "$file" -o "$outfile" "$file" || echo "Something went wrong with file " $file >> $LOG
    echo $file ' : ' $(date) >> $LOG

    ## Move the original file to the backup folder
    base=${file##*/} 
    mv $file $BAK/$base
done

echo 'Done : ' $(date) >> $LOG
