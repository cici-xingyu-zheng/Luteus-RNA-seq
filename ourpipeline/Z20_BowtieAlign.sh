## Align reads with bowtie2. Multiple fastq files -> Single fasta file.  
## Fasta file must be indexed ahead of time.  Fastq files may be gzipped or not,
## but should contain 'fastq' (e.g. not 'fq') in the file name.  No other files 
## in the fastq folder should have 'fastq' in their file name.
##
##
## ----------------------------------------------------------------------


## ----- Settings ----------------------------
## Number of cores bowtie will use.  Careful!
NUM_CORES=10

## Path to main directory
MAIN=/Users/puzey/Desktop/Cici

## Path to fasta file (without file extension)
REF=$MAIN/Lut_Genome/Mlut

## Path to fastq files
FASTQ=$MAIN/fastq

## Output directory for SAM files
OUT=$MAIN/SAMBAM

## Name of log file
LOG=$MAIN/Log_Z20.txt
## ----- End of settings ---------------------


## Ensure output directory exists
mkdir $OUT
## Start new log file (overwrite existing file)
echo 'Starting : ' $(date) > $LOG

## Loop throgh fastq files
for file in $FASTQ/*fastq*; do
    
    ## Create new file name.  Print new and old names to screen.
    echo -e '\n'
    echo -e '\tCurrent file:\t\t' $file
    base=${file##*/}
    outfile=${base/.fastq.gz/.sam}
    outfile=${outfile/.fastq/.sam}
    outfile=$OUT/$outfile
    echo -e '\tOutput will be:\t\t' $outfile
    
    ## Run bowtie
    [ -e "$file" ] && bowtie2 -p $NUM_CORES --local --very-sensitive-local -x $REF -U "$file" -S "$outfile" || echo "Something went wrong with file " $file >> $LOG
    echo $file ' : ' $(date) >> $LOG
done

echo 'Process complete : ' $(date) >> $LOG
