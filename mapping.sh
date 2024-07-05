# mapping script includes loops for mapping mNET-seq and ChIP seq samples

# user inputs
INPUT_DIR=""          # path to folder containing .fastq files, change accordingly depending on the data
OUT_DIR=""            # path to save output files
STAR_INDEX=""         # path to STAR index
BOWTIE2_INDEX=""    # path to bowtie2 index

# variables created and used
STAR_OUT_DIR=$OUT_DIR/$PROJECT_NAME/$SEQ/processed_data/star_out_files
SORTED_BAM_DIR=$OUT_DIR/$PROJECT_NAME/$SEQ/processed_data/sorted_bam_files

mkdir -p $OUT_DIR
STAR_OUT_DIR=$OUT_DIR/"star_out_dir" #rename accordingly
BOWTIE2_OUT_DIR=$OUT_DIR/"bowtie2_out_dir"

declare -a INPUT_FASTQ_FILES
for f in $INPUT_DIR/*.fastq*; do
    INPUT_FASTQ_FILES[${#INPUT_FASTQ_FILES[@]}]=$(echo "$f");
done

N_FASTQ=${#INPUT_FASTQ_FILES[@]}

for i in $(seq 1 $N_FASTQ); do
    echo ${INPUT_FASTQ_FILES[$i-1]}
done

SAMPLES=()
for i in $(seq 1 $(($N_FASTQ / 2)));do
    temp=${INPUT_FASTQ_FILES[2*$i-1]##*/}
    temp=${temp::-9}
    SAMPLES+=("$temp")
done

for i in $(seq 1 $(($N_FASTQ / 2))); do
    echo ${SAMPLES[$i-1]}
done

# mapping with STAR
mkdir -p $STAR_OUT_DIR

for i in $(seq 1 $(($N_FASTQ / 2))); do
    echo $STAR_OUT_DIR/${SAMPLES[$i-1]}
    STAR --genomeDir $STAR_INDEX --readFilesIn \
    $INPUT_DIR/${SAMPLES[$i-1]}_1.fastq $INPUT_DIR/${SAMPLES[$i-1]}_2.fastq --runThreadN 8 \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $STAR_OUT_DIR/${SAMPLES[$i-1]} --outFilterMultimapNmax 1 \
    --outFilterMismatchNoverLmax 0.02 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 --alignIntronMax 500000
done
wait

cp $STAR_OUT_DIR/*.bam $SORTED_BAM_DIR
# rename to remove the Aligned.sortedByCoord.out
for file in $STAR_OUT_DIR/*.bam; do
    mv "$file" "${file/Aligned.sortedByCoord.out/}"
done

for file in $STAR_OUT_DIR/*.bam; do
    samtools index "$file" &
done
wait



# mapping with Bowtie2
mkdir -p $BOWTIE2_OUT_DIR

for i in $(seq 1 $(($N_FASTQ / 2))); do
    echo ${FastqFileList[2*$i-1]}
    echo ${FastqFileList[2*$i-2]}
    # Paired end
    bowtie2 -x $BOWTIE2_INDEX -1 $INPUT_DIR/${SAMPLES[$i-1]}_1.fastq -2 $INPUT_DIR/${SAMPLES[$i-1]}_2.fastq -S $BOWTIE2_OUT_DIR/${SAMPLES[$i-1]}.sam &
done
wait

for i in $(seq 1 $(($N_FASTQ / 2))); do
    samtools view -bS $BOWTIE2_OUT_DIR/${SAMPLES[$i-1]}.sam > $BOWTIE2_OUT_DIR/${SAMPLES[$i-1]}.bam &
done
wait


for file in $BOWTIE2_OUT_DIR/*.bam; do
    samtools index "$file" &
done
wait
