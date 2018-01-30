#!/bin/bash
#######################################################################
# ======= PBS OPTIONS ======= 
### Specify queue to run - mamba, Cobra, or copperhead
#PBS -q copperhead
### Set the job name
#PBS -N Metaphlan2
### Specify the # of cpus for your job.
#PBS -l nodes=1:ppn=16
### Adjust walltime below (default walltime = 7 days, or 168 hours)
### if you require > 7 days, INCREASE to estimated # hours needed
### if you DON'T require 7 days DECREASE to estimated # hours needed
### (hint: jobs with smaller walltime value tend to run sooner)
#PBS -l walltime=24:00:00
### pass the full environment(its like adding the path, so its useful to pass the enviornment)
#PBS -V
# send PBS output to /dev/null  (we redirect it below)
#PBS -o /dev/null
#PBS -e /dev/null
# ===== END PBS OPTIONS =====

#######################################################################
SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`
exec 1>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out 2>&1


cd /scratch/alulla/tp1/

module load python
module load bowtie2/2.1.0
export PATH=$PATH:~/.local/bin/:/users/alulla/metaphlan2

for file1 in ./*R1*;
    do file2=${file1/R1/R2};
    cat ${file1} ${file2} > ${file1%R1*}merged.fastq.gz;
    done

for f in *merged.fastq.gz;
    do metaphlan2.py $f --input_type multifastq --nproc 16 > ${f%merged.fastq.gz}_profile.txt;
        done

users/alulla/metaphlan2/utils/merge_metaphlan_tables.py *profile.txt > merged_abundance_table.txt

grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
grep -E "(g__)|(^ID)" merged_abundance_table.txt | grep -v "s__" | sed 's/^.*g__//g' > merged_abundance_table_genus.txt
grep -E "(f__)|(^ID)" merged_abundance_table.txt | grep -v "g__" | sed 's/^.*f__//g' > merged_abundance_table_family.txt
grep -E "(o__)|(^ID)" merged_abundance_table.txt | grep -v "f__" | sed 's/^.*o__//g' > merged_abundance_table_order.txt
grep -E "(c__)|(^ID)" merged_abundance_table.txt | grep -v "o__" | sed 's/^.*c__//g' > merged_abundance_table_class.txt
grep -E "(p__)|(^ID)" merged_abundance_table.txt | grep -v "c__" | sed 's/^.*p__//g' > merged_abundance_table_phylum.txt

#/users/alulla/hclust2/hclust2.py -i merged_abundance_table_species.txt -o abundance_heatmap_species_top100nomin.png --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 3 --slabel_size 3 --max_flabel_len 100 --max_slabel_len 100 --dpi 1000 --ftop 100







#Run code like on the command line. Ex:
#bash script.sh
#python script.py
#perl script.pl
#module load R
#R --no-save < script.R
#module load matlab
#matlab -nodisplay <script.m
