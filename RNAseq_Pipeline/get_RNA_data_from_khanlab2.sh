#!/usr/bin/bash
for i in `cat rna_stack`;
do mkdir -p /data/khanlab/projects/ChIP_seq/RNA_DATA/$i/pbs_log
ln -s /data/khanlab2/DATA/$i/*.fastq.gz /data/khanlab/projects/ChIP_seq/RNA_DATA/$i
chgrp -R khanlab /data/khanlab/projects/ChIP_seq/RNA_DATA/$i
chmod -R g+w /data/khanlab/projects/ChIP_seq/RNA_DATA/$i
echo "linked fastq from khanlab2 for $i"
done
