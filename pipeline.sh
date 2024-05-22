#!/bin/bash

for filename in /*_1.fastq.gz
do
	base=$(basename $filename _1.fastq.gz)
	echo $base
	java -jar /home/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 90 -phred33 ${base}_1.fastq.gz ${base}_2.fastq.gz ${base}.pe1.fq.gz ${base}.se1.fq.gz ${base}.pe2.fq.gz ${base}.se2.fq.gz ILLUMINACLIP:/home/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:6:30 MINLEN:100

done

for i in /*.paired_1.fastq
do
	base=$(basename $i .paired_1.fastq)
	echo $base

	time metawrap assembly -1 ${base}.paired_1.fastq -2 ${base}.paired_2.fastq -o ./${base}_assembly -t 80

done


for i in /*_assembly
do
	base=$(basename $i _assembly)
	echo $base

	metawrap binning -a ${base}_assembly/final_assembly.fasta -o ${base}_bin -t 80 -m 220 --metabat2 --maxbin2 --concoct /${base}.paired_1.fastq /${base}.paired_2.fastq

done


metawrap bin_refinement -o /${base}_binrefinement -t 60 -m 350 -A /metabat2_bins -B /maxbin2_bins -C /concoct_bins -c 50 -x 10

dRep dereplicate -p 60 dRep_out -g all_bin/*.fa -comp 50 -con 10 -sa 0.95

gtdbtk classify_wf --cpus 100 --genome_dir /metawrap_50_10_bins/ --out_dir ./${base}_results -x fa

VIBRANT_run.py -i /metawrap_50_10_bins/${base}.fa -f nucl -t 30 -l 5000 -no_plot -folder ${base}_VB_out -d /VIBRANT_db/databases

virsorter run -w ${base}_VS_out -i /metawrap_50_10_bins/${base}.fa --min-length 5000 --rm-tmpdir --verbose -j 100 all

cd-hit-est -c 0.95 -aS 0.8 -T 0 -M 0 -d 200 -i virus.fasta -o virus.clstr.fasta

checkv end_to_end virus.clstr.fasta checkv-output -t 60

prodigal -a pro.faa -d nucl.fna -f gff -i nuc.fasta -p meta -q -m

vcontact2_gene2genome -p clstr.faa -o clstr_mapping_output.csv -s Prodigal-FAA

vcontact2 --rel-mode 'Diamond' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /cluster_one-1.0.jar --db 'ProkaryoticViralRefSeq201-Merged' --verbose --threads 90 --raw-proteins clstr.faa --proteins-fp clstr_mapping_output.csv --output-dir vcontact2_output

source activate coverm

coverm genome --coupled /ERR712369.pe1.fq.gz /ERR712369.pe2.fq.gz /ERR712370.pe1.fq.gz /ERR712370.pe2.fq.gz /ERR712371-73.pe1.fq.gz /ERR712371-73.pe2.fq.gz /ERR712374-76.pe1.fq.gz /ERR712374-76.pe2.fq.gz /ERR712377-79.pe1.fq.gz /ERR712377-79.pe2.fq.gz /ERR712380.pe1.fq.gz /ERR712380.pe2.fq.gz /ERR712389.pe1.fq.gz /ERR712389.pe2.fq.gz /ERR712390.pe1.fq.gz /ERR712390.pe2.fq.gz /ERR712391.pe1.fq.gz /ERR712391.pe2.fq.gz /SRR11673965.pe1.fq.gz /SRR11673965.pe2.fq.gz /SRR11673976.pe1.fq.gz /SRR11673976.pe2.fq.gz /SRR11673987.pe1.fq.gz /SRR11673987.pe2.fq.gz /SRR11673989.pe1.fq.gz /SRR11673989.pe2.fq.gz /SRR11673990.pe1.fq.gz /SRR11673990.pe2.fq.gz /SRR11673991.pe1.fq.gz /SRR11673991.pe2.fq.gz /SRR11673992.pe1.fq.gz /SRR11673992.pe2.fq.gz /SRR11673993.pe1.fq.gz /SRR11673993.pe2.fq.gz /SRR11673994.pe1.fq.gz /SRR11673994.pe2.fq.gz /SRR11673995.pe1.fq.gz /SRR11673995.pe2.fq.gz /SRR11673996.pe1.fq.gz /SRR11673996.pe2.fq.gz /SRR11673997.pe1.fq.gz /SRR11673997.pe2.fq.gz /SRR11673998.pe1.fq.gz /SRR11673998.pe2.fq.gz /SRR11673999.pe1.fq.gz /SRR11673999.pe2.fq.gz /SRR11674000.pe1.fq.gz /SRR11674000.pe2.fq.gz /SRR11674001.pe1.fq.gz /SRR11674001.pe2.fq.gz /SRR11674002.pe1.fq.gz /SRR11674002.pe2.fq.gz /SRR11674003.pe1.fq.gz /SRR11674003.pe2.fq.gz /SRR11674004.pe1.fq.gz /SRR11674004.pe2.fq.gz /SRR11674005.pe1.fq.gz /SRR11674005.pe2.fq.gz /SRR11674007.pe1.fq.gz /SRR11674007.pe2.fq.gz /SRR11674008.pe1.fq.gz /SRR11674008.pe2.fq.gz /SRR11674009.pe1.fq.gz /SRR11674009.pe2.fq.gz /SRR11674011.pe1.fq.gz /SRR11674011.pe2.fq.gz /SRR11674012.pe1.fq.gz /SRR11674012.pe2.fq.gz /SRR11674013.pe1.fq.gz /SRR11674013.pe2.fq.gz /SRR11674014.pe1.fq.gz /SRR11674014.pe2.fq.gz /SRR11674015.pe1.fq.gz /SRR11674015.pe2.fq.gz/SRR11674016.pe1.fq.gz /SRR11674016.pe2.fq.gz /SRR11674017.pe1.fq.gz /SRR11674017.pe2.fq.gz /SRR11674018.pe1.fq.gz /SRR11674018.pe2.fq.gz /SRR11674019.pe1.fq.gz /SRR11674019.pe2.fq.gz /SRR11674020.pe1.fq.gz /SRR11674020.pe2.fq.gz /SRR11674021.pe1.fq.gz /SRR11674021.pe2.fq.gz /SRR11674022.pe1.fq.gz /SRR11674022.pe2.fq.gz /SRR11674023.pe1.fq.gz /SRR11674023.pe2.fq.gz /SRR11674024.pe1.fq.gz /SRR11674024.pe2.fq.gz /SRR11674025.pe1.fq.gz /SRR11674025.pe2.fq.gz /SRR11674026.pe1.fq.gz /SRR11674026.pe2.fq.gz /SRR11674027.pe1.fq.gz /SRR11674027.pe2.fq.gz /SRR11674028.pe1.fq.gz /SRR11674028.pe2.fq.gz /SRR11674029.pe1.fq.gz /SRR11674029.pe2.fq.gz /SRR11674030.pe1.fq.gz /SRR11674030.pe2.fq.gz /SRR11674031.pe1.fq.gz /SRR11674031.pe2.fq.gz /SRR11674032.pe1.fq.gz /SRR11674032.pe2.fq.gz /SRR11674034.pe1.fq.gz /SRR11674034.pe2.fq.gz /SRR11674035.pe1.fq.gz /SRR11674035.pe2.fq.gz /SRR11674036.pe1.fq.gz /SRR11674036.pe2.fq.gz /SRR11674037.pe1.fq.gz /SRR11674037.pe2.fq.gz /SRR11674038.pe1.fq.gz /SRR11674038.pe2.fq.gz /SRR11674039.pe1.fq.gz /SRR11674039.pe2.fq.gz /SRR11674040.pe1.fq.gz /SRR11674040.pe2.fq.gz /SRR11674041.pe1.fq.gz /SRR11674041.pe2.fq.gz /SRR11674042.pe1.fq.gz /SRR11674042.pe2.fq.gz /SRR11674043.pe1.fq.gz /SRR11674043.pe2.fq.gz /SRR11674044.pe1.fq.gz /SRR11674044.pe2.fq.gz /SRR11674045.pe1.fq.gz /SRR11674045.pe2.fq.gz /SRR11674046.pe1.fq.gz/SRR11674046.pe2.fq.gz /SRR11674047.pe1.fq.gz /SRR11674047.pe2.fq.gz /SRR11674048.pe1.fq.gz /SRR11674048.pe2.fq.gz /SRR11674049.pe1.fq.gz /SRR11674049.pe2.fq.gz /SRR11674050.pe1.fq.gz /SRR11674050.pe2.fq.gz /SRR11674051.pe1.fq.gz /SRR11674051.pe2.fq.gz /SRR11674052.pe1.fq.gz /SRR11674052.pe2.fq.gz /SRR11674053.pe1.fq.gz /SRR11674053.pe2.fq.gz /SRR11674054.pe1.fq.gz /SRR11674054.pe2.fq.gz /SRR2107210-11.pe1.fq.gz /SRR2107210-11.pe2.fq.gz /SRR2107212-13.pe1.fq.gz /SRR2107212-13.pe2.fq.gz /SRR2107215-16.pe1.fq.gz /SRR2107215-16.pe2.fq.gz /SRR2107218-19.pe1.fq.gz /SRR2107218-19.pe2.fq.gz /SRR3501849-51.pe1.fq.gz /SRR3501849-51.pe2.fq.gz /SRR3501850-52.pe1.fq.gz /SRR3501850-52.pe2.fq.gz /SRR3501853-62.pe1.fq.gz /SRR3501853-62.pe2.fq.gz /SRR3501854-73.pe1.fq.gz /SRR3501854-73.pe2.fq.gz /SRR3501855-84.pe1.fq.gz /SRR3501855-84.pe2.fq.gz /SRR3501856-85.pe1.fq.gz /SRR3501856-85.pe2.fq.gz /SRR3501857-86.pe1.fq.gz /SRR3501857-86.pe2.fq.gz /SRR3501858-87.pe1.fq.gz /SRR3501858-87.pe2.fq.gz /SRR3501859-88.pe1.fq.gz /SRR3501859-88.pe2.fq.gz /SRR3501861-89.pe1.fq.gz /SRR3501861-89.pe2.fq.gz /SRR5570991-92-1003.pe1.fq.gz /SRR5570991-92-1003.pe2.fq.gz /SRR5571008-09-11.pe1.fq.gz /SRR5571008-09-11.pe2.fq.gz /SRR6032601.pe1.fq.gz /SRR6032601.pe2.fq.gz /SRR6032602.pe1.fq.gz /SRR6032602.pe2.fq.gz /SRR6032603.pe1.fq.gz /SRR6032603.pe2.fq.gz /SRR9260993.pe1.fq.gz /SRR9260993.pe2.fq.gz /SRR9260994.pe1.fq.gz /SRR9260994.pe2.fq.gz /SRR9260995.pe1.fq.gz /SRR9260995.pe2.fq.gz /SRR9260996.pe1.fq.gz/SRR9260996.pe2.fq.gz /SRR9260997.pe1.fq.gz /SRR9260997.pe2.fq.gz /SRR9260998.pe1.fq.gz /SRR9260998.pe2.fq.gz /SRR9260999.pe1.fq.gz /SRR9260999.pe2.fq.gz /SRR9261000.pe1.fq.gz /SRR9261000.pe2.fq.gz /SRR9261001.pe1.fq.gz /SRR9261001.pe2.fq.gz /SRR9261002.pe1.fq.gz /SRR9261002.pe2.fq.gz /SRR9261003.pe1.fq.gz /SRR9261003.pe2.fq.gz /SRR9261004.pe1.fq.gz /SRR9261004.pe2.fq.gz /SRR9261005.pe1.fq.gz /SRR9261005.pe2.fq.gz /SRR9261006.pe1.fq.gz /SRR9261006.pe2.fq.gz /SRR9261007.pe1.fq.gz /SRR9261007.pe2.fq.gz /SRR9261008.pe1.fq.gz /SRR9261008.pe2.fq.gz /SRR9261009.pe1.fq.gz /SRR9261009.pe2.fq.gz /SRR9261010.pe1.fq.gz /SRR9261010.pe2.fq.gz /SRR9261011.pe1.fq.gz /SRR9261011.pe2.fq.gz /SRR9261012.pe1.fq.gz /SRR9261012.pe2.fq.gz /SRR9261013.pe1.fq.gz /SRR9261013.pe2.fq.gz /SRR9261014.pe1.fq.gz /SRR9261014.pe2.fq.gz /SRR9261015.pe1.fq.gz /SRR9261015.pe2.fq.gz /SRR9261016.pe1.fq.gz /SRR9261016.pe2.fq.gz/SRR9261017.pe1.fq.gz /SRR9261017.pe2.fq.gz /SRR9261018.pe1.fq.gz /SRR9261018.pe2.fq.gz /SRR9261019.pe1.fq.gz /SRR9261019.pe2.fq.gz /SRR9261020.pe1.fq.gz /SRR9261020.pe2.fq.gz /SRR9261021.pe1.fq.gz /SRR9261021.pe2.fq.gz /SRR9261022.pe1.fq.gz /SRR9261022.pe2.fq.gz /SRR9261023.pe1.fq.gz /SRR9261023.pe2.fq.gz /SRR9261024.pe1.fq.gz /SRR9261024.pe2.fq.gz /SRR9261025.pe1.fq.gz /SRR9261025.pe2.fq.gz /SRR9261026.pe1.fq.gz /SRR9261026.pe2.fq.gz /SRR9261027.pe1.fq.gz /SRR9261027.pe2.fq.gz /SRR9261028.pe1.fq.gz /SRR9261028.pe2.fq.gz /SRR9261029.pe1.fq.gz /SRR9261029.pe2.fq.gz /SRR9261030.pe1.fq.gz /SRR9261030.pe2.fq.gz /SRR9261031.pe1.fq.gz /SRR9261031.pe2.fq.gz /SRR9261032.pe1.fq.gz /SRR9261032.pe2.fq.gz /SRR9261033.pe1.fq.gz /SRR9261033.pe2.fq.gz /SRR9261034.pe1.fq.gz /SRR9261034.pe2.fq.gz /SRR9261035.pe1.fq.gz /SRR9261035.pe2.fq.gz /SRR9261036.pe1.fq.gz /SRR9261036.pe2.fq.gz /SRR9261037.pe1.fq.gz /SRR9261037.pe2.fq.gz /SRR9261038.pe1.fq.gz /SRR9261038.pe2.fq.gz /SRR9261039.pe1.fq.gz /SRR9261039.pe2.fq.gz /SRR9261040.pe1.fq.gz /SRR9261040.pe2.fq.gz /SRR9261041.pe1.fq.gz /SRR9261041.pe2.fq.gz /SRR9261042.pe1.fq.gz /SRR9261042.pe2.fq.gz /SRR9261043.pe1.fq.gz /SRR9261043.pe2.fq.gz /SRR9261044.pe1.fq.gz /SRR9261044.pe2.fq.gz /SRR9261045.pe1.fq.gz /SRR9261045.pe2.fq.gz /SRR9261046.pe1.fq.gz /SRR9261046.pe2.fq.gz /SRR9261047.pe1.fq.gz /SRR9261047.pe2.fq.gz /SRR9261048.pe1.fq.gz /SRR9261048.pe2.fq.gz /SRR9261049.pe1.fq.gz /SRR9261049.pe2.fq.gz /SRR9261050.pe1.fq.gz/SRR9261050.pe2.fq.gz /SRR9261051.pe1.fq.gz /SRR9261051.pe2.fq.gz /SRR9261052.pe1.fq.gz /SRR9261052.pe2.fq.gz /SRR9852164.pe1.fq.gz /SRR9852164.pe2.fq.gz /SRR9852165.pe1.fq.gz /SRR9852165.pe2.fq.gz /SRR9852166.pe1.fq.gz /SRR9852166.pe2.fq.gz /SRR9852167.pe1.fq.gz /SRR9852167.pe2.fq.gz /SRR9852168.pe1.fq.gz /SRR9852168.pe2.fq.gz /SRR9852169.pe1.fq.gz /SRR9852169.pe2.fq.gz /SRR9852170.pe1.fq.gz /SRR9852170.pe2.fq.gz /SRR9852171.pe1.fq.gz /SRR9852171.pe2.fq.gz /SRR9852172.pe1.fq.gz /SRR9852172.pe2.fq.gz /SRR9852173.pe1.fq.gz /SRR9852173.pe2.fq.gz /SRR9852174.pe1.fq.gz /SRR9852174.pe2.fq.gz /SRR9852175.pe1.fq.gz /SRR9852175.pe2.fq.gz /SRR9852176.pe1.fq.gz /SRR9852176.pe2.fq.gz /SRR9852177.pe1.fq.gz /SRR9852177.pe2.fq.gz /SRR9852178.pe1.fq.gz /SRR9852178.pe2.fq.gz /SRR9852179.pe1.fq.gz /SRR9852179.pe2.fq.gz /SRR9852180.pe1.fq.gz /SRR9852180.pe2.fq.gz /SRR9852181.pe1.fq.gz /SRR9852181.pe2.fq.gz /SRR9852182.pe1.fq.gz /SRR9852182.pe2.fq.gz /SRR9852183.pe1.fq.gz /SRR9852183.pe2.fq.gz /SRR9852184.pe1.fq.gz /SRR9852184.pe2.fq.gz /SRR9852185.pe1.fq.gz /SRR9852185.pe2.fq.gz /SRR9852186.pe1.fq.gz /SRR9852186.pe2.fq.gz /SRR9852187.pe1.fq.gz /SRR9852187.pe2.fq.gz /SRR9852188.pe1.fq.gz /SRR9852188.pe2.fq.gz /SRR9852189.pe1.fq.gz /SRR9852189.pe2.fq.gz /SRR9852190.pe1.fq.gz /SRR9852190.pe2.fq.gz /SRR9852191.pe1.fq.gz /SRR9852191.pe2.fq.gz /SRR9852192.pe1.fq.gz /SRR9852192.pe2.fq.gz /SRR9852193.pe1.fq.gz /SRR9852193.pe2.fq.gz /SRR9852194.pe1.fq.gz /SRR9852194.pe2.fq.gz/SRR9852195.pe1.fq.gz /SRR9852195.pe2.fq.gz /SRR9852196.pe1.fq.gz /SRR9852196.pe2.fq.gz /SRR9852197.pe1.fq.gz /SRR9852197.pe2.fq.gz /SRR9852198.pe1.fq.gz /SRR9852198.pe2.fq.gz /SRR9852199.pe1.fq.gz /SRR9852199.pe2.fq.gz /SRR9852200.pe1.fq.gz /SRR9852200.pe2.fq.gz /SRR9852201.pe1.fq.gz /SRR9852201.pe2.fq.gz /SRR9852202.pe1.fq.gz /SRR9852202.pe2.fq.gz /SRR9852203.pe1.fq.gz /SRR9852203.pe2.fq.gz /SRR9852204.pe1.fq.gz /SRR9852204.pe2.fq.gz /SRR9852205.pe1.fq.gz /SRR9852205.pe2.fq.gz /SRR9852206.pe1.fq.gz /SRR9852206.pe2.fq.gz /SRR9852207.pe1.fq.gz /SRR9852207.pe2.fq.gz /SRR9852208.pe1.fq.gz /SRR9852208.pe2.fq.gz /SRR9852209.pe1.fq.gz /SRR9852209.pe2.fq.gz /SRR9852210.pe1.fq.gz /SRR9852210.pe2.fq.gz /SRR9852211.pe1.fq.gz /SRR9852211.pe2.fq.gz /SRR9852212.pe1.fq.gz /SRR9852212.pe2.fq.gz /SRR9852213.pe1.fq.gz /SRR9852213.pe2.fq.gz /SRR9852214.pe1.fq.gz /SRR9852214.pe2.fq.gz /SRR9852215.pe1.fq.gz /SRR9852215.pe2.fq.gz /SRR9852216.pe1.fq.gz /SRR9852216.pe2.fq.gz /SRR9852217.pe1.fq.gz /SRR9852217.pe2.fq.gz /SRR9852218.pe1.fq.gz /SRR9852218.pe2.fq.gz /SRR9852219.pe1.fq.gz /SRR9852219.pe2.fq.gz /SRR9852220.pe1.fq.gz /SRR9852220.pe2.fq.gz -d /data/LJ/drep/input2 -x fa -t 60 --min-read-percent-identity 95 --min-read-aligned-percent 75 -m covered_fraction -o all.MAG.abundance.csv


exec_annotation -o ko.detail.tsv contigs.faa --cpu 80 -E 1e-5 -f detail-tsv -p /profiles -k /ko_list

hmmsearch -o pfam.out --tblout pfam.tblout -E 1e-5 --incE 1e-5 --cpu 80 /pfam/Pfam-A.hmm contigs.faa

Propagate -f contig.fa -r /${base}_R1.fastq.gz /${base}_R2.fastq.gz -o ${base}_output -v /prophage_coordinates.txt -t 90 --clean

blastn -query contig.fa -db database.db -out blastn.output.csv -evalue 1e-5 -outfmt "6 qseqid sseqid pident qcovs qcovhsp evalue bitscore length mismatch gapopen qstart qend sstart send" -perc_identity 95 -qcov_hsp_perc 90 -max_target_seqs 1 -num_threads 90
