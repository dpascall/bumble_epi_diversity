#!/bin/bash

###convert and trim sequences

~/Desktop/seqtk/seqtk seq -Q64 -V  B_lucorum_l1_1.fq > B_lucorum_l1_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_lucorum_l1_2.fq > B_lucorum_l1_2.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_lucorum_run2_1.fq > B_lucorum_run2_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_lucorum_run2_2.fq > B_lucorum_run2_2.sanger.fq

rm B_lucorum_l1_1.fq B_lucorum_l1_2.fq B_lucorum_run2_1.fq B_lucorum_run2_2.fq

sickle pe -f B_lucorum_l1_1.sanger.fq -r B_lucorum_l1_2.sanger.fq -t sanger -o trimmed_BL_1_1.fastq -p trimmed_BL_1_2.fastq -s trimmed_BL_1_S.fastq
sickle pe -f B_lucorum_run2_1.sanger.fq -r B_lucorum_run2_2.sanger.fq -t sanger -o trimmed_BL_2_1.fastq -p trimmed_BL_2_2.fastq -s trimmed_BL_2_S.fastq

rm B_lucorum_l1_1.sanger.fq B_lucorum_l1_2.sanger.fq B_lucorum_run2_1.sanger.fq B_lucorum_run2_2.sanger.fq

~/Desktop/seqtk/seqtk seq -Q64 -V  B_terrestris_l1_1.fq > B_terrestris_l1_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_terrestris_l1_2.fq > B_terrestris_l1_2.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_terrestris_run2_1.fq > B_terrestris_run2_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_terrestris_run2_2.fq > B_terrestris_run2_2.sanger.fq

rm B_terrestris_l1_1.fq B_terrestris_l1_2.fq B_terrestris_run2_1.fq B_terrestris_run2_2.fq

sickle pe -f B_terrestris_l1_1.sanger.fq -r B_terrestris_l1_2.sanger.fq -t sanger -o trimmed_BT_1_1.fastq -p trimmed_BT_1_2.fastq -s trimmed_BT_1_S.fastq
sickle pe -f B_terrestris_run2_1.sanger.fq -r B_terrestris_run2_2.sanger.fq -t sanger -o trimmed_BT_2_1.fastq -p trimmed_BT_2_2.fastq -s trimmed_BT_2_S.fastq

rm B_terrestris_l1_1.sanger.fq B_terrestris_l1_2.sanger.fq B_terrestris_run2_1.sanger.fq B_terrestris_run2_2.sanger.fq

~/Desktop/seqtk/seqtk seq -Q64 -V  B_pascuorum_l1_1.fq > B_pascuorum_l1_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_pascuorum_l1_2.fq > B_pascuorum_l1_2.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_pascuorum_run2_1.fq > B_pascuorum_run2_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  B_pascuorum_run2_2.fq > B_pascuorum_run2_2.sanger.fq

rm B_pascuorum_l1_1.fq B_pascuorum_l1_2.fq B_pascuorum_run2_1.fq B_pascuorum_run2_2.fq

sickle pe -f B_pascuorum_l1_1.sanger.fq -r B_pascuorum_l1_2.sanger.fq -t sanger -o trimmed_BP_1_1.fastq -p trimmed_BP_1_2.fastq -s trimmed_BP_1_S.fastq
sickle pe -f B_pascuorum_run2_1.sanger.fq -r B_pascuorum_run2_2.sanger.fq -t sanger -o trimmed_BP_2_1.fastq -p trimmed_BP_2_2.fastq -s trimmed_BP_2_S.fastq

rm B_pascuorum_l1_1.sanger.fq B_pascuorum_l1_2.sanger.fq B_pascuorum_run2_1.sanger.fq B_pascuorum_run2_2.sanger.fq

~/Desktop/seqtk/seqtk seq -Q64 -V  Bombus_sp_l1_1.fq > Bombus_sp_l1_1.sanger.fq
~/Desktop/seqtk/seqtk seq -Q64 -V  Bombus_sp_l1_2.fq > Bombus_sp_l1_2.sanger.fq

rm Bombus_sp_l1_1.fq Bombus_sp_l1_2.fq

sickle pe -f Bombus_sp_l1_1.sanger.fq -r Bombus_sp_l1_2.sanger.fq -t sanger -o B_1.fastq -p B_2.fastq -s B_S.fastq

rm Bombus_sp_l1_1.sanger.fq Bombus_sp_l1_2.sanger.fq

###merge

cat trimmed_BL_1_1.fastq trimmed_BL_2_1.fastq > ./Merged/BL_1.fastq
cat trimmed_BL_1_2.fastq trimmed_BL_2_2.fastq > ./Merged/BL_2.fastq
cat trimmed_BL_1_S.fastq trimmed_BL_2_S.fastq > ./Merged/BL_S.fastq

rm trimmed_BL_1_1.fastq trimmed_BL_2_1.fastq trimmed_BL_1_2.fastq trimmed_BL_2_2.fastq trimmed_BL_1_S.fastq trimmed_BL_2_S.fastq

cd Merged

flash --threads=6 -o BL BL_1.fastq BL_2.fastq

cat BL_S.fastq BL.extendedFrags.fastq > BL_Smerged.fastq

rm BL_S.fastq BL_1.fastq BL_2.fastq BL.extendedFrags.fastq

mv BL.notCombined_1.fastq BL_1.fastq

mv BL.notCombined_2.fastq BL_2.fastq

cd ..

cat trimmed_BP_1_1.fastq trimmed_BP_2_1.fastq > ./Merged/BP_1.fastq
cat trimmed_BP_1_2.fastq trimmed_BP_2_2.fastq > ./Merged/BP_2.fastq
cat trimmed_BP_1_S.fastq trimmed_BP_2_S.fastq > ./Merged/BP_S.fastq

rm trimmed_BP_1_1.fastq trimmed_BP_2_1.fastq trimmed_BP_1_2.fastq trimmed_BP_2_2.fastq trimmed_BP_1_S.fastq trimmed_BP_2_S.fastq

cd Merged

flash --threads=6 -o BP BP_1.fastq BP_2.fastq

cat BP_S.fastq BP.extendedFrags.fastq > BP_Smerged.fastq

rm BP_S.fastq BP_1.fastq BP_2.fastq BP.extendedFrags.fastq

mv BP.notCombined_1.fastq BP_1.fastq

mv BP.notCombined_2.fastq BP_2.fastq

cd ..

cat trimmed_BT_1_1.fastq trimmed_BT_2_1.fastq > ./Merged/BT_1.fastq
cat trimmed_BT_1_2.fastq trimmed_BT_2_2.fastq > ./Merged/BT_2.fastq
cat trimmed_BT_1_S.fastq trimmed_BT_2_S.fastq > ./Merged/BT_S.fastq

rm trimmed_BT_1_1.fastq trimmed_BT_2_1.fastq trimmed_BT_1_2.fastq trimmed_BT_2_2.fastq trimmed_BT_1_S.fastq trimmed_BT_2_S.fastq

cd Merged

flash --threads=6 -o BT BT_1.fastq BT_2.fastq

cat BT_S.fastq BT.extendedFrags.fastq > BT_Smerged.fastq

rm BT_S.fastq BT_1.fastq BT_2.fastq BT.extendedFrags.fastq

mv BT.notCombined_1.fastq BT_1.fastq

mv BT.notCombined_2.fastq BT_2.fastq

flash --threads=6 -o B ../B_1.fastq ../B_2.fastq

cat ../B_S.fastq B.extendedFrags.fastq > B_Smerged.fastq

rm ../B_S.fastq ../B_1.fastq ../B_2.fastq B.extendedFrags.fastq

mv B.notCombined_1.fastq B_1.fastq

mv B.notCombined_2.fastq B_2.fastq

##align

/Applications/MOSAIK_2.1.73/MosaikBuild -q BL_1.fastq -q2 BL_2.fastq -st illumina -mfl 200 -out lucorum.mkb
rm BL_1.fastq BL_2.fastq
/Applications/MOSAIK_2.1.73/MosaikBuild -q BL_Smerged.fastq -st illumina -out lucorum_s.mkb
rm BL_Smerged.fastq

/Applications/MOSAIK_2.1.73/MosaikBuild -q BP_1.fastq -q2 BP_2.fastq -st illumina -mfl 200 -out pascuorum.mkb
rm BP_1.fastq BP_2.fastq
/Applications/MOSAIK_2.1.73/MosaikBuild -q BP_Smerged.fastq -st illumina -out pascuorum_s.mkb
rm BP_Smerged.fastq

/Applications/MOSAIK_2.1.73/MosaikBuild -q BT_1.fastq -q2 BT_2.fastq -st illumina -mfl 200 -out terrestris.mkb
rm BT_1.fastq BT_2.fastq
/Applications/MOSAIK_2.1.73/MosaikBuild -q BT_Smerged.fastq -st illumina -out terrestris_s.mkb
rm BT_Smerged.fastq

/Applications/MOSAIK_2.1.73/MosaikBuild -q B_1.fastq -q2 B_2.fastq -st illumina -mfl 200 -out bombus.mkb
rm B_1.fastq B_2.fastq
/Applications/MOSAIK_2.1.73/MosaikBuild -q B_Smerged.fastq -st illumina -out bombus_s.mkb
rm B_Smerged.fastq

/Applications/MOSAIK_2.1.73/MosaikBuild -fr ReftrimmedtoRdRpsdisamb.fasta -oa ref.dat

/Applications/MOSAIK_2.1.73/MosaikAligner -in bombus.mkb -out bombus.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7
/Applications/MOSAIK_2.1.73/MosaikAligner -in bombus_s.mkb -out bombus_s.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7

rm bombus.mkb bombus_s.mka

/Applications/MOSAIK_2.1.73/MosaikAligner -in lucorum.mkb -out lucorum.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7
/Applications/MOSAIK_2.1.73/MosaikAligner -in lucorum_s.mkb -out lucorum_s.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7

rm lucorum.mkb lucorum_s.mkb

/Applications/MOSAIK_2.1.73/MosaikAligner -in pascuorum.mkb -out pascuorum.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7
/Applications/MOSAIK_2.1.73/MosaikAligner -in pascuorum_s.mkb -out pascuorum_s.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7

rm pascuorum.mkb pascuorum_s.mkb

/Applications/MOSAIK_2.1.73/MosaikAligner -in terrestris.mkb -out terrestris.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7
/Applications/MOSAIK_2.1.73/MosaikAligner -in terrestris_s.mkb -out terrestris_s.mka -ia ref.dat -annpe /Applications/MOSAIK_2.1.73/2.1.26.pe.100.0065.ann -annse /Applications/MOSAIK_2.1.73/2.1.26.se.100.005.ann -statmq 30 -p 7

rm terrestris.mkb terrestris_s.mkb

samtools merge -@ 6 lucorummerged.bam lucorum.mka.bam lucorum_s.mka.bam
samtools merge -@ 6 terrestrismerged.bam terrestris.mka.bam terrestris_s.mka.bam
samtools merge -@ 6 pascuorummerged.bam pascuorum.mka.bam pascuorum_s.mka.bam
samtools merge -@ 6 bombusmerged.bam bombus.mka.bam bombus_s.mka.bam

samtools sort -@ 6 bombusmerged.bam -o bombusmergedsorted.bam
samtools sort -@ 6 lucorummerged.bam -o lucorummergedsorted.bam
samtools sort -@ 6 terrestrismerged.bam -o terrestrismergedsorted.bam
samtools sort -@ 6 pascuorummerged.bam -o pascuorummergedsorted.bam

rm bombusmerged.bam lucorummerged.bam terrestrismerged.bam pascuorummerged.bam

samtools index -@ 6 bombusmergedsorted.bam
samtools index -@ 6 lucorummergedsorted.bam
samtools index -@ 6 terrestrismergedsorted.bam
samtools index -@ 6 pascuorummergedsorted.bam

samtools view -@ 6 -b bombusmergedsorted.bam Loch_Morlich_virus Mayfield_virus_1 Mayfield_virus_2 River_Liunaeg_virus AF150629.1_Acute_bee_paralysis_virus,_complete_genome EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds > bombusmergedsortedreduced.bam
samtools view -@ 6 -b lucorummergedsorted.bam Loch_Morlich_virus Mayfield_virus_1 Mayfield_virus_2 River_Liunaeg_virus AF150629.1_Acute_bee_paralysis_virus,_complete_genome EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds > lucorummergedsortedreduced.bam
samtools view -@ 6 -b terrestrismergedsorted.bam Loch_Morlich_virus Mayfield_virus_1 Mayfield_virus_2 River_Liunaeg_virus AF150629.1_Acute_bee_paralysis_virus,_complete_genome EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds > terrestrismergedsortedreduced.bam
samtools view -@ 6 -b pascuorummergedsorted.bam Loch_Morlich_virus Mayfield_virus_1 Mayfield_virus_2 River_Liunaeg_virus AF150629.1_Acute_bee_paralysis_virus,_complete_genome EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds > pascuorummergedsortedreduced.bam

samtools merge -@ 6 allmergedreduced.bam bombusmergedsortedreduced.bam lucorummergedsortedreduced.bam terrestrismergedsortedreduced.bam pascuorummergedsortedreduced.bam

samtools sort -@ 6 allmergedreduced.bam -o allmergedreducedsorted.bam

samtools index allmergedreducedsorted.bam

~/Downloads/lofreq_star-2.1.2/bin/lofreq call -f ReftrimmedtoRdRpsdisamb.fasta -o vars.vcf allmergedreducedsorted.bam

java -jar /Users/davidpascall/Downloads/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=ReftrimmedtoRdRpsdisamb.fasta OUTPUT=ReftrimmedtoRdRpsdisamb.dict

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedsorted.bam -knownSites vars.vcf -o recal_data.table

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedsorted.bam -knownSites vars.vcf -BQSR recal_data.table -o post_recal_data.table 

java -jar ~/Desktop/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ReftrimmedtoRdRpsdisamb.fasta -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf

java -jar ~/Desktop/GenomeAnalysisTK.jar -T PrintReads -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedsorted.bam -BQSR recal_data.table -o allmergedreducedrecal1.bam

samtools sort -@ 6 allmergedreducedrecal1.bam -o allmergedreducedrecal1sorted.bam

samtools index allmergedreducedrecal1sorted.bam

~/Downloads/lofreq_star-2.1.2/bin/lofreq call -f ReftrimmedtoRdRpsdisamb.fasta -o vars2.vcf allmergedreducedrecal1sorted.bam

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal1sorted.bam -knownSites vars2.vcf -o recal_data2.table

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal1sorted.bam -knownSites vars2.vcf -BQSR recal_data2.table -o post_recal_data2.table 

java -jar ~/Desktop/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ReftrimmedtoRdRpsdisamb.fasta -before recal_data2.table -after post_recal_data2.table -plots recalibration_plots2.pdf

java -jar ~/Desktop/GenomeAnalysisTK.jar -T PrintReads -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal1sorted.bam -BQSR recal_data2.table -o allmergedreducedrecal2.bam

samtools sort -@ 6 allmergedreducedrecal2.bam -o allmergedreducedrecal2sorted.bam

samtools index allmergedreducedrecal2sorted.bam

~/Downloads/lofreq_star-2.1.2/bin/lofreq call -f ReftrimmedtoRdRpsdisamb.fasta -o vars3.vcf allmergedreducedrecal2sorted.bam

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal2sorted.bam -knownSites vars3.vcf -o recal_data3.table

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal2sorted.bam -knownSites vars3.vcf -BQSR recal_data3.table -o post_recal_data3.table 

java -jar ~/Desktop/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ReftrimmedtoRdRpsdisamb.fasta -before recal_data3.table -after post_recal_data3.table -plots recalibration_plots3.pdf

java -jar ~/Desktop/GenomeAnalysisTK.jar -T PrintReads -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal2sorted.bam -BQSR recal_data3.table -o allmergedreducedrecal3.bam

samtools sort -@ 6 allmergedreducedrecal3.bam -o allmergedreducedrecal3sorted.bam

samtools index allmergedreducedrecal3sorted.bam

~/Downloads/lofreq_star-2.1.2/bin/lofreq call -f ReftrimmedtoRdRpsdisamb.fasta -o vars4.vcf allmergedreducedrecal3sorted.bam

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal3sorted.bam -knownSites vars4.vcf -o recal_data4.table

java -jar ~/Desktop/GenomeAnalysisTK.jar -T BaseRecalibrator -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal3sorted.bam -knownSites vars4.vcf -BQSR recal_data4.table -o post_recal_data4.table 

java -jar ~/Desktop/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ReftrimmedtoRdRpsdisamb.fasta -before recal_data4.table -after post_recal_data4.table -plots recalibration_plots4.pdf

java -jar ~/Desktop/GenomeAnalysisTK.jar -T PrintReads -R ReftrimmedtoRdRpsdisamb.fasta -I allmergedreducedrecal3sorted.bam -BQSR recal_data4.table -o allmergedreducedrecal4.bam

samtools sort -@ 6 allmergedreducedrecal4.bam -o allmergedreducedrecal4sorted.bam

samtools index allmergedreducedrecal4sorted.bam

~/Downloads/lofreq_star-2.1.2/bin/lofreq call -f ReftrimmedtoRdRpsdisamb.fasta -o vars5.vcf allmergedreducedrecal4sorted.bam
