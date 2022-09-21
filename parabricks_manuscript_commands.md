## Bioinformatic workflow commands

### 1. _Germline workflows_ 

GATK Haplotypecaller:
```
bwa mem -t 32 -Y -K 10000000 -r $READGROUP $REF $FORWARD $REVERSE | gatk SortSam -I out.sam -O out.bam -SO coordinate
gatk MarkDuplicates I=out.bam O=dup.bam M=metrics.txt
gatk BaseRecalibrator -I dup.bam -O bqsr_out.txt --known-sites Homo_sapiens_assembly38.known_indels.vcf.gzs -R $REF
gatk ApplyBQSR -R $REF -I dup.bam -O combined.bam -bqsr-recal-file bqsr_out.txt
gatk HaplotypeCaller -I combined.bam -o combined.vcf -r $REF ---native-pair-hmm-threads 32 --output-mode EMIT_VARIANTS_ONLY
```

DeepVariant:
```
bwa mem -t 32 -Y -K 10000000 -r $READGROUP $REF $FORWARD $REVERSE | gatk SortSam -I out.sam -O out.bam -SO coordinate
gatk MarkDuplicates I=out.bam O=dup.bam M=metrics.txt
gatk BuildBamIndex -I dup.bam -O dup.bai

#run as a shell script
INPUT_DIR="${PWD}/data/dedup/"
REF_DIR="${PWD}/data/reference"
REF="Homo_sapiens_assembly38.fasta"
BAM="${S}_dedup.bam"
OUTPUT_DIR="${PWD}/data/calls/"
OUTPUT_VCF="${S}.dv.vcf.gz"
nproc=32

sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}":"/output" \
    -v "${REF_DIR}":"/refdir" \
    google/deepvariant:"1.1.0"  \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/refdir/${REF}" \
    --reads="/input/${BAM}" \
    --output_vcf="/output/${OUTPUT_VCF}" \
    --num_shards=${nproc} \
    --intermediate_results_dir /output/intermediate_results_dir
```

### 2. _Somatic workflows_
 
Muse:
```
sudo ./MuSE call -f {input.ref} -O {params.prefix} -n 1 {input.T} {input.N} |& tee logs/muse_v1.log
sudo ./MuSE sump -I {params.prefix} -G -O muse/muse_v1.vcf -D {input.dref}
```

Mutect2:
```
gatk Mutect2 \
-R {input.fasta} \
--input {input.T} \
--tumor-sample HG002  \
--input {input.N} \
--normal-sample HG002 \
--output {output} \
--native-pair-hmm-threads {THREADS}")
```

LoFreq:

`lofreq somatic -n {input.N} -t {input.T} -o lofreq/ -f {input.ref} --threads {THREADS} -d {input.dref} --baq-off --no-src-qual --call-rlx-extra-args '@d 2147483647' |& tee logs/lofreq.log`

SomaticSniper:
```
bam-somaticsniper -q {params.q} -G -L -F vcf -f {input.ref} {input.T} {input.N} {output} 2> {log}
bcftools mpileup -A -B -d {params.d} -Ou -f {input.ref} {input.T} | bcftools call -c | perl /home/usa_roscampbell_deloitte_com/miniconda3/envs/snakemake-tutorial/bin/vcfutils.pl varFilter -Q 20 | awk NR > 55 {{print}} > var_calling/somaticsniper/output.indel_pileup_Tum.pileup 2> {log}
#cp {output} var_calling/somaticsniper/backup.vcf 2> {log}
perl somatic-sniper/src/scripts/snpfilter.pl --snp-file {output} --indel-file var_calling/somaticsniper/output.indel_pileup_Tum.pileup 2> {log})
perl somatic-sniper/src/scripts/prepare_for_readcount.pl --snp-file {output}.SNPfilter 2> {log}
bam-readcount/build/bin/bam-readcount -b 15 -f {input.ref} -l {output}.SNPfilter.pos {input.T} > var_calling/somaticsniper/output.readcounts.rc 2> {log}
perl somatic-sniper/src/scripts/fpfilter.pl -snp-file {output}.SNPfilter -readcount-file var_calling/somaticsniper/output.readcounts.rc 2> {log}
cp {output}.SNPfilter.fp_pass {output}.SNPfilter.fp_pass.vcf 2> {log}
perl somatic-sniper/src/scripts/highconfidence.pl -snp-file {output}.SNPfilter.fp_pass.vcf 2> {log}
```

Strelka:
```
mkdir -p var_calling/manta_work
/usr/bin/python2.7 /home/usa_roscampbell_deloitte_com/manta-1.6.0.centos6_x86_64/bin/configManta.py --referenceFasta {input.ref} --normalBam {input.N} --tumorBam {input.T} --runDir var_calling/manta_work 2> {log}
cd var_calling/manta_work
/usr/bin/python2.7 /home/usa_roscampbell_deloitte_com/varcalling/var_calling/manta_work/runWorkflow.py -m local -j 32 2> {log}
cd ..
mkdir -p var_calling/strelka_work
/usr/bin/python2.7 strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --referenceFasta {input.ref} --normalBam {input.N} --tumorBam {input.T} --indelCandidates var_calling/manta_work/results/variants/candidateSmallIndels.vcf.gz --runDir var_calling/strelka_work 2> {log}
cd var_calling/strelka_work
/usr/bin/python2.7 varcalling/var_calling/strelka_work/runWorkflow.py -m local -j 32 2> {log}
```
