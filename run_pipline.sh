### analysi pipeline
source activate hicexplorer
# HiC-Pro 稀疏矩阵 + BED → .cool（原始 & ICE）
hicConvertFormat -m sample_20kb.raw_contact_map.txt.gz \
  --inputFormat hicpro --bedFileHicpro genomic_intervals.20kb.bed \
  --outputFormat cool -o sample_20kb.raw.cool

hicConvertFormat -m sample_20kb.iced_contact_map.txt.gz \
  --inputFormat hicpro --bedFileHicpro genomic_intervals.20kb.bed \
  --outputFormat cool -o sample_20kb.iced.cool

# 用标准化后的 Hi-C（.cool）计算 TAD 分隔分数
hicFindTADs -m sample_20kb.iced.cool --outPrefix sample_20kb \
  --numberOfProcessors 8 --correctForMultipleTesting fdr
# 产物包括：
#   sample_20kb_tad_score.bm      ← bedgraph-matrix，可直接画出“分隔分数曲线”
#   sample_20kb_domains.bed       ← TAD 区域
#   sample_20kb_boundaries.bed    ← TAD 边界

### Chip-seq and ATAC-seq
  ###alignment 
bwa-mem2 mem -t 8 -R $'@RG\\tID:SRR13274560\\tSM:SRR13274560\\tPL:ILLUMINA' ~/data/hs37d5.fa SRR13274560_1.clean.fq.gz SRR13274560_2.clean.fq.gz | samtools sort -@ 8 -O BAM -o SRR13274560.hs37d5.sorted.bam -
  ###Remove_deplicates
java -jar ~/bin/picard.jar MarkDuplicates I=SRR13274568.hs37d5.sorted.bam O=SRR13274568.hs37d5.sorted.dedup.bam M=SRR13274568.metrics.txt REMOVE_DUPLICATES=true
  ###去除blacklist区域，降低假阳性
samtools view -b -q 30 SRR11016356.hs37d5.sorted.dedup.bam|bedtools intersect -v -abam stdin -b hs37d5-blacklist.bed > SRR11016356.filt.bam
samtools index SRR11016356.filt.bam
  ###bam2bw ###hg37 genome size: 2864785220 
source activate deeptools
bamCoverage -b SRR11016356.shift.sort.bam -o SRR11016356.RPGC.bw \
  -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2864785220
  ###ATAC-seq 需要Tn5 切割位点偏移校正（强烈推荐）ATAC-seq 常将 reads 按 Tn5 切割位点做 +4/-5 移位。deepTools 的 alignmentSieve --ATACshift 可直接对 BAM 做位移：
alignmentSieve -b SRR11016356.filt.bam --ATACshift -p 8 -o SRR11016356.shift.bam
samtools index SRR11016356.shift.bam
###生成/编辑绘图配置
source activate hicexplorer
make_tracks_file  --trackFiles GSM3905152_G_rep1_20kb.iced.cool sample_20kb_tad_score.bm ~/data/RNA-seq/Homo_sapiens.GRCh37.87.gtf data/SRR13274568.RPGC.bw data/SRR13274560.RPGC.bw data/SRR11016356.RPGC.bw  data/SRR9601022.RNA.hisat2.bw data/HC12.RNA.bw data/M21.RNA.bw data/MSN.RNA.bw data/SDD-1.RNA.bw data/SDD-2.RNA.bw data/SDD-3.RNA.bw  sample_20kb_boundaries.bed -o tracks.3.ini

pyGenomeTracks --tracks tracks.3.ini  --region chr2:31700000-33500000 --dpi 330 -o hic_multi_panel4.pdf
