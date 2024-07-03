 ####### 16S V4-V5 amplicon sequencing analysis pipeline #######

mkdir -p Analysis_TEMP/Cluster_TEMP  rawdata/qc cleandata Analysis_TEMP/tree
cd rawdata
conda activate fastp
for i in $(ls *.R1.fq.gz); do i=${i/.R1.fq.gz/}; \
nohup fastp -i ${i}.R1.fq.gz -o ../cleandata/${i}_1.fastp.fq.gz \
-I ${i}.R2.fq.gz -O ../cleandata/${i}_2.fastp.fq.gz \
-l 40 -q 20 —compression=6 -R ${i}  \
-h qc/${i}.fastp.html -j qc/${i}.fastp.json; done

cd ..
conda deactivate
conda activate qiime2-2020.6
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.txt \
--output-path Analysis_TEMP/Cluster_TEMP/paired-seqs.qza \
--input-format PairedEndFastqManifestPhred33

# summarize
qiime demux summarize --i-data Analysis_TEMP/Cluster_TEMP/paired-seqs.qza \
--o-visualization Analysis_TEMP/Cluster_TEMP/paired-seqs.qzv

####remove adapter ####
time qiime cutadapt trim-paired \
--i-demultiplexed-sequences Analysis_TEMP/Cluster_TEMP/paired-seqs.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r CCGYCAATTYMTTTRAGTTT  \
--o-trimmed-sequences Analysis_TEMP/Cluster_TEMP/paired-end-demux.qza \
--verbose \
&> primer_trimming.log

####visualization
time qiime demux summarize \
--i-data Analysis_TEMP/Cluster_TEMP/paired-end-demux.qza \
--o-visualization Analysis_TEMP/Cluster_TEMP/demux.qzv



https://view.qiime2.org/


# export data
qiime tools export --input-path Analysis_TEMP/Cluster_TEMP/paired-seqs.qzv \
--output-path Analysis_TEMP/Quality

qiime dada2 denoise-paired \
--i-demultiplexed-seqs Analysis_TEMP/Cluster_TEMP/paired-end-demux.qza \
--p-trunc-len-f 220 \
--p-trunc-len-r 200 \
--o-representative-sequences Analysis_TEMP/Cluster_TEMP/rep-seqs.qza \
--o-table Analysis_TEMP/Cluster_TEMP/table.qza \
--o-denoising-stats Analysis_TEMP/Cluster_TEMP/denoising.qza \
--p-n-threads 10 \
--verbose

qiime feature-table summarize --i-table Analysis_TEMP/Cluster_TEMP/table.qza \
--o-visualization Analysis_TEMP/Cluster_TEMP/table.qzv \
--m-sample-metadata-file metadata.txt --quiet
# reference seqs visualization
qiime feature-table tabulate-seqs --i-data Analysis_TEMP/Cluster_TEMP/rep-seqs.qza \
--o-visualization Analysis_TEMP/Cluster_TEMP/rep-seqs.qzv

#V4V5
qiime feature-classifier classify-sklearn \
--i-classifier /home/penghaoran/16s_classifier/silva-138-ssu-nr99-515f-926r-classifier.qza \
--i-reads Analysis_TEMP/Cluster_TEMP/rep-seqs.qza \
--o-classification Analysis_TEMP/Cluster_TEMP/taxonomy.qza \
--p-n-jobs 3

######pre-processing
qiime feature-table filter-features --i-table Analysis_TEMP/Cluster_TEMP/table.qza \
--p-min-frequency 1 --p-min-samples 1 \
--o-filtered-table Analysis_TEMP/Cluster_TEMP/feature-frequency-filtered-table.qza

#######remove mitochondria,chloroplast
qiime taxa filter-table \
--i-table Analysis_TEMP/Cluster_TEMP/feature-frequency-filtered-table.qza \
--i-taxonomy Analysis_TEMP/Cluster_TEMP/taxonomy.qza \
--p-exclude mitochondria,chloroplast,Arthropoda \
--o-filtered-table Analysis_TEMP/Cluster_TEMP/table-with-phyla-no-metochondria-no-chloroplast.qza

## barplot
qiime taxa barplot \
  --i-table Analysis_TEMP/Cluster_TEMP/table-with-phyla-no-metochondria-no-chloroplast.qza \
  --i-taxonomy Analysis_TEMP/Cluster_TEMP/taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization Analysis_TEMP/Cluster_TEMP/taxa-bar-plots.qzv

## Phylogenetic tree
time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences Analysis_TEMP/Cluster_TEMP/rep-seqs.qza \
--o-alignment Analysis_TEMP/tree/aligned-rep-seqs.qza \
--o-masked-alignment Analysis_TEMP/tree/masked-aligned-rep-seqs.qza \
--o-tree Analysis_TEMP/tree/unrooted-tree.qza \
--o-rooted-tree Analysis_TEMP/tree/rooted-tree.qza




## Export table
mkdir phyloseq
qiime tools export \
--input-path Analysis_TEMP/Cluster_TEMP/table-with-phyla-no-metochondria-no-chloroplast.qza \
--output-path phyloseq

biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/otu_table.tsv \
--to-tsv
cd phyloseq; sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID/ASV ID/' otu_table.tsv
cd ../

qiime tools export \
--input-path Analysis_TEMP/Cluster_TEMP/taxonomy.qza \
--output-path phyloseq

qiime tools export \
--input-path Analysis_TEMP/Cluster_TEMP/rep-seqs.qza \
--output-path phyloseq

qiime tools export \
--input-path Analysis_TEMP/tree/unrooted-tree.qza \
--output-path phyloseq
cd phyloseq; mv tree.nwk unrooted_tree.nwk; cd ..

qiime tools export \
--input-path Analysis_TEMP/tree/rooted-tree.qza \
--output-path phyloseq
cd phyloseq; mv tree.nwk rooted_tree.nwk;cd .. 

#time qiime tools import \
#--type 'FeatureData[Sequence]' \
#--input-path rep-seqs.fasta \
#--output-path rep-seqs.qza

