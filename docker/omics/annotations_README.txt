Curated T2T-CHM13v2.0 (UCSC hs1) NCBI RefSeq annotations
========================================================

Source
  URL            https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz
  Last-Modified  2023-05-29 (UCSC server)
  md5 of .gz     bd7063efcbfed4774d9be42db8173a6e
  md5 of the decompressed GTF (t2t.ncbiRefSeq.gtf)  95c46ca9573b2ed12fd03eec5e949164
The image build verifies both md5 values and fails on any mismatch. The uncurated
GTF is not kept.

Curation (awk on the decompressed t2t.ncbiRefSeq.gtf; $12 is the transcript_id value)
  awk '$1 == "chrM" || $12 ~ /"N[MR]_/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.gtf
      md5 e5106937f58cefaac61483720ae44b74   chrM + NM_ and NR_ transcripts
  awk '$1 == "chrM" || $12 ~ /"NM_/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.proteincoding.gtf
      md5 39249e1d0b047822d321a22d18c40a8e   chrM + NM_ (protein coding) transcripts
  awk '($1 == "chrM" || $12 ~ /\"N[MR]_/) && $1 !~ /chrUn|_random/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.norandom.gtf
      md5 e5106937f58cefaac61483720ae44b74   as curated.gtf; no chrUn / *_random contigs present
The md5 values are of the UNCOMPRESSED files; the three .gtf.gz files here are
`gzip -9n` of them (zcat <file> | md5sum reproduces the value). See checksums.txt.

TxDb
  txdb.hs1.sqlite  built from the norandom GTF (env r, R 4.4.3):
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format="gtf", organism="Homo sapiens",
              dataSource="ncbiRefSeq excluding random contigs and keeping only NM and NR annotations")
    AnnotationDbi::saveDb(txdb, file="txdb.hs1.sqlite")
  Load with AnnotationDbi::loadDb("/opt/omics/annotations/t2t/txdb.hs1.sqlite").
