#!/bin/bash

nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate1_S1_L001_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate1_S1_L001_R2_001.fastq.gz \
--metadata plate1.csv --read-structures 8B12S+T 8S+T \
--output test/Plate1_L001 --metrics Plate1_L001_demux2_metrics.txt --threads 10 --output-type Fastq &

nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate2_S2_L001_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate2_S2_L001_R2_001.fastq.gz \
--metadata plate2.csv --read-structures 8B12S+T 8S+T \
--output Plate2_L001 --metrics Plate2_L001_demux2_metrics.txt --threads 10 --output-type Fastq &

nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate3_S3_L001_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate3_S3_L001_R2_001.fastq.gz \
--metadata plate3.csv --read-structures 8B12S+T 8S+T \
--output Plate3_L001 --metrics Plate3_L001_demux2_metrics.txt --threads 10 --output-type Fastq &


nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate1_S1_L002_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate1_S1_L002_R2_001.fastq.gz \
--metadata plate1.csv --read-structures 8B12S+T 8S+T \
--output Plate1_L002 --metrics Plate1_L002_demux2_metrics.txt --threads 10 --output-type Fastq &

nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate2_S2_L002_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate2_S2_L002_R2_001.fastq.gz \
--metadata plate2.csv --read-structures 8B12S+T 8S+T \
--output Plate2_L002 --metrics Plate2_L002_demux2_metrics.txt --threads 10 --output-type Fastq &

nohup fgbio -Xmx16G DemuxFastqs --inputs /bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate3_S3_L002_R1_001.fastq.gz \
/bulk/instrument_data/imc_mg/2023/imp_syc_dans/novaseq_run/data/Plate3_S3_L002_R2_001.fastq.gz \
--metadata plate3.csv --read-structures 8B12S+T 8S+T \
--output Plate3_L002 --metrics Plate3_L002_demux2_metrics.txt --threads 10 --output-type Fastq &
