# 멀티오믹스 암 게놈 분석 파이프라인 설치 가이드

## 📋 목차
1. [기본 환경 설치](#1-기본-환경-설치)
2. [Conda 환경 생성](#2-conda-환경-생성)
3. [컨테이너 이미지 다운로드](#3-컨테이너-이미지-다운로드)
4. [리소스 파일 준비](#4-리소스-파일-준비)
5. [설치 확인](#5-설치-확인)

## 1. 기본 환경 설치

### 필수 요구사항
- **OS**: Linux (Ubuntu 18.04+ 또는 CentOS 7+ 권장)
- **메모리**: 최소 64GB RAM (128GB+ 권장)
- **저장공간**: 최소 2TB (참조 데이터 + 분석 결과)
- **CPU**: 48+ cores 권장

### Conda/Mamba 설치
```bash
# Miniforge 설치 (mamba 포함)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

# 또는 Miniconda 설치 후 mamba 설치
conda install -c conda-forge mamba
```

## 2. Conda 환경 생성

```bash
# 1. 파이프라인 디렉토리로 이동
cd snakemake_v1.1_250605_denovo

# 2. Conda 환경 생성 (시간이 오래 걸릴 수 있음)
mamba env create -f env.yaml

# 3. 환경 활성화
conda activate long_read_pipeline

# 4. 설치 확인
snakemake --version
samtools --version
```

## 3. 컨테이너 이미지 다운로드

### 필요한 Singularity 이미지들

```bash
# 컨테이너 디렉토리 생성
mkdir -p resource/container

# Clair3 이미지
singularity pull resource/container/clair3_latest.sif docker://hkubal/clair3:latest

# DeepSomatic 이미지  
singularity pull resource/container/deepsomatic_1.8.0.sif docker://google/deepsomatic:1.8.0

# SAVANA 이미지
singularity pull resource/container/savana_latest.sif docker://quay.io/biocontainers/savana

# Wakhan 이미지
singularity pull resource/container/wakhan_latest.sif docker://mkolmogo/wakhan:dev_c717baa
```

### 이미지 용량 확인
```bash
ls -lh resource/container/
# 예상 총 용량: ~10-15GB
```

## 4. 리소스 파일 준비

### 4.1 참조 게놈 (필수)
```bash
mkdir -p resource/ref

# GRCh38 참조 게놈 다운로드
wget -O resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"

gunzip resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz

# 인덱스 생성
samtools faidx resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta
```

### 4.2 VEP 캐시 (필수)
```bash
mkdir -p resource/vep_cache

# VEP 캐시 다운로드 (약 15GB)
wget -O resource/vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz \
  "https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz"
```

### 4.3 AnnotSV 캐시 (필수)
```bash
mkdir -p resource/annotsv

# AnnotSV 주석 데이터 다운로드
wget -O resource/annotsv/annotsv_cache.tar.gz \
  "https://www.lbgi.fr/~geoffroy/Annotations/AnnotSV_annotations_3.4.tar.gz"
```

### 4.4 기타 필수 리소스
```bash
# VNTR BED 파일
mkdir -p resource/severus
wget -O resource/severus/human_GRCh38_no_alt_analysis_set.trf.bed \
  "https://github.com/KolmogorovLab/Severus/raw/main/resources/human_GRCh38_no_alt_analysis_set.trf.bed"

# Contig BED 파일
mkdir -p resource/tabix
echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY" > resource/tabix/chr.bed

# 암 유전자 목록
mkdir -p resource/intogen_genelist
# IntOGen Compendium Cancer Genes 파일 필요 (사용자가 직접 다운로드)
```

## 5. 설치 확인

### 5.1 기본 도구 확인
```bash
# 환경 활성화
conda activate long_read_pipeline

# 주요 도구들 확인
pbmm2 --version
samtools --version
bcftools --version
hiphase --version
severus --version
aligned_bam_to_cpg_scores --version
Rscript --version
```

### 5.2 컨테이너 확인
```bash
# 컨테이너 이미지 확인
singularity run resource/container/clair3_latest.sif /opt/bin/run_clair3.sh --version
singularity run resource/container/deepsomatic_1.8.0.sif run_deepsomatic --version
```

### 5.3 R 패키지 확인
```bash
Rscript -e "library(DSS); library(annotatr); sessionInfo()"
```

### 5.4 드라이런 테스트
```bash
# 샘플 데이터로 드라이런 실행
snakemake -n --configfile config.yaml
```

## 📝 주의사항

1. **메모리 사용량**: 일부 단계(특히 매핑과 변이 검출)에서 많은 메모리를 사용합니다.
2. **디스크 공간**: 중간 파일들이 많이 생성되므로 충분한 저장공간이 필요합니다.
3. **실행 시간**: 전체 파이프라인 실행에는 하루 이상 소요될 수 있습니다.
4. **네트워크**: 리소스 파일 다운로드에 안정적인 인터넷 연결이 필요합니다.

## 🆘 문제 해결

### 일반적인 문제들
1. **conda 패키지 충돌**: `mamba env create --force -f env.yaml`로 재설치
2. **컨테이너 오류**: Singularity 버전 확인 및 권한 설정
3. **메모리 부족**: config.yaml에서 스레드 수 조정
4. **R 패키지 오류**: `conda install -c bioconda r-*` 개별 설치

### 도움 요청
- GitHub Issues
- Snakemake 공식 문서: https://snakemake.readthedocs.io/ 