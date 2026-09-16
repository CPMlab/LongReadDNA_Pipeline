# 사용법 (범용 변이 검출 워크플로우)

TSV 파일에서 샘플 정보를 읽어들여 자동으로 여러 환자의 변이 검출을 수행하는 Snakemake 파이프라인입니다.

## 주요 특징

- **TSV 기반 샘플 관리**: 하드코딩된 샘플 정보 대신 TSV 파일로 관리
- **Severus Multimode 지원**: PRIMARY, META 등 여러 종양 샘플을 동시에 분석
- **다중 환자 지원**: 하나의 파이프라인으로 여러 환자 동시 처리
- **완전 자동화**: 매핑부터 주석까지 전체 파이프라인 자동화

## 파일 구조

```
LONG_READ_DNA_WGS/
├── Snakefile      # 메인 Snakefile
├── config.yaml    # 설정 파일
├── samples.tsv             # 샘플 정보 TSV 파일
└── docs/USAGE.md            # 사용법 (이 파일)
```

## samples.tsv 형식

TSV 파일은 다음과 같은 형식으로 작성해야 합니다:

```tsv
patient_id	sample_type	bam_path
PT101	NORMAL	/path/to/PT101_N.hifi_reads.bam
PT101	PRIMARY	/path/to/PT101_Primary.hifi_reads.bam
PT101	META	/path/to/PT101_Meta_1.hifi_reads.bam
PT101	META	/path/to/PT101_Meta_2.hifi_reads.bam
PT101	META	/path/to/PT101_Meta.hifi_reads.bam
PT102	NORMAL	/path/to/PT102_N.hifi_reads.bam
PT102	PRIMARY	/path/to/PT102_Primary.hifi_reads.bam
PT102	META	/path/to/PT102_Meta.hifi_reads.bam
```

### 필수 컬럼:
- **patient_id**: 환자 식별자
- **sample_type**: 샘플 타입 (NORMAL, PRIMARY, META 등)
- **bam_path**: BAM 파일의 절대 경로

### 주의사항:
- 각 환자는 반드시 **NORMAL** 샘플을 가져야 함
- PRIMARY, META 등은 종양 샘플로 취급됨
- 같은 patient_id + sample_type 조합이 여러 개 있으면 자동으로 병합됨

## 실행 방법

### 1. 기본 실행
```bash
snakemake --configfile config.yaml --cores 48
```

### 2. 특정 환자만 실행
```bash
snakemake --configfile config.yaml --cores 48 \
    results/PT101/sv/PT101.sv.annotsv_intogenCCG.tsv
```

### 3. 특정 단계까지만 실행 (예: 매핑까지)
```bash
snakemake --configfile config.yaml --cores 48 \
    results/PT101/mapping/PT101.NORMAL.aligned.bam \
    results/PT101/mapping/PT101.PRIMARY.aligned.bam \
    results/PT101/mapping/PT101.META.aligned.bam
```

### 4. Dry-run (실제 실행하지 않고 계획만 확인)
```bash
snakemake --configfile config.yaml --cores 48 -n
```

## 출력 구조

```
results/
├── PT101/                    # 환자별 디렉토리
│   ├── mapping/           # 매핑 결과
│   │   ├── PT101.NORMAL.aligned.bam
│   │   ├── PT101.PRIMARY.aligned.bam
│   │   └── PT101.META.aligned.bam
│   ├── qc/               # QC 결과
│   ├── snv/              # SNV 검출 결과
│   ├── phasing/          # 위상 결정 결과
│   ├── annotation/       # VEP 주석 결과
│   ├── sv/               # SV 분석 결과
│   └── logs/             # 로그 파일들
└── PT102/                   # 다른 환자 디렉토리
    └── ...
```

## 주요 분석 단계

1. **매핑 (Mapping)**: pbmm2를 사용한 HiFi reads 매핑
2. **QC**: mosdepth를 사용한 커버리지 분석
3. **SNV 검출**: Clair3를 사용한 변이 검출
4. **체세포 변이**: DeepSomatic을 사용한 체세포 변이 검출
5. **위상 결정**: HiPhase를 사용한 위상 결정
6. **정규화**: bcftools를 사용한 VCF 정규화
7. **주석**: VEP를 사용한 변이 주석
8. **SV 분석**: Severus multimode를 사용한 구조적 변이 검출
9. **SV 후처리**: SVpack, AnnotSV, IntOGen을 사용한 SV 주석

## Severus Multimode

이 파이프라인의 핵심 특징은 Severus multimode 지원입니다:

- Normal 샘플을 control로 사용
- PRIMARY, META 등 모든 종양 샘플을 동시에 target으로 분석
- 한 번의 실행으로 모든 종양 샘플 간 비교 분석 가능

## 설정 수정

`config.yaml` 파일에서 다음을 수정할 수 있습니다:

- **samples_file**: 샘플 정보 TSV 파일 경로
- **threads**: 사용할 스레드 수
- **output_dir**: 결과 출력 디렉토리
- **리소스 파일 경로들**: 참조 게놈, VEP 캐시, 컨테이너 등

## 새로운 환자 추가

1. `samples.tsv` 파일에 새로운 환자의 샘플 정보 추가
2. 파이프라인 재실행 - 기존 결과는 그대로 유지되고 새로운 환자만 분석됨

## 문제 해결

### 일반적인 오류들:
- **"환자에게 NORMAL 샘플이 없습니다"**: samples.tsv에 NORMAL 샘플 추가 필요
- **"필수 컬럼이 없습니다"**: TSV 파일의 헤더 확인
- **BAM 파일 경로 오류**: samples.tsv의 bam_path가 올바른지 확인

### 로그 확인:
```bash
# 특정 환자의 로그 확인
ls results/PT101/logs/

# 특정 단계의 로그 확인
cat results/PT101/logs/severus_PT101_multimode.log
```

## 예시

PT101 환자의 NORMAL, PRIMARY, META 샘플을 분석하는 예시:

1. `samples.tsv` 작성:
```tsv
patient_id	sample_type	bam_path
PT101	NORMAL	/data/PT101_N.hifi_reads.bam
PT101	PRIMARY	/data/PT101_Primary.hifi_reads.bam
PT101	META	/data/PT101_Meta1.hifi_reads.bam
PT101	META	/data/PT101_Meta2.hifi_reads.bam
```

2. 실행:
```bash
snakemake --configfile config.yaml --cores 48
```

3. 결과 확인:
```bash
ls results/PT101/sv/PT101.sv.annotsv_intogenCCG.tsv
```

이제 하드코딩 없이 TSV 파일만 수정하여 어떤 환자든 분석할 수 있습니다! 



추가변경사항 
1. DSS_tumor_normal.R 이 pb-cpg-tools 가 head 를 남기도록 최근에 패치(V3.0.0)가 되어서 
    -> DSS_tumor_normal.R 에서 
    ```
    tumor  <- fread(tumor_file,  header = FALSE, sep = "\t", comment.char = "#")
    normal <- fread(normal_file, header = FALSE, sep = "\t", comment.char = "#")
    ```
    이렇게 바굼 원래는 sep = "\t" 까지만 있었음.
    -> 철회하고 다시 원복, 헤드 없애는걸로 함

2. savana 패치를 해서(v1.3.0) --phased_vcf 가 -> --snp_vcf 로 바뀜