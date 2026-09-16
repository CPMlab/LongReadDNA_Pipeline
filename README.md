# LONG_READ_DNA_WGS

PacBio HiFi long-read WGS 기반 **암 게놈 (tumor/normal) 통합 분석 Snakemake 파이프라인**.
SNV/indel · 구조변이 · copy number · 메틸화/DMR 까지 한 번에 돌리고, 환자별로 결과를 정리한다.

> A Snakemake workflow for PacBio HiFi whole-genome tumor/normal cancer analysis:
> SNV/indel calling, structural variants, copy number, CpG methylation and DMR, with per-patient outputs.

## 이 파이프라인의 출발점

PacBio 공식 워크플로우 **[HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL)**
(tumor-only / matched tumor-normal somatic variant calling for HiFi reads)를 참고해서, 같은 분석 설계를
우리 클러스터(SLURM + conda + Singularity) 환경에 맞게 **Snakemake 로 재구성한 것**이다.

원본에서 그대로 가져온 설계:

- 체세포 SNV/indel: DeepSomatic → VEP 주석 → IntOGen Cancer Gene(CCG) 교차
- 구조변이: Severus → 필터링 → AnnotSV 주석 → IntOGen CCG 교차
- Mitelman database 기반 fusion 표시
- CNV/purity·ploidy: Wakhan
- 메틸화: pb-CpG-tools → DSS 로 tumor vs normal DMR

우리 쪽에서 달라진 점:

| 항목 | HiFi-somatic-WDL | 이 저장소 |
|---|---|---|
| 실행 엔진 | WDL (miniwdl / Cromwell) | Snakemake + SLURM 제출 스크립트 |
| 종양 샘플 | tumor-only 또는 tumor/normal 1쌍 | 환자당 **종양 샘플 여러 개** (PRIMARY/META 등)를 Severus multimode 로 동시 분석 |
| 샘플 관리 | 샘플별 입력 JSON | `samples.tsv` 한 장으로 다환자·다샘플 일괄 |
| CNV | cnvkit / PURPLE / Wakhan | SAVANA + Wakhan |
| 시각화 | Severus cluster plot (HTML) | 자체 circos 스크립트(`scripts/circosplot.py`)로 SV·fusion circos |
| 생식세포 변이 | - | Clair3 + HiPhase 페이징 결과를 함께 산출 |

원본 워크플로우의 라이선스와 인용은 [HiFi-somatic-WDL 저장소](https://github.com/PacificBiosciences/HiFi-somatic-WDL)를 따른다.

---

## 특징

- **samples.tsv 한 장으로 제어** — 환자 수, 샘플 타입, BAM 경로만 적으면 나머지는 자동.
- **종양 샘플 타입 자유** — `TUMOR/NORMAL` 페어든 `PRIMARY/META` 다중 종양이든 그대로 동작한다.
  `Snakefile` 이 samples.tsv 에서 종양 타입을 읽어 wildcard 제약(`TUMOR_TYPE_CONSTRAINT`)을 자동 생성하므로,
  데이터셋마다 파이프라인을 복제할 필요가 없다.
- **Severus multimode** — NORMAL 을 control 로 두고 여러 종양 샘플을 동시에 SV 분석.
- **여러 환자 동시 처리** — 환자별 디렉토리로 결과가 분리된다. 새 환자는 samples.tsv 에 줄만 추가.
- 같은 `patient_id + sample_type` 이 여러 줄이면(= 여러 SMRT cell) 자동으로 병합한다.

## 워크플로우

| 단계 | 도구 | rule |
|---|---|---|
| 매핑 | pbmm2 | `mapping` |
| 커버리지 QC | mosdepth | `mosdepth` |
| 생식세포 변이 | Clair3 (container) | `clair3` |
| 체세포 변이 | DeepSomatic (container) | `deepsomatic` |
| 페이징 | HiPhase | `hiphase_normal`, `hiphase_tumor` |
| VCF 정규화 | bcftools | `normalize_*_vcf` |
| 변이 주석 | VEP | `vep_annotate_normal`, `vep_annotate_somatic` |
| 구조변이 | Severus (multimode) | `severus_multimode` |
| SV 필터/후처리 | SVpack, mate-BND 복구 | `tabix_filter`, `svpack`, `recover_mate_bnd` |
| SV 주석 | AnnotSV + IntOGen CCG | `annotsv`, `sv_intogen` |
| CNV / purity·ploidy | SAVANA, Wakhan (container) | `savana_sv`, `wakhan_cnv` |
| 메틸화 | pb-CpG-tools | `cpg_methylation` |
| 차등 메틸화 | DSS + annotatr | `dss_dmr`, `annotate_dmr` |
| 시각화 | circosplot (자체 스크립트) | `circosplot`, `collect_sample_circos` |

## 빠른 시작

```bash
# 1. 환경 (자세한 내용은 docs/INSTALLATION.md)
mamba env create -f env.yaml
conda activate long_read_pipeline

# 2. 샘플 정보 작성
cp samples.example.tumor_normal.tsv samples.tsv   # 또는 samples.example.multi_tumor.tsv
vi samples.tsv

# 3. 참조 데이터 경로 확인
vi config.yaml

# 4. dry-run 으로 계획 확인
snakemake -n --cores 8

# 5. 실행 (SLURM)
sbatch slurm/run_snakemake.sh
# 특정 환자만
bash slurm/run_patient.sh PT001 96
```

## samples.tsv 형식

```tsv
patient_id	sample_type	bam_path
PT001	NORMAL	/path/to/PT001_N.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_1.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_2.hifi_reads.bam
```

- `sample_type` 에서 **`NORMAL` 은 예약어**(control), 나머지는 전부 종양 샘플로 취급된다.
- 따라서 `TUMOR`, `PRIMARY`, `META`, `RELAPSE` 등 원하는 이름을 쓰면 되고, 환자마다 NORMAL 은 반드시 있어야 한다.
- 실제 `samples.tsv` 는 환자 식별 정보를 담을 수 있어 `.gitignore` 로 커밋에서 제외되어 있다. 예시 파일만 저장소에 포함된다.

## 디렉토리 구조

```
LONG_READ_DNA_WGS/
├── Snakefile                        # 메인 워크플로우 (샘플 로딩 + rule all)
├── config.yaml                      # 경로/파라미터 설정
├── env.yaml                         # conda 환경
├── samples.example.*.tsv            # 샘플 시트 예시 2종
├── rules/                           # 단계별 rule (mapping, variant_calling, ...)
├── scripts/                         # 자체 스크립트 (circosplot.py, DSS/annotatr R, svpack.py)
├── resources/                       # 소형 참조 파일 (chr.bed, hg38.bed, cytoband)
├── slurm/                           # SLURM 제출 스크립트
└── docs/                            # 설치 · 사용법 · 운영 가이드
```

## 별도로 준비해야 하는 대용량 리소스

저장소에는 코드와 소형 참조 파일만 들어 있다. 아래는 `resource/` 아래에 직접 준비해야 하며(`.gitignore` 대상),
설치 절차는 [docs/INSTALLATION.md](docs/INSTALLATION.md) 참조.

| config 키 | 내용 | 대략 크기 |
|---|---|---|
| `reference_fasta` | GRCh38 no-alt analysis set | ~3 GB |
| `vep_cache` | VEP 112 RefSeq 캐시 | ~26 GB |
| `annotsv_cache` | AnnotSV annotation | ~5 GB |
| `*_container` | Clair3 / DeepSomatic / SAVANA / Wakhan 이미지 | ~11 GB |
| `vntr_bed` | Severus VNTR BED | ~7 MB |
| `svpack_match_vcf`, `reference_gff` | SVpack 대조 VCF, GFF3 | ~700 MB |
| `compendium_file` | IntOGen Compendium Cancer Genes | ~1 MB |
| `mitelman_mcgene` | Mitelman DB MCGENE 덤프 (fusion 표시용) | ~4 MB |

## 입력 데이터 규모 (참고)

입력은 PacBio Revio 에서 나온 **unaligned BAM (`*.hifi_reads.bam`, uBAM)** 이다. 같은 샘플이 여러 SMRT cell 로
나뉘어 있으면 `samples.tsv` 에 줄을 여러 개 적으면 `mapping` rule 이 합쳐서 정렬한다.

실제 운영 중인 유방암 코호트(WGS, 30x 내외) 기준 실측치:

| 항목 | 샘플당 | 비고 |
|---|---|---|
| uBAM (`*.hifi_reads.bam`) | **약 200~370 GB / SMRT cell** | 종양은 보통 2 cell → 600~700 GB |
| HiFi fastq.gz | 약 30~50 GB / cell | uBAM 이 있으면 파이프라인에 불필요 |
| fail_reads BAM | 약 40~60 GB | 분석에 사용하지 않음 |
| 정렬 후 haplotagged BAM | 약 240~760 GB | 파이프라인 산출물 |

- 환자 1명(정상 1 + 종양 1) 기준 **입력 uBAM 만 약 1 TB**, 중간·최종 산출물까지 합치면 **2~3 TB** 의 작업 공간이 필요하다.
- 종양이 원발+전이 2개인 환자는 입력만 1.4~1.7 TB 수준.
- 14 샘플 코호트 실측 합계: uBAM 약 7.0 TB, 전체(raw + 납품 BAM 포함) 약 16.8 TB.

## 검증 이력

유방암 long-read 코호트(정상–종양 페어 3명, 원발–전이 다중 종양 2명)에서 전 단계 완주로 검증했다.
`TUMOR/NORMAL` 구성과 `PRIMARY/META` 구성 두 갈래로 따로 유지하던 파이프라인을 이 저장소에서 한 벌로 합쳤고,
통합 후 Snakemake 7.32.4 dry-run 으로 두 구성 모두 DAG 가 정상 생성되는 것을 확인했다.

| 구성 | samples.tsv | 생성된 job |
|---|---|---|
| 정상 1 + 종양 1 | `samples.example.tumor_normal.tsv` 형식 | 28 |
| 정상 1 + 원발 + 전이 | `samples.example.multi_tumor.tsv` 형식 | 40 |

```bash
# 실제로 돌리기 전에 항상 이걸로 확인
snakemake -n --cores 4
```

## 알려진 이슈 / 패치 노트

- **SAVANA v1.3.0**: 옵션명이 `--phased_vcf` → `--snp_vcf` 로 변경됨. rule 이 새 옵션명을 사용한다.
- **pb-CpG-tools v3.0.0**: 출력에 헤더 라인이 생겨 DSS 입력에서 깨졌다. 헤더를 제거하고 넘기는 방식으로 정리했다.
- DMR 주석 산출물은 파일명이 버전에 따라 달라져서 glob 으로 찾아 이동시킨다.

## 외부 코드 / 데이터 출처

- `scripts/svpack.py` — PacBio [svpack](https://github.com/PacificBiosciences/svpack) 에서 가져옴.
- `scripts/DSS_tumor_normal.R`, `scripts/annotatr_dmr.R`, `scripts/circosplot.py` — 본 파이프라인 자체 스크립트.
- Mitelman Database, IntOGen Compendium 은 각 배포처 라이선스를 따르며 저장소에 포함하지 않는다.

## 문서

- [docs/USAGE.md](docs/USAGE.md) — 샘플 시트 작성, 단계별 실행, 출력 구조
- [docs/OPERATIONS.md](docs/OPERATIONS.md) — 실행/모니터링/트러블슈팅
- [docs/INSTALLATION.md](docs/INSTALLATION.md) — 환경 · 컨테이너 · 참조 데이터 설치
