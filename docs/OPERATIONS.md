# 운영 가이드 (실행 · 모니터링 · 트러블슈팅)

## 📋 빠른 시작

### 1. 샘플 정보 준비

실제 BAM 파일 경로로 `samples.tsv`를 수정하세요:

```bash
cd LONG_READ_DNA_WGS
vi samples.tsv
```

예시:
```
patient_id	sample_type	bam_path
PT001	NORMAL	/실제/경로/PT001_normal.hifi_reads.bam
PT001	TUMOR	/실제/경로/PT001_tumor.hifi_reads.bam
PT002	NORMAL	/실제/경로/PT002_normal.hifi_reads.bam
PT002	TUMOR	/실제/경로/PT002_tumor.hifi_reads.bam
```

### 2. 설정 확인

resource 경로들이 올바른지 확인하세요:

```bash
vi config.yaml
```

주요 확인 사항:
- `reference_fasta`: 참조 게놈 경로
- `vep_cache`: VEP 캐시 경로  
- `*_container`: 컨테이너 이미지 경로

### 3. 파이프라인 실행

#### 전체 분석 실행
```bash
# SLURM 스케줄러 사용 (권장)
sbatch slurm/run_snakemake.sh

# 또는 직접 실행 (터미널이 끊어지지 않도록 주의)
bash slurm/run_snakemake.sh
```

#### 특정 환자만 분석
```bash
# 환자 PT001만 분석
bash slurm/run_patient.sh PT001

# 코어 수 지정 (기본 96)
bash slurm/run_patient.sh PT001 48
```

#### 테스트 실행 (실제로 실행하지 않고 계획만 확인)
```bash
snakemake -n --cores 96
```

## 🔍 진행 상황 모니터링

### 로그 확인
```bash
# SLURM 로그 확인
tail -f longread_wgs.o[JOB_ID]
tail -f longread_wgs.e[JOB_ID]

# 또는 실시간 로그
watch -n 10 'tail -20 longread_wgs.o[JOB_ID]'
```

### 결과 파일 확인
```bash
# 생성된 파일들 확인
find results/ -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.svg" | head -20

# 특정 환자 결과 확인
ls -la results/PT001/*/
```

### 진행률 확인
```bash
# 완료된 작업 수 확인
snakemake --summary | grep -c "done"

# 전체 작업 계획 확인
snakemake --dryrun | grep "Job counts"
```

## 📊 결과 파일 해석

### 주요 결과 파일들

#### 1. 변이 호출 결과
```
results/PT001/annotation/
├── PT001.NORMAL.germline.vep.vcf.gz    # 생식세포 변이 (유전적)
└── PT001.TUMOR.somatic.vep.vcf.gz      # 체세포 변이 (암 특이적)
```

#### 2. 구조변이 분석
```
results/PT001/sv/
├── PT001.sv.annotsv_intogenCCG.tsv     # 주석된 구조변이 목록
└── circos_PT001/
    └── PT001.sv.circos.svg             # 구조변이 시각화
```

#### 3. 카피수 변이
```
results/PT001/cnv/
├── savana_PT001_TUMOR/
│   ├── PT001.TUMOR.classified.somatic.vcf    # SAVANA CNV 결과
│   └── PT001.TUMOR_fitted_purity_ploidy.tsv  # 순도/배수성
└── wakhan_PT001_TUMOR/
    ├── PT001.TUMOR.copynumbers_segments.bed  # Wakhan CNV 결과
    └── purity_ploidy.tsv                       # 순도/배수성
```

#### 4. 메틸화 분석
```
results/PT001/methylation/
├── PT001.NORMAL.cpg.combined.bed.gz    # 정상 조직 메틸화
├── PT001.TUMOR.cpg.combined.bed.gz     # 종양 조직 메틸화
└── PT001.NORMAL.cpg.combined.bw        # BigWig 시각화 파일
```

#### 5. 차등 메틸화
```
results/PT001/dmr/
├── PT001.TUMOR_vs_NORMAL.DMR.tsv              # 차등 메틸화 영역
└── PT001.TUMOR_vs_NORMAL.annotated_DMR.tsv.gz # 주석된 DMR
```

## 🛠️ 문제 해결

### 일반적인 문제들

#### 1. 메모리 부족 오류
```bash
# 코어 수 줄이기
snakemake --cores 48  # 96에서 48로

# 또는 config.yaml에서 스레드 수 조정
threads: 48
threads_low: 12
```

#### 2. 디스크 공간 부족
```bash
# 중간 파일 정리
snakemake --delete-temp-output

# 결과만 남기고 임시 파일 삭제
find results/ -name "*.tmp" -delete
```

#### 3. 컨테이너 오류
```bash
# Singularity 캐시 정리
rm -rf ${SINGULARITY_CACHEDIR:-$HOME/singularity_cache}/*

# 컨테이너 다시 다운로드
snakemake --use-singularity --singularity-download-images
```

#### 4. 특정 작업 재실행
```bash
# 특정 환자의 변이 호출만 재실행
snakemake --forcerun annotation --cores 96 results/PT001/annotation/PT001.TUMOR.somatic.vep.vcf.gz

# 모든 차등 메틸화 분석 재실행
snakemake --forcerun dmr_analysis --cores 96
```

#### 5. 락 파일 문제
```bash
# 락 해제
snakemake --unlock

# 강제 재시작
snakemake --rerun-incomplete --cores 96
```

## 📈 성능 최적화

### 리소스 사용량 모니터링
```bash
# CPU 사용률
htop

# 메모리 사용률  
free -h

# 디스크 사용률
df -h
```

### 병렬 처리 최적화
```bash
# 최대 동시 작업 수 조정
snakemake --cores 96 --jobs 16  # 16개 작업 동시 실행

# 네트워크 파일시스템용 지연 시간 증가
snakemake --latency-wait 300  # 5분
```

## 💡 추가 팁

### 특정 분석만 실행
```bash
# 매핑만 실행
snakemake --cores 96 results/PT001/mapping/PT001.TUMOR.aligned.bam

# 변이 호출만 실행  
snakemake --cores 96 results/PT001/annotation/PT001.TUMOR.somatic.vep.vcf.gz

# 구조변이 분석만 실행
snakemake --cores 96 results/PT001/sv/PT001.sv.annotsv_intogenCCG.tsv
```

### 결과 요약 생성
```bash
# 모든 환자의 변이 수 요약
for patient in $(cut -f1 samples.tsv | tail -n +2 | sort -u); do
    echo "=== $patient ==="
    if [ -f "results/$patient/annotation/$patient.TUMOR.somatic.vep.vcf.gz" ]; then
        echo "체세포 변이: $(zcat results/$patient/annotation/$patient.TUMOR.somatic.vep.vcf.gz | grep -v '^#' | wc -l)"
    fi
    if [ -f "results/$patient/sv/$patient.sv.annotsv_intogenCCG.tsv" ]; then
        echo "구조변이: $(tail -n +2 results/$patient/sv/$patient.sv.annotsv_intogenCCG.tsv | wc -l)"
    fi
done
```

### 백업 권장사항
```bash
# 중요한 결과만 백업
rsync -av results/*/annotation/ /backup/tnbc_variants/
rsync -av results/*/sv/*.tsv /backup/tnbc_sv/
rsync -av results/*/cnv/ /backup/tnbc_cnv/
```

## 📞 지원

문제가 발생하면 다음을 확인하세요:

1. **로그 파일**: 에러 메시지 확인
2. **리소스**: 메모리/디스크 공간 충분한지
3. **입력 파일**: BAM 파일 경로가 올바른지
4. **권한**: 파일 읽기/쓰기 권한 있는지

추가 도움이 필요하면 생성된 로그 파일과 함께 문의하세요. 