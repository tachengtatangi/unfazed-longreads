# Unfazed-longreads：phasing tool for *de novo* SNVs（third-generation sequencing）

Unfazed-longreads identifies the parent of origin for *de novo* mutations from the third-generation sequencing bam file, accepting input from a bed file of variant information. Unfazed-longreads is applicable to *de novo* SNVs. We sincerely thank Jonathan R Belyeu for his second-generation sequencing phasing software unfazed(https://github.com/jbelyeu/unfazed), which allowed us to generate unfazed-longreads, see (and cite) his paper at https://academic.oup.com/bioinformatics/article/37/24/4860/6306405 and our paper at [<u>unpublished</u>]( 'When will it be published?').

## How it works

Similar to the workflow of [unfazed SNV](https://github.com/jbelyeu/unfazed#extended-read-backed-phasing-snvindeldeldupinv) phasing, we first obtain the bam file with *de novo* mutations throughout all the reads in the bam file, and then look for SNPs on these reads as informative sites, and determine whether the *de novo* mutation originates from the parent or the mother based on all SNPs.

## How to use it

Because we are improving on Jonathan R Belyeu's unfazed, we use the original unfazed installation, via conda, which recommends a new environment and requires python version 3.5 and above, see (https://github.com/jbelyeu/unfazed) for details.

```shell
conda install -c bioconda unfazed
```

After conda is installed, execute

```shell
python setup.py
```

If you see the following output after running, the configuration of unfazed-longreads is successful

```
Congratulations!
You have configured successfully!
```

<details>
    <summary>The options for unfazed-longreads are：</summary>
    <pre><code>
        UNFAZED v1.0.2
usage: unfazed [-h] [-v] -d DNMS -s SITES -p PED [-b BAM_DIR]
               [--bam-pairs [BAM_PAIRS [BAM_PAIRS ...]]] [-t THREADS]
               [-o {vcf,bed}] [--include-ambiguous] [--verbose]
               [--outfile OUTFILE] [-r REFERENCE] -g {37,38,na}
               [--no-extended] [--multiread-proc-min MULTIREAD_PROC_MIN] [-q]
               [--min-gt-qual MIN_GT_QUAL] [--min-depth MIN_DEPTH]
               [--ab-homref AB_HOMREF] [--ab-homalt AB_HOMALT]
               [--ab-het AB_HET] [--evidence-min-ratio EVIDENCE_MIN_RATIO]
               [--search-dist SEARCH_DIST]
               [--insert-size-max-sample INSERT_SIZE_MAX_SAMPLE]
               [--min-map-qual MIN_MAP_QUAL] [--stdevs STDEVS]
               [--readlen READLEN] [--split-error-margin SPLIT_ERROR_MARGIN]
               [--max-reads MAX_READS]
        optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed version (1.0.2)
  -d DNMS, --dnms DNMS  valid VCF OR BED file of the DNMs of interest> If BED,
                        must contain chrom, start, end, kid_id, var_type
                        columns (default: None)
  -s SITES, --sites SITES
                        sorted/bgzipped/indexed VCF/BCF file of SNVs to
                        identify informative sites. Must contain each kid and
                        both parents (default: None)
  -p PED, --ped PED     ped file including the kid and both parent IDs
                        (default: None)
  -b BAM_DIR, --bam-dir BAM_DIR
                        directory where bam/cram files (named {sample_id}.bam
                        or {sample_id}.cram) are stored for offspring. If not
                        included, --bam-pairs must be set (default: None)
  --bam-pairs [BAM_PAIRS [BAM_PAIRS ...]]
                        space-delimited list of pairs in the format
                        {sample_id}:{bam_path} where {sample_id} matches an
                        offspring id from the dnm file. Can be used with
                        --bam-dir arg, must be used in its absence (default:
                        None)
  -t THREADS, --threads THREADS
                        number of threads to use (default: 2)
  -o {bed}, --output-type {bed}
                        choose output type. If --dnms is not a VCF/BCF, output
                        must be to BED format. Defaults to match --dnms input
                        file (default: None)
  --include-ambiguous   include ambiguous phasing results (default: False)
  --verbose             print verbose output including sites and reads used
                        for phasing. Only applies to BED output (default:
                        False)
  --outfile OUTFILE     name for output file. Defaults to stdout (default:
                        /dev/stdout)
  -r REFERENCE, --reference REFERENCE
                        reference fasta file (required for crams) (default:
                        None)
  -g {37,38,na}, --build {37,38,na}
                        human genome build, used to determine sex chromosome
                        pseudoautosomal regions. If `na` option is chosen, sex
                        chromosomes will not be auto-phased. HG19/GRCh37
                        interchangeable (default: None)
  --no-extended         do not perform extended read-based phasing (default
                        True) (default: False)
  --multiread-proc-min MULTIREAD_PROC_MIN
                        min number of variants required to perform multiple
                        parallel reads of the sites file (default: 1000)
  -q, --quiet           no logging of variant processing data (default: False)
  --min-gt-qual MIN_GT_QUAL
                        min genotype and base quality for informative sites
                        (default: 20)
  --min-depth MIN_DEPTH
                        min coverage for informative sites (default: 10)
  --ab-homref AB_HOMREF
                        allele balance range for homozygous reference
                        informative sites (default: 0.0:0.2)
  --ab-homalt AB_HOMALT
                        allele balance range for homozygous alternate
                        informative sites (default: 0.8:1.0)
  --ab-het AB_HET       allele balance range for heterozygous informative
                        sites (default: 0.2:0.8)
  --evidence-min-ratio EVIDENCE_MIN_RATIO
                        minimum ratio of evidence for a parent to provide an
                        unambiguous call. Default 10:1 (default: 1)
  --search-dist SEARCH_DIST
                        maximum search distance from variant for informative
                        sites (in bases) (default: 30000)
  --insert-size-max-sample INSERT_SIZE_MAX_SAMPLE
                        maximum number of read inserts to sample in order to
                        estimate concordant read insert size (default:
                        1000000)
  --min-map-qual MIN_MAP_QUAL
                        minimum map quality for reads (default: 1)
  --stdevs STDEVS       number of standard deviations from the mean insert
                        length to define a discordant read (default: 3)
  --readlen READLEN     expected length of input reads (default: 15000)
  --split-error-margin SPLIT_ERROR_MARGIN
                        margin of error for the location of split read
                        clipping in bases (default: 5)
  --max-reads MAX_READS
                        maximum number of reads to collect for phasing a
                        single variant (default: 100)
                        </code></pre>
</details>

### For the third-generations sequencing of phasing, our suggested usage is:

```shell
unfazed\
    -d alltri.bed\
    -s 729.vcf.gz\
    -p 1.ped\
    -o bed\
    --build na\
    --verbose\
    --bam-pairs '729':729tiger.sort.bam\
    -r Tiger_VersionChr.genome.fa\
    --outfile result.bed
```

The contents of the bed, vcf.gz, and ped files are exemplified below, and `'729'` in `-bam-pairs '729':7222tiger.sort.bam` is the child name.

## Input and output of unfazed-longreads

### bed input

The input bed file must have the following tab-separated columns: chrom, start, end, kid_id, var_type, where vartype is SNV.

```
#chrom  start   end     kid_id  var_type
A1      3843125 3843126 729     SNV
A1      80958650        80958651        729     SNV
A1      190290922       190290923       729     SNV
A1      232566922       232566923       729     SNV
B1      137225256       137225257       729     SNV
B1      190227222       190227223       729     SNV
A2      2558972 2558973 729     SNV
A2      18860243        18860244        729     SNV
A2      58547780        58547781        729     SNV
```

### vcf input

For detailed information see [Jonathan R Belyeu](https://github.com/jbelyeu/unfazed#vcf-annotations 'vcf')

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  yangyang        huihui  729
A1      508     .       C       A       293.25  PASS    AC=0;AF=0.333;AN=2;BaseQRankSum=2.22;DP=71;ExcessHet=4.7712;FS=5.092;MLEAC=3;MLEAF=0.5;MQ=41.25;MQRankSum=0.674;QD=8.15;ReadPosRankSum=3.65;SOR=1.776   GT:AD:DP:GQ:PL  ./.:2,8:10:.:0,0,0      ./.:0,7:7:.:0,0,0       0/0:6,3:9:0:0,0,80
A1      650     .       T       C       37.65   PASS    AC=0;AF=0.1;AN=6;BaseQRankSum=0;DP=56;ExcessHet=3.0103;FS=3.09;MLEAC=1;MLEAF=0.1;MQ=56.43;MQRankSum=-2.838;QD=3.42;ReadPosRankSum=-0.068;SOR=2.494      GT:AD:DP:GQ:PL  0/0:9,0:9:21:0,21,315   0/0:16,1:17:5:0,5,543   0/0:7,2:9:0:0,0,165
A1      697     .       C       A       51.3    PASS    AC=1;AF=0.1;AN=6;BaseQRankSum=0.696;DP=61;ExcessHet=3.0103;FS=1.705;MLEAC=1;MLEAF=0.1;MQ=57.82;MQRankSum=-1.521;QD=5.13;ReadPosRankSum=1.07;SOR=0.321   GT:AD:DP:GQ:PGT:PID:PL:PS       0|0:11,0:14:33:1|0:686_T_TAAAAATAGAGATACCCTATGATCCGG:0,33,495:686       0/0:16,1:17:14:.:.:0,14,624:.   0|1:8,2:13:60:0|1:697_C_A:60,0,330:697
```

### ped input

The input ped file must have the following tab-separated columns: Family-ID, Individual-ID, Paternal-ID, Maternal-ID, Gender, where Family-ID can take an arbitrary number to indicate the pedigree, and Gender where 1 represents males and 2 represents females.

```
#Family-ID	Individual-ID	Paternal-ID	Maternal-ID	Gender
1463	729	yangyang	huihui	1
1463	yangyang	0	0	1
1463	huihui	0	0	2
```

### bed output

```
#chrom	start	end	vartype	kid	origin_parent	other_parent	evidence_count	evidence_types	origin_parent_sites	origin_parent_reads	other_parent_sites	other_parent_reads
A1	3843125	3843126	POINT	729	yangyang	huihui	177	READBACKED	3784986...	m64067_220520_105315/99288489/37898_60841_np12...	3784986...	m64067_220520_105315/166135198/103529_127402_np12...
A1	232566922	232566923	POINT	729	yangyang	huihui	138	READBACKED	232508189...	m64067_220520_105315/54199284/38561_40508_np12...	232508189...	m64067_220520_105315/141493925/145329_167749_np12...
A2	2558972	2558973	POINT	729	yangyang	huihui	18	READBACKED	2505731...	m64067_220520_105315/131204883/139696_161199_np12...	2505731...	m64067_220520_105315/179964199/18076_37423_np12...
A2	18860243	18860244	POINT	729	yangyang	huihui	17	READBACKED	18811021...	m64067_220520_105315/108333713/87912_105679_np12...	18849838...	m64067_220520_105315/70255600/251052_253839_np12...
```

Due to sequencing errors (e.g. CLR sequencing) or other potential reasons, the program may not be accurate for some sites, and it is almost impossible to manually check for these ambiguous results, so we eliminated the ambiguous results.

## Note

- Make sure your vcf file is consistent with the names of the offspring and their parents in the ped file or an error will be reported!
- Please make sure your bed file, bam file, vcf file and the chrom names in the reference genome are consistent, or an error will be reported!

## Performance

The phasing  results achieve a successful rate higher than 60%, which improved one fold compared with the old version (a little under than 30%) using our unfazed-longreads basing on CLR bams (or other long-read sequencing bam files).  Consistent with [Jonathan R Belyeu](https://github.com/jbelyeu/unfazed 'unfazed'), unfazed-longreads cannot achieve 100% phasing rate  because of  lacking valid SNP sites for some DNMs.
