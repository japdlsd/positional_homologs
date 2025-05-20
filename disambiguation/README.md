## to run dummy data

```bash
snakemake -j16 --snakefile ../workflows/disambiguation.smk --configfile dummy-anopheles2.yml --printshellcmds
```

## to run anopheles data

```bash
snakemake -j48 --snakefile ../workflows/disambiguation.smk --configfile anopheles2.yml --printshellcmds --scheduler=greedy
```