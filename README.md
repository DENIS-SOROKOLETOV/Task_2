# **Интегративный функциональный геномный анализ** 

В данном репозитории представлены результаты, отчеты и програмные решения для лабораторной работы "Интегративный функциональный геномный анализ" по дисциплине "Биоинформатика", 4 вариант

## **Информация об организме**

*   Вид: *Fusarium graminearum* (*Gibberella zeae*) — патоген пшеницы и ячменя.
*   Размер генома:  ~36 Мб (гаплоидный).

## **Биологический контекст**
Регуляция посредством метилирования ДНК кластеров биосинтетических генов (BGC) и вирулентности в ответ на доступность питательных веществ. Исследование роли ДНК-метилтрансфераз (DNMT) в восприятии сигналов среды и вторичном метаболизме

## **Дизайн эксперимента**
Факторный план 2$\times$2 (WT/$\Delta$$\Delta$$\times$PN/NPN)

**Контрасты**

* Дикий тип (WT) vs. двойной мутант по DNMT ($\Delta$ FgDim-2/$\Delta$ FgRid)  (потеря функции ДНК-метилтрансферазы)
* Предпочтительная питательная среда (PN) vs. непредпочтительная питательная среда (NPN) 

**Тип образца**

Аксенические мицелиальные культуры
* Условия PN: рост в течение 24 часов
* Условия NPN: рост в течение 6 часов

## **Материалы и методы**


### **RNA-Seq**

**Набор данных**

* GEO: [GSE140030](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140030)
* BioProject: [PRJNA588066](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA588066)
* Референс генома: [Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.dna.toplevel.fa](https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-62/fasta/fungi_ascomycota3_collection/fusarium_graminearum_ph_1_gca_000240135/dna/)
* Аннотация: [Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.62.gff3](https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-62/gff3/fungi_ascomycota3_collection/fusarium_graminearum_ph_1_gca_000240135/)
* Silva

| Запуск | Среда | Тип |
| ------ | ----- | --- |
| SRR10408671 | PN | WT | 
| SRR10408672 | PN | WT |
| SRR10408673 | PN | WT |
| SRR10408674 | PN | $\Delta$$\Delta$ |
| SRR10408675 | PN | $\Delta$$\Delta$ |
| SRR10408676 | PN | $\Delta$$\Delta$ |
| SRR10408683 | NPN | WT |
| SRR10408684 | NPN | WT |
| SRR10408685 | NPN | WT |
| SRR10408686 | NPN | $\Delta$$\Delta$ |
| SRR10408687 | NPN | $\Delta$$\Delta$ |
| SRR10408688 | NPN | $\Delta$$\Delta$ |

**Инструменты**
* sra-toolkit (v. 3.4.1)
```bash
conda install bioconda::sra-tools
```
* HISAT2 (v. 2.2.2)
```bash
conda install bioconda::hisat2
```
* bowtie2 (v. 2.5.5)
```bash
conda install bioconda::bowtie2
```
* samtools (v. 1.23.1)
```bash
conda install bioconda::samtools
```
* fastp (v. 1.3.3)
```bash
conda install bioconda::fastp
```
* fastqc (v. 0.12.1)
```bash
conda install bioconda::fastqc
```
* subread (v. 2.1.1)
```bash
conda install bioconda::subread
```
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) 
