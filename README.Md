# Biotagging

Various scripts used in the Biotagging project carried out in the laboratory of T. Sauka-Spengler at the Weatherall Institute of Molecular Medicine at the University of Oxford, in the UK.

#### To reproduce analysis in manuscript, required software: 

* bowtie v.1.1.2
* bedtools v.2.15.0
* Kent tools  bedGraphToBigWig (downloadable from UCSC genome browser website -[here]())
* python v.2.7.5

#### Files :

* bowtie v.1.1.2 index for danRer7 zebrafish genome downloadable from [here](http://userweb.molbiol.ox.ac.uk:/public/dariag/damir/danRer7_bowtie1_indices.tar.gz)
* chromosome size file : danRer7.sizes

#### Example commands for ATAC-seq parsing: 

```{r eval = FALSE}
bowtie -S -p 4 -X 2000 -m 2 $GENOME -1 $A1 -2 $A2  --chunkmb 500 $AN\.sam

samtools view -bS $AN\.sam > $AN\.bam
samtools sort $AN\.bam $AN\.sort
genomeCoverageBed -bg -split -ibam $AN\.sort.bam -g $CHROM > $AN\.bg
bedGraphToBigWig $AN\.bg $CHROM $AN\.bw

sam2bwPE_zf.pl -sam $AN\.sam -build danRer10 -name smo_$AN

samtools sort -n $AN\.bam $AN\.nsort
bedtools bamtobed -bedpe -i $AN\.nsort.bam > $AN\.nsort.bam.bed
atac_bedpe_parse2.py $CHROM $AN\.nsort.bam.bed
macs2 callpeak -t $AN\.nsort.bam.pebed -f BED --name macs2_$AN --shiftsize=100 --nomodel --slocal 1000
```


Original sam2bwPE_PE.pl script was kindly provided by [Pr. Jim Hughes](https://github.com/Hughes-Genome-Group) (WIMM, University of Oxford). Version presented here has been modified to work with zebrafish genome. 


#### To split RNA-seq data by strand use split-strand-rna.py python script: 

```{r eval = FALSE}
samtools sort -n myFile.bam myFile.nsort
samtools view -h myFile.nsort.bam | python split-strand-rna.py PreFixToOutputFile
```


This will yield two files : PreFixToOutputFile_+.sam and PreFixToOutputFile_-.sam
To visualise these files can be processed with samtools, bedtools genomeCoverageBed and bedGraphToBigWig. Bigwig files can then be uploaded and visualised in UCSC genome browser. 

#### Hub hosting all data : 

http://userweb.molbiol.ox.ac.uk/public/dariag/biotagging/hub.txt


#### Additional tracks of potentional intrest: 

track type=bigBed name=" Zv9_introns"  bigDataUrl=http://userweb.molbiol.ox.ac.uk/public/dariag/biotag_ss/Zv9_introns.bb


For additional information consult github hosted by [tsslab/biotagging](https://github.com/tsslab/biotagging)

#### Manuscript information 


Trinh, Le A., Vanessa Chong-Morrison, Daria Gavriouchkina, Tatiana Hochgreb-Hägele, Upeka Senanayake, Scott E. Fraser, and Tatjana Sauka-Spengler. "Biotagging of specific cell populations in zebrafish reveals gene regulatory logic encoded in the nuclear transcriptome." Cell reports 19, no. 2 (2017): 425-440.
[doi:10.1016/j.celrep.2017.03.045](https://www.sciencedirect.com/science/article/pii/S2211124717303911)
