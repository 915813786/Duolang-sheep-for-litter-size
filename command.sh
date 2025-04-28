#quality evalutation before quality control
fastqc -t 4 -o /home/chao/sheep/litterSize/00.qualityControl/D96 /home/chao/sheep/litterSize/D96_1.fq.gz /home/chao/sheep/litterSize/D96_2.fq.gz
seqkit stat /home/chao/sheep/litterSize/D96_1.fq.gz > /home/chao/sheep/litterSize/00.qualityControl/D96.stat
seqkit stat /home/chao/sheep/litterSize/D96_2.fq.gz >> /home/chao/sheep/litterSize/00.qualityControl/D96.stat

#quality control using Trimmomatic v0.39
java -jar /home/chao/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 /home/chao/sheep/litterSize/D96_1.fq.gz \
    /home/chao/sheep/litterSize/D96_2.fq.gz \
	/home/chao/sheep/litterSize/00.qualityControl/D96_1.clean.fq.gz \
	/home/chao/sheep/litterSize/00.qualityControl/D96_1.unpaired.fq.gz \
	/home/chao/sheep/litterSize/00.qualityControl/D96_2.clean.fq.gz \
	/home/chao/sheep/litterSize/00.qualityControl/D96_2.unpaired.fq.gz \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90

#quality evalutation after quality control
fastqc -t 4 -o /home/chao/sheep/litterSize/00.qualityControl/D96 /home/chao/sheep/litterSize/00.qualityControl/D96_1.clean.fq.gz /home/chao/sheep/litterSize/00.qualityControl/D96_2.clean.fq.gz
seqkit stat /home/chao/sheep/litterSize/00.qualityControl/D96_1.clean.fq.gz >> /home/chao/sheep/litterSize/00.qualityControl/D96.stat
seqkit stat /home/chao/sheep/litterSize/00.qualityControl/D96_2.clean.fq.gz >> /home/chao/sheep/litterSize/00.qualityControl/D96.stat

#unzip the output of fastqc
unzip /home/chao/sheep/litterSize/00.qualityControl/D1/D1_1_fastqc.zip
unzip /home/chao/sheep/litterSize/00.qualityControl/D1/D1_2_fastqc.zip
unzip /home/chao/sheep/litterSize/00.qualityControl/D1/D1_1.clean_fastqc.zip
unzip /home/chao/sheep/litterSize/00.qualityControl/D1/D1_2.clean_fastqc.zip

#perl script for extract key information from the results of qualimap bamqc
perl /home/chao/sheep/litterSize/00.qualityControl/fastq.stat.pl /home/chao/sheep/litterSize/00.qualityControl/list /home/chao/sheep/litterSize/00.qualityControl/fastq.stat #list is a file with absolute directory of the output of seqkit and fastqc for all individuals (e.g. /home/chao/sheep/litterSize/00.qualityControl/D96.stat	/home/chao/sheep/litterSize/00.qualityControl/D96/D96_1_fastqc/fastqc_data.txt	/home/chao/sheep/litterSize/00.qualityControl/D96/D96_2_fastqc/fastqc_data.txt	/home/chao/sheep/litterSize/00.qualityControl/D96/D96_1.clean_fastqc/fastqc_data.txt	/home/chao/sheep/litterSize/00.qualityControl/D96/D96_2.clean_fastqc/fastqc_data.txt)

####################################R code of density plot for bam statistics
library(ggplot2)
a<-read.table("fastq.stat",header=F,sep="\t")
p1<-ggplot(a, aes(x = a[,3])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,3])
p1<-p1+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p1<-p1+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p2<-ggplot(a, aes(x = a[,4])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,4])
p2<-p2+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p2<-p2+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p3<-ggplot(a, aes(x = a[,5])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,5])
p3<-p3+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p3<-p3+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p4<-ggplot(a, aes(x = a[,6])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,6])
p4<-p4+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p4<-p4+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p5<-ggplot(a, aes(x = a[,7])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,7])
p5<-p5+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p5<-p5+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
library(gridExtra)
grid.arrange(p1,p2,p3,p4,p5, nrow = 2, ncol = 3)
#######################################

#alignment using bwa v0.7.17-r1188, and format transformation from sam to bam and sort by chromosome coordinate using samtools v1.17
bwa mem -t 4 -R '@RG\\tID:D96\\tPL:illumina\\tLB:D96\\tSM:D96' /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    /home/chao/sheep/litterSize/00.qualityControl/D96_1.clean.fq.gz \
	/home/chao/sheep/litterSize/00.qualityControl/D96_2.clean.fq.gz | samtools view -@ 4 -S -b - | samtools sort -@ 4 -T /home/chao/tmp -O bam -o /home/chao/sheep/litterSize/01.bam/D96.sort.bam -

#mark duplicate reads using gatk v4.4.0.0
gatk MarkDuplicates -I /home/chao/sheep/litterSize/01.bam/D96.sort.bam \
    -O /home/chao/sheep/litterSize/01.bam/D96.mark.bam -M /home/chao/sheep/litterSize/01.bam/D96.metrics.txt

#synchronize mate-pair information using gatk 	
gatk FixMateInformation -I /home/chao/sheep/litterSize/01.bam/D96.mark.bam -O /home/chao/sheep/litterSize/01.bam/D96.marked_fixed.bam -SO coordinate

#create bam index file using samtools
samtools index /home/chao/sheep/litterSize/01.bam/D96.marked_fixed.bam

#generate statistic information about coverge, mapping quality, etc using qualimap v2.2.1
qualimap bamqc -bam /home/chao/sheep/litterSize/01.bam/D96.marked_fixed.bam -outdir /home/chao/sheep/litterSize/01.bam/D96 -outformat HTML --java-mem-size=100G

#perl script for extract key information from the results of qualimap bamqc
perl /home/chao/sheep/litterSize/01.bam/bam.stat.pl /home/chao/sheep/litterSize/01.bam/list /home/chao/sheep/litterSize/01.bam/bam.stat #list is a file with absolute directory of genome_results.txt for all individuals (e.g. /home/chao/sheep/litterSize/01.bam/D96/genome_results.txt)

####################################R code of density plot for bam statistics
library(ggplot2)
a<-read.table("bam.stat",header=F,sep="\t")
a[,4] <- as.numeric(gsub("%", "", a[,4])) / 100
p4<-ggplot(a, aes(x = a[,4])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,4])
p4<-p4+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p4<-p4+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p5<-ggplot(a, aes(x = a[,7])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,7])
p5<-p5+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p5<-p5+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
p6<-ggplot(a, aes(x = a[,10])) + geom_density(fill = "blue",adjust = 10)+theme_minimal()
median_value <- median(a[,10])
p6<-p6+geom_vline(xintercept = median_value, linetype = "dashed", color = "red")
p6<-p6+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
library(gridExtra)
grid.arrange(p4, p5, p6, nrow = 1, ncol = 3)
#######################################

#call variants (SNPs and indels) using gatk
gatk HaplotypeCaller --native-pair-hmm-threads 10 -ERC GVCF -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    -I /home/chao/sheep/litterSize/01.bam/D96.marked_fixed.bam -L 26 -O /home/chao/sheep/litterSize/02.gvcf/D96/26.g.vcf.gz
	
#merge multiple GVCF (Genomic VCF) files using gatk
java -Xmx200g -jar /home/chao/software/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CombineGVCFs -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
	-V /home/chao/sheep/litterSize/02.gvcf/D10/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D100/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D101/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D102/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D103/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D104/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D105/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D106/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D107/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D108/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D109/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D110/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D111/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D112/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D113/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D114/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D115/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D116/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D117/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D118/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D119/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D12/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D120/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D121/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D122/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D123/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D124/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D125/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D126/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D127/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D128/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D129/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D13/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D130/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D131/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D132/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D133/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D134/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D135/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D136/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D137/26.g.vcf.gz \
	-V /home/chao/sheep/litterSize/02.gvcf/D138/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D139/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D14/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D140/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D141/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D142/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D143/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D144/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D145/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D146/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D147/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D148/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D149/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D15/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D150/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D151/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D152/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D153/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D154/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D155/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D156/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D157/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D158/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D159/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D16/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D17/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D18/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D19/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D2/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D20/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D21/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D22/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D23/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D24/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D25/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D27/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D28/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D29/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D3/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D30/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D31/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D32/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D33/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D34/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D35/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D36/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D37/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D38/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D39/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D4/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D40/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D96/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D42/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D43/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D44/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D45/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D46/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D47/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D48/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D49/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D5/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D50/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D51/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D52/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D53/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D54/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D55/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D56/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D57/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D58/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D59/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D6/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D60/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D61/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D62/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D63/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D64/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D65/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D66/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D67/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D68/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D69/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D7/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D70/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D71/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D72/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D74/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D75/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D76/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D77/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D78/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D79/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D8/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D80/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D81/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D82/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D83/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D84/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D85/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D86/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D88/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D89/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D9/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D90/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D92/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D93/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D94/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D95/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/D96/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S1/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S10/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S100/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S101/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S102/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S103/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S104/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S105/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S106/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S107/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S108/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S109/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S11/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S110/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S111/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S112/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S115/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S116/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S117/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S118/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S119/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S12/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S120/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S121/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S122/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S123/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S124/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S125/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S126/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S127/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S128/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S129/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S13/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S131/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S132/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S133/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S134/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S135/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S136/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S138/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S139/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S14/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S140/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S141/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S142/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S143/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S144/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S145/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S146/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S147/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S148/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S149/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S150/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S151/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S152/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S153/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S154/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S16/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S17/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S18/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S2/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S20/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S21/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S22/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S23/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S24/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S25/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S26/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S27/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S28/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S29/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S3/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S30/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S31/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S32/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S33/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S34/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S35/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S36/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S37/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S38/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S39/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S4/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S40/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S41/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S42/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S43/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S44/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S45/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S46/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S47/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S48/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S49/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S5/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S50/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S51/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S52/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S53/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S54/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S55/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S56/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S57/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S58/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S59/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S6/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S60/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S61/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S62/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S63/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S64/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S65/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S66/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S67/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S68/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S69/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S7/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S70/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S71/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S72/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S73/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S74/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S75/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S76/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S77/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S78/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S79/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S8/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S80/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S81/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S82/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S83/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S84/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S85/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S86/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S87/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S88/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S9/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S90/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S91/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S92/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S93/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S94/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S95/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S96/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S97/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S98/26.g.vcf.gz \
    -V /home/chao/sheep/litterSize/02.gvcf/S99/26.g.vcf.gz \
	-O /home/chao/sheep/litterSize/03.mergeGvcf/26.g.vcf.gz

# joint genotyping using gatk
gatk GenotypeGVCFs -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/03.mergeGvcf/1.g.vcf.gz -O /home/chao/sheep/litterSize/04.snp/1.vcf.gz

#Select SNP type
gatk SelectVariants -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/04.snp/1.vcf.gz -O /home/chao/sheep/litterSize/04.snp/1.raw.SNP.vcf.gz --select-type-to-include SNP

#hard filtering for SNP
gatk VariantFiltration -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/04.snp/1.raw.SNP.vcf.gz \
	-O /home/chao/sheep/litterSize/04.snp/1.filter.SNP.vcf.gz \
	--filter-name "snp_filter" --filter-expression "QD < 3.0 || FS > 30.0 || SOR > 4.0 || MQ < 30.0 || MQRankSum < -10.0 || QUAL < 50.0 || ReadPosRankSum < -5.0"

#exclude any variants that have failed filtering criteria, and select biallelic SNPs using gatk
gatk SelectVariants -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/04.snp/1.filter.SNP.vcf.gz \
	-O /home/chao/sheep/litterSize/04.snp/1.clean.SNP.vcf.gz \
	--select-type-to-include SNP --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC -select "AC > 0" -select "AF < 1.00"
	
#select INDEL type
gatk SelectVariants -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/04.snp/1.vcf.gz \
	-O /home/chao/sheep/litterSize/04.snp/1.raw.indel.vcf.gz --select-type-to-include INDEL

#hard filtering for INDEL
gatk VariantFiltration -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa \
    --variant /home/chao/sheep/litterSize/04.snp/1.raw.indel.vcf.gz \
	-O /home/chao/sheep/litterSize/04.snp/1.filter.indel.vcf.gz \
	--filter-name "INDEL_filter" --filter-expression "QD < 3.0 || FS > 100.0 || ReadPosRankSum < -10.0"

#exclude any variants that have failed filtering criteria, and select biallelic INDELs using gatk
gatk SelectVariants -R /home/chao/software/annovar/ovinedb/ARS-UI.v3.0.fa --variant /home/chao/sheep/litterSize/04.snp/1.filter.indel.vcf.gz -O /home/chao/sheep/litterSize/04.snp/1.clean.indel.vcf.gz --select-type-to-include INDEL --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC -select "AC > 0" -select "AF < 1.00"

#merge SNP and INDELs for each chromosome
/home/chao/software/gatk-4.4.0.0/gatk MergeVcfs -I /home/chao/sheep/litterSize/04.snp/26.clean.indel.vcf.gz \
    -I /home/chao/sheep/litterSize/04.snp/26.clean.SNP.vcf.gz -O /home/chao/sheep/litterSize/04.snp/26.clean.vcf.gz
	
#merge variants for all chromosomes
/home/chao/software/gatk-4.4.0.0/gatk MergeVcfs -I /home/chao/sheep/litterSize/04.snp/1.clean.vcf.gz \
    -I /home/chao/sheep/litterSize/04.snp/2.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/3.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/4.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/5.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/6.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/7.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/8.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/9.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/10.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/11.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/12.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/13.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/14.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/15.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/16.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/17.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/18.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/19.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/20.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/21.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/22.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/23.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/24.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/25.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/26.clean.vcf.gz \
	-I /home/chao/sheep/litterSize/04.snp/X.clean.vcf.gz \
	-O /home/chao/sheep/litterSize/04.snp/All.clean.vcf.gz

##Further filtering
perl /home/chao/sheep/litterSize/04.snp/filter.pl /home/chao/sheep/litterSize/04.snp/bam.stat /home/chao/sheep/litterSize/04.snp/All.clean.vcf.gz /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz


#gene annotation using annovar v2016-02-01
perl /home/chao/software/annovar/table_annovar.pl /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz /home/chao/software/annovar/ovinedb/ -out /home/chao/sheep/litterSize/04.snp/All -buildver v3UseName -protocol refGene -operation g -nastring . -vcfinput

#the distribution of variants across different regions of the genome, such as exonic, intronic, UTR, upstream, and downstream regions.
awk '{print $1}' All.refGene.variant_function  |sort | uniq -c
##########################################R code for plot the functional classification of the detected variants.
data <- data.frame(
  value = c(424614, 379213, 118, 27454293, 18163662, 703043, 95, 1857774, 774, 6, 3196, 388227, 15477, 635855, 440674, 5154),
  category = c("downstream", "exonic", "exonic;splicing", "intergenic", "intronic", "ncRNA_exonic", "ncRNA_exonic;splicing",
               "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5", "splicing", "upstream", "upstream;downstream", "UTR3", "UTR5", "UTR5;UTR3")
)
data$merged_category <- with(data, ifelse(category %in% c("UTR3", "UTR5", "UTR5;UTR3"), "UTR",
                                ifelse(category %in% c("downstream", "upstream", "upstream;downstream"), "upstream and downstream",
                                ifelse(category %in% c("exonic", "exonic;splicing", "splicing", "ncRNA_splicing", "ncRNA_exonic", "ncRNA_exonic;splicing"), "exonic",
                                ifelse(category %in% c("intronic", "ncRNA_intronic"), "intronic", "intergenic")))))
merged_data <- aggregate(value ~ merged_category, data = data, sum)
gap <- 0.08  
merged_data$fraction <- merged_data$value / sum(merged_data$value)
merged_data$fraction_adjusted <- merged_data$fraction * (1 - gap)  
merged_data$ymax <- cumsum(merged_data$fraction_adjusted) + gap * (1:nrow(merged_data)) / nrow(merged_data)
merged_data$ymin <- c(0, head(merged_data$ymax, -1))

library(ggplot2)
ggplot(merged_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = merged_category)) +
  geom_rect(color = "white") +  
  coord_polar(theta = "y") +
  labs(fill = "Category", title = "Pie Chart of Genomic Categories (With Gaps)") +
  theme_void() +
  theme(legend.position = "right")
##########################################################################

#the specific types of variants within the exonic region, such as synonymous variants, nonsynonymous variants, frameshift insertions/deletions, and stop-gain or stop-loss mutations.
perl -ne '{chomp;@F=split/\t/,$_;print "$F[1]\n"}' All.refGene.exonic_variant_function | sort | uniq -c
#######################################################################R code for plot the functional classification of the detected variants in exonic region.
data <- data.frame(
  value = c(32578, 3295, 5080, 1541, 145815, 3402, 343, 174533, 12744),
  category = c("frameshift deletion", "frameshift insertion", "nonframeshift deletion", 
               "nonframeshift insertion", "nonsynonymous SNV", "stopgain", 
               "stoploss", "synonymous SNV", "unknown")
)
data$fraction <- data$value / sum(data$value)
gap <- 0.08  
data$ymax <- cumsum(data$fraction) + gap * (1:nrow(data)) / nrow(data)
data$ymin <- c(0, head(data$ymax, n = -1))
library(ggplot2)
ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = category)) +
  geom_rect(color = "white") +  
  coord_polar(theta = "y") +
  labs(fill = "Category", title = "Pie Chart of Mutation Types (With Gaps)") +
  theme_void() +
  theme(legend.position = "right")
#######################################################################

#the Transition / Transversion ratio using vcftools v0.1.17
vcftools --gzvcf /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz --TsTv-summary --out /home/chao/sheep/litterSize/05.tstv/tstv

#the total number of SNPs and INDELs
vcftools --gzvcf /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz --remove-indels --recode --stdout | grep -v "^#"  |wc -l
vcftools --gzvcf /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz --keep-only-indels --recode --stdout | grep -v "^#"  |wc -l

#plot variants density using CMplot v4.5.1
zcat /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz | grep -v "^#" | perl -lane 'print "$F[0]:$F[1]\t$F[0]\t$F[1]\t0.1"' | sed '1i SNP\tChromosome\tPosition\tP.value' > /home/chao/sheep/litterSize/04.snp/07.densityPlot/density
############################################R code for plotting variants density
data <- read.table("density", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library(CMplot)
 CMplot(data,plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="tif",file.name="",dpi=300,
    main="Variants density",file.output=TRUE,verbose=TRUE,width=9,height=6)
#################################################
	
#additive effect: Generate MAP, PED, BIM, FAM, and BED files using VCFtools and PLINK v1.90b6.21, applying filters with a minor allele frequency (MAF) threshold of 0.05, allowing a maximum of 10% missing data per site, and selecting sites with mean depth between one-third and three times the overall average sequencing coverage
vcftools --gzvcf /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz --max-missing 0.9 --maf 0.05  --plink --out /home/chao/sheep/litterSize/08.plink/All
plink --file /home/chao/sheep/litterSize/08.plink/All --chr-set 26 --make-bed --out /home/chao/sheep/litterSize/08.plink/All --snps-only

#prepare phenotype file
perl -ne '{chomp;@F=split/\s+/,$_;if (/^D/){print "$F[0] $F[0] 1\n";}else{print "$F[0] $F[0] 0\n";}}' /home/chao/sheep/litterSize/08.plink/All.fam > /home/chao/sheep/litterSize/10.gcta/phenotype.txt

#gcta v1.94.1
gcta --bfile /home/chao/sheep/litterSize/08.plink/All --autosome-num 26 --make-grm  --out /home/chao/sheep/litterSize/10.gcta/All
gcta --grm /home/chao/sheep/litterSize/10.gcta/All --make-bK-sparse 0.05 --out /home/chao/sheep/litterSize/10.gcta/sp_grm --autosome-num 26
gcta --grm /home/chao/sheep/litterSize/10.gcta/All --pca 5 --out /home/chao/sheep/litterSize/10.gcta/All
gcta --bfile /home/chao/sheep/litterSize/08.plink/All --autosome-num 26 --grm-sparse /home/chao/sheep/litterSize/10.gcta/sp_grm --joint-covar --fastGWA-mlm-binary --pheno /home/chao/sheep/litterSize/10.gcta/phenotype.txt --qcovar /home/chao/sheep/litterSize/10.gcta/All.eigenvec --thread-num 10 --out /home/chao/sheep/litterSize/10.gcta/assoc

#plot principal component analysis
perl -lane 'if (/^S/){print "$F[0]\tSingle\t$F[3]\t$F[4]";}else{print "$F[0]\tDouble\t$F[3]\t$F[4]";}' /home/chao/sheep/litterSize/10.gcta/All.eigenvec > /home/chao/sheep/litterSize/10.gcta/pca.txt
######################################################R code for principal component analysis
a<-read.table("E:/research/sheep/litterSize/pca.txt",header=F,sep="\t")
library(ggplot2)
p<-ggplot(a,aes(a[,3],a[,4],col=a[,2]))+geom_point(cex=2)+theme_set(theme_bw())+theme(panel.grid.minor=element_blank())+theme(panel.grid.major=element_blank())+theme(panel.background=element_blank())
######################################################

##manhattan plot and QQ plot
awk '{print $2,$1,$3,$13}' OFS="\t" /home/chao/sheep/litterSize/10.gcta/assoc.fastGWA | sed '1d' | sed '1i SNP\tChromosome\tPosition\tadd' > /home/chao/sheep/litterSize/10.gcta/gwas.plot
###########################################R code for manhattan plot, QQ plot, and Inflation factor (lambda)
data <- read.table("gwas.plot", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data$Chromosome <- ifelse(data$Chromosome == 27, "X", data$Chromosome)
data <- data[!(data$add == 0 | data$add == 1), ]
library("CMplot")
CMplot(data, plot.type="m", col=c("grey60"), LOG10=TRUE,threshold=c(1e-6),
        threshold.lty=c(1), threshold.lwd=c(1), threshold.col=c("black"), amplify=TRUE,
        chr.den.col=NULL, signal.col=c("red"), signal.cex=c(1.5),signal.pch=c(19),
        file="tiff",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
CMplot(data,plot.type="q",box=FALSE,file="tif",file.name="",dpi=300,
    conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
    file.output=TRUE,verbose=TRUE,width=5,height=5)
#calculate Inflation factor (lambda) 
data$Chi2 <- qchisq(1 - data$add, df = 1)
median_chi2_obs <- median(data$Chi2, na.rm = TRUE)
lambda <- median_chi2_obs / 0.455
cat("Inflation factor (lambda):", lambda, "\n")
#####################################################

#dominant effect:Generate MAP, PED, BIM, FAM, and BED files 
vcftools --gzvcf /home/chao/sheep/litterSize/04.snp/All.filtered.clean.vcf.gz --max-missing 0.9 --maf 0.05 --recode --out /home/chao/sheep/litterSize/09.plinkDominant/All
gzip /home/chao/sheep/litterSize/09.plinkDominant/All.recode.vcf
perl /home/chao/sheep/litterSize/09.plinkDominant/dominantGenotype.pl /home/chao/sheep/litterSize/09.plinkDominant/All.recode.vcf.gz /home/chao/sheep/litterSize/09.plinkDominant/All.dominant.vcf.gz
vcftools --gzvcf /home/chao/sheep/litterSize/09.plinkDominant/All.dominant.vcf.gz --plink --out /home/chao/sheep/litterSize/09.plinkDominant/All
plink --file /home/chao/sheep/litterSize/09.plinkDominant/All --chr-set 26 --make-bed --out /home/chao/sheep/litterSize/09.plinkDominant/All  --snps-only

#gcta v1.94.1
gcta --bfile /home/chao/sheep/litterSize/09.plinkDominant/All --autosome-num 26 --make-grm  --out /home/chao/sheep/litterSize/11.gctaDominant/All
gcta --grm /home/chao/sheep/litterSize/11.gctaDominant/All --make-bK-sparse 0.05 --out /home/chao/sheep/litterSize/11.gctaDominant/sp_grm --autosome-num 26
gcta  --grm /home/chao/sheep/litterSize/11.gctaDominant/All --pca 5 --out /home/chao/sheep/litterSize/11.gctaDominant/All
gcta --bfile /home/chao/sheep/litterSize/09.plinkDominant/All --autosome-num 26 --grm-sparse /home/chao/sheep/litterSize/11.gctaDominant/sp_grm --joint-covar --fastGWA-mlm-binary --pheno /home/chao/sheep/litterSize/10.gcta/phenotype.txt --qcovar /home/chao/sheep/litterSize/11.gctaDominant/All.eigenvec --thread-num 10 --out /home/chao/sheep/litterSize/11.gctaDominant/assoc

##manhattan plot and QQ plot
awk '{print $2,$1,$3,$13}' OFS="\t" /home/chao/sheep/litterSize/11.gctaDominant/assoc.fastGWA | sed '1d' | sed '1i SNP\tChromosome\tPosition\tdom' > /home/chao/sheep/litterSize/11.gctaDominant/gwas.plot
###########################################R code for manhattan plot, QQ plot, and Inflation factor (lambda)
data <- read.table("gwas.plot", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data$Chromosome <- ifelse(data$Chromosome == 27, "X", data$Chromosome)
data <- data[!(data$dom == 0 | data$dom == 1), ]
library("CMplot")
CMplot(data, plot.type="m", col=c("grey60"), LOG10=TRUE,threshold=c(1e-6),
        threshold.lty=c(1), threshold.lwd=c(1), threshold.col=c("black"), amplify=TRUE,
        chr.den.col=NULL, signal.col=c("red"), signal.cex=c(1.5),signal.pch=c(19),
        file="tiff",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
CMplot(data,plot.type="q",box=FALSE,file="tif",file.name="",dpi=300,
    conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
    file.output=TRUE,verbose=TRUE,width=5,height=5)
#calculate Inflation factor (lambda) 
data$Chi2 <- qchisq(1 - data$dom, df = 1)
median_chi2_obs <- median(data$Chi2, na.rm = TRUE)
lambda <- median_chi2_obs / 0.455
cat("Inflation factor (lambda):", lambda, "\n")
##############################################

##extract and annotate gwas results with P < 1e-6
perl -e 'open (IN1,"$ARGV[0]");open (IN2,"$ARGV[1]");%add;%pos;while(<IN1>){chomp;@F=split/\t/,$_; next if ($F[3] == 0);if ($F[3] < 0.000001){$add{$F[1]}{$F[2]}=$F[3];$pos{$F[1]}{$F[2]}=1}};%dom;while(<IN2>){chomp;@F=split/\t/,$_;next if ($F[3] == 0);if ($F[3] < 0.000001){$dom{$F[1]}{$F[2]}=$F[3];$pos{$F[1]}{$F[2]}=1}};foreach $key1 (sort { $a <=> $b } keys %pos){foreach $key2 (sort { $a <=> $b } keys %{$pos{$key1}}){print "$key1\t$key2";if (exists $add{$key1}{$key2}){print "\t$add{$key1}{$key2}";}else{print "\tNA"};if (exists $dom{$key1}{$key2}){print "\t$dom{$key1}{$key2}\n";}else{print "\tNA\n"}}}' ../06.gcta/gwas.plot gwas.plot > gwas.1e-6.txt
awk '{print $1,$2}' OFS="\t" gwas.1e-6.txt > gwas.1e-6.pos
vcftools --gzvcf ../All.filtered.clean.vcf.gz --positions gwas.1e-6.pos --recode --out gwas.1e-6
perl /home/chao/software/annovar/table_annovar.pl gwas.1e-6.recode.vcf /home/chao/software/annovar/ovinedb/ -out gwas.1e-6 -buildver v3UseName -protocol refGene -operation g -nastring . -vcfinput
perl -e 'open (IN1,"$ARGV[0]");open (IN2,"$ARGV[1]");%hash;while(<IN1>){chomp;@F=split/\t/,$_,3;$hash{$F[0]}{$F[1]}=$F[2];};while(<IN2>){chomp;@F=split/\t/,$_;if (/^Chr/){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\tadditive\tdominant\n";}else{if (exists $hash{$F[0]}{$F[1]}){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$hash{$F[0]}{$F[1]}\n";}elsif(exists $hash{$F[0]}{$F[1]-1}){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$hash{$F[0]}{$F[1]-1}\n";}}}' gwas.1e-6.txt gwas.1e-6.v3UseName_multianno.txt > gwas.1e-6.anno

#LD heatmap
vcftools --vcf gwas.1e-6.recode.vcf --plink --out gwas.1e-6
plink --file gwas.1e-6 --make-bed --out gwas.1e-6
plink --bfile gwas.1e-6 --recode hv --out gwas.1e-6
plink --noweb --bfile gwas.1e-6 --r2 --ld-window-r2 0.1 --out gwas.1e-6 --chr-set 26 --allow-no-sex

#check the associated SNPs
grep -v "^##" gwas.1e-6.recode.vcf | perl -p -e "s/#//;s/:(.*?)\t/\t/g;s/:(.*?)$//" | perl -ne '{chomp;@F=split/\t/,$_,10;print "$F[0]\t$F[1]\t$F[9]\n"}' | perl -p -e "s/\.\/\./NA/g;s/0\/1/1/g;s/1\/1/2/g;s/0\/0/0/g;s/0\|1/1/g;s/1\|1/2/g" > gwas.1e-6.genotype.tab
Rscript add.r
Rscript dom.r