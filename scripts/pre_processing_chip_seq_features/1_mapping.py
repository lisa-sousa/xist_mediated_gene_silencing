import os
import argparse

# before using this script download sra files with SRA toolkit
# for information on the parameters type python mapping.py -h
# example call: python mapping.py -i SRR1991254.sra -g GSE68195 -n input -s 0 -d mapping/

# set the correct path to the tools --> set path to tools to "" if you can call them directly!!!
ucsc_tools = "/scratch/ngsvin/bin/chip-seq/bigwig/" 
bedTools = "/home/lisasous/tools/bedtools2/bin/"
bowtie2 = "/home/lisasous/tools/bowtie2-2.3.4.3/"
samtools = "/home/lisasous/tools/samtools-1.9/"

# set the correct path to the chromosome files!!!
chrom_file = "/project/lncrna/Xist/data/annotation_files/mouse_genome/mm9.chrom.sizes"
ORG_ASSEMBLY="/project/epigenome/index/mm9"

#============================================================
#read terminal parameters
#============================================================

def args():
  args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  args_parser.add_argument("-i","--input_file",help="sra or fastq file",type=str,required=True)
  args_parser.add_argument("-g","--geo_number",help="GEO accession number",type=str,required=True)
  args_parser.add_argument("-n","--name",help="name of control or experiment, e.g. TF or HM",type=str,required=True)
  args_parser.add_argument("-s","--input_source",help="sra or fastq file -> 0:sra; 1: fastq",type=int,required=True)
  args_parser.add_argument("-np","--number_processors",help="number of cores that should be used by bowtie for mapping (default: 4)",type=int,default=4)
  args_parser.add_argument("-d","--output_dir",help="output directory",type=str,required=True)
  return args_parser.parse_args()
 
#===============================================================
#Methods
#===============================================================

# function to convert from sra to fastq
def sra_to_fastq(sra_file,fastq_dir):
    cmd = "fastq-dump --split-3 " + sra_file + " --outdir " + fastq_dir
    print('\n' + cmd)
    os.system(cmd)

#============================================================

# Mapping reads with bowite 2, allow up to 1 mismatch and keep only high quality reads. Convert sam to bam
def bowtie_mapping(fastq_file,mapping_file_name,number_processors,bowtie2):

    # Bowtie parameters
    V="1" # number of mismatches
    cmd = bowtie2 + "bowtie2 -p " + str(number_processors) + " -N " + V + " -x " + ORG_ASSEMBLY + " -U " + fastq_file + " -S " + mapping_file_name + ".fastq.bed.sam"
    print('\n' + cmd)
    os.system(cmd)

    # convert to bam & take only reads above a certain quality specified by option -q
    cmd = samtools + "samtools view -b -q 10 -S "+ mapping_file_name + ".fastq.bed.sam" + " > " + mapping_file_name + ".fastq.bed.bam"
    print('\n' + cmd)
    os.system(cmd)

    # count mapped reads: samtools view -F 0x904 -c *.fastq.bed.bam
    # sort bam file
    cmd = samtools + "samtools sort -o " + mapping_file_name + ".fastq.bed.sorted.bam" + " " + mapping_file_name + ".fastq.bed.bam"
    print('\n' + cmd)
    os.system(cmd)

    # create the .bai index, it adds ".bai" automatically
    cmd = samtools + "samtools index " + mapping_file_name + ".fastq.bed.sorted.bam"
    print('\n' + cmd)
    os.system(cmd)

    cmd = "rm -rf " + mapping_file_name + ".fastq.bed.sam"
    print('\n' + cmd)
    os.system(cmd)
    
#============================================================

def fastq_to_bw(fastq_file, output_file_name, mapping_dir, bedgraph_dir, bigwig_dir, number_processors,ucsc_tools, bedTools, bowtie2):

    # mapping
    mapping_file_name = mapping_dir + output_file_name
    bowtie_mapping(fastq_file,mapping_file_name,number_processors,bowtie2)

    # convert bam to BedGraph
    bedgraph_file = bedgraph_dir + output_file_name + ".fastq.bed.sorted.bedgraph"
    cmd = bedTools + "genomeCoverageBed -split -bg -ibam " + mapping_file_name + ".fastq.bed.sorted.bam"  + " -g " + chrom_file + " > " + bedgraph_file
    print('\n' + cmd)
    os.system(cmd)

    # Convert the BedGraph file to BigWig
    cmd = ucsc_tools + "bedGraphToBigWig " + bedgraph_file + " " + chrom_file + " " + bigwig_dir + output_file_name + ".fastq.bed.sorted.bw"
    print('\n' + cmd)
    os.system(cmd)


#============================================================     
#File Parameters   
#============================================================     
      
parameters = args()

input_file = parameters.input_file
geo_number = parameters.geo_number
name = parameters.name
input_source = parameters.input_source
number_processors = parameters.number_processors
output_dir = parameters.output_dir

#=============================================================
#executive part
#=============================================================

if not output_dir.endswith("/"):
  output_dir = output_dir + "/"

fastq_dir = output_dir + "fastq/"
mapping_dir = output_dir + "bam/"
bedgraph_dir = output_dir + "bedgraph/"
bigwig_dir = output_dir + "bigwig/"

if not os.path.isdir(fastq_dir):
    cmd = "mkdir " + fastq_dir
    os.system(cmd)

if not os.path.isdir(mapping_dir):
    cmd = "mkdir " + mapping_dir
    os.system(cmd)

if not os.path.isdir(bedgraph_dir):
    cmd = "mkdir " + bedgraph_dir
    os.system(cmd)

if not os.path.isdir(bigwig_dir):
    cmd = "mkdir " + bigwig_dir
    os.system(cmd)

output_file_name = name + "_" + geo_number

if input_source == 0:
    sra_file = input_file
    sra_to_fastq(sra_file,fastq_dir)
    fastq_file = sra_file.split("/")
    fastq_file = fastq_dir + fastq_file[len(fastq_file)-1].replace(".sra",".fastq")
else:
    fastq_file = input_file


fastq_to_bw(fastq_file, output_file_name, mapping_dir, bedgraph_dir, bigwig_dir, number_processors,ucsc_tools, bedTools, bowtie2)
