#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $gtf, $threads, $script, $script2, $script3, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "g|gtf:s"       =>  \$gtf,
    "t|threads:s"   =>  \$threads,
    "s|script:s"    =>  \$script,
    "s2|script2:s"  =>  \$script2, 
    "s3|script3:s"  =>  \$script3,              
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_count = "$out/stringtie/"; 
system "mkdir -p $dir_count" unless(-d $dir_count);

# script
my $dir_script = "$out/script/stringtie/"; 
system "mkdir -p $dir_script" unless(-d $dir_script);

open(IN, $file) or die "can't open $file";
my %file_name;
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    push (@{$file_name{$tmp[0]}}, $tmp[2]);
}
close(IN);

my $gtf_list = join("/", $dir_count, "gtf_list.tsv");
my $tsv_list = join("/", $dir_count, "tsv_list.tsv");
open(GTF, "> $gtf_list") or die "can't open $gtf_list\n"; 
open(TSV, "> $tsv_list") or die "can't open $tsv_list\n"; 
open(OT, "> $prefix") or die "can't open $prefix\n";
foreach my $key (keys %file_name){
    my $bam = join("", "$out/align_hisat2/", $key, ".sorted.bam");
    my $bash = join("", $dir_script, $key, ".stringtie.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "stringtie --rf -p $threads -G $gtf -e -B -o $dir_count/$key\.gtf -A $dir_count/$key\.tsv $bam\n";
    print OT2 "perl $script -f $dir_count/$key\.tsv -p $dir_count/$key\.extract.tsv\n";
    close(OT2);

    print GTF "$key $dir_count/$key\.gtf\n";
    print TSV "$key $dir_count/$key\.extract.tsv\n";  
    print OT "sh $bash\n";
}

# FPKM & TPM
my $coverage = join("/", $dir_count, "stringtie.coverage");
my $FPKM = join("/", $dir_count, "stringtie.FPKM");
my $TPM = join("/", $dir_count, "stringtie.TPM");
print OT "perl $script2 $tsv_list 1 $coverage  1\n";
print OT "perl $script2 $tsv_list 2 $FPKM  1\n";
print OT "perl $script2 $tsv_list 3 $TPM  1\n";

# count
my $gene_counts = join("/", $dir_count, "stringtie_gene_counts.csv");
my $transcript_counts = join("/", $dir_count, "stringtie_transcript_counts.csv");
print OT "source activate py27\n";
print OT "python $script3 -i $gtf_list -g $gene_counts -t $transcript_counts\n";

close(GTF);
close(TSV);
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -g <gtf> -t <threads> -s <script> -s2 <script> -o <out> 
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure   
    -g|gtf      :[essential]. gtf file
    -t|threads  :[essential]. threads
    -s|script   :[essential]. script
    -s2|script2 :[essential]. script
    -s3|script3 :[essential]. script
	-o|out      :[essential]. output directory

USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.0
    update:     20201231 - 20201231
    author:     zouhua1\@outlook.com
VERSION
};
