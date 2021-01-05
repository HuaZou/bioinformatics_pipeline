#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $index, $threads, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "i|index:s"     =>  \$index,
    "t|threads:s"   =>  \$threads,            
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_align = "$out/align_hisat2/"; 
system "mkdir -p $dir_align" unless(-d $dir_align);

# script
my $dir_script = "$out/script/align_hisat2/"; 
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

my ($fq1, $fq2);
open(OT, "> $prefix") or die "can't open $prefix\n";
foreach my $key (keys %file_name){
    $fq1 = join("", "$out/flexbar/", $key, "_1.fastq.gz");
    $fq2 = join("", "$out/flexbar/", $key, "_2.fastq.gz");    
    my $bash = join("", $dir_script, $key, ".align_hisat2.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "hisat2 -x $index -1 $fq1 -2 $fq2 -S $dir_align/$key\.sam -p $threads --rg-id=$key\n";
    print OT2 "samtools view -@ 16 -b -S $dir_align/$key\.sam -o $dir_align/$key\.bam\n";
    print OT2 "samtools sort -@ $threads $dir_align/$key\.bam $dir_align/$key\.sorted\n";
    print OT2 "rm -rf $dir_align/$key\.sam $dir_align/$key\.bam\n";     
    close(OT2);
    
    print OT "sh $bash\n";        
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -o <out>
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure   
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
