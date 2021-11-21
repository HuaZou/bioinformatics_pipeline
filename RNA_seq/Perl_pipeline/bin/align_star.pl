#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $index, $threads, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "i|index:s"   =>  \$index,
    "t|threads:s"   =>  \$threads,            
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_align = "$out/align_star/"; 
system "mkdir -p $dir_align" unless(-d $dir_align);

# script
my $dir_script = "$out/script/align_star/"; 
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
    $fq1 = join("", "$out/trim_galore/", $key, ".r1_val_1.fq.gz");
    $fq2 = join("", "$out/trim_galore/", $key, ".r2_val_2.fq.gz");    
    my $bash = join("", $dir_script, $key, ".align_star.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "STAR --genomeDir $index --readFilesIn $fq1 $fq2 --runThreadN $threads --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $dir_align/$key --quantMode GeneCounts --readFilesCommand zcat\n";       
    close(OT2);
    
    print OT "sh $bash\n";        
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -i <index> -o <out>
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure
    -i|index    :[essential]. genome index builded by STAR  
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
