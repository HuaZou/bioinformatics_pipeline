#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $real_dir,  $out, $help, $version);
GetOptions(
    "f|file:s"  =>  \$file, 
    "d|real_dir:s"  => \$real_dir,  
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

my $dir = "$real_dir/result/02.merge"; 
system "mkdir -p $dir" unless(-d $dir);

# script
my $dir_script = "$real_dir/result/script/02.merge/"; 
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
open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    $fq1 = join("", "./result/01.kneaddata/", $key, "_paired_1.fastq");
    $fq2 = join("", "./result/01.kneaddata/", $key, "_paired_2.fastq");
    #if(-e $fq1 && -e $fq2){
        #print OT "fastp -i $fq1 -I $fq2 -h $dir/$key\_merge.html -j $dir/$key\_merge.json -m --merged_out $dir/$key\_merge.fastq.gz --failed_out  $dir/$key\_failed.fastq.gz --include_unmerged --overlap_len_require 6 --overlap_diff_percent_limit 20 --detect_adapter_for_pe -5 -r -l 20 -y --thread 5\n";
        my $bash = join("", $dir_script, $key, ".merge.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";        
        print OT2 "fastp -i $fq1 -I $fq2 -h $dir/$key\_merge.html -j $dir/$key\_merge.json -m --merged_out $dir/$key\_merge.fastq.gz --failed_out  $dir/$key\_failed.fastq.gz --include_unmerged --overlap_len_require 6 --overlap_diff_percent_limit 20 --detect_adapter_for_pe -5 -r -l 20 -y --thread 5\n";
        close(OT2);

        print OT "sh $bash\n";         
    #}
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -d <real_dir> -o <out> 
options:
        -f|file	    :[essential].
        -d|real_dir :[essential].    
        -o|out      :[essential].
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.1
    update:     20201223 - 20201224
    author:     zouhua1\@outlook.com
VERSION
};
