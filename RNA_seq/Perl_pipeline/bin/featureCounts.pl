#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $gtf, $threads, $rscript, $occurrence, $number, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "g|gtf:s"       =>  \$gtf,
    "t|threads:s"   =>  \$threads,
    "s|rscript:s"   =>  \$rscript,
    "c|occurrence:s"   =>  \$occurrence,
    "nu|number:s"   =>  \$number,               
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_count = "$out/featureCounts/"; 
system "mkdir -p $dir_count" unless(-d $dir_count);

# script
my $dir_script = "$out/script/featureCounts/"; 
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

my @sam;
open(OT, "> $prefix") or die "can't open $prefix\n";
foreach my $key (keys %file_name){
    my $temp = join("", "$out/align_star/", $key, "Aligned.sortedByCoord.out.bam");
    push(@sam, $temp);         
}
my $sam_list = join(" ", @sam);
my $bash = join("", $dir_script, "featureCounts.sh");
open(OT2, "> $bash") or die "can't open $bash\n";
print OT2 "featureCounts -a $gtf -o $dir_count/gene_count_table.tsv -t exon -T $threads -g gene_id --primary $sam_list\n"; 
print OT2 "Rscript $rscript -f $dir_count/gene_count_table.tsv -c $occurrence -n $number -o $dir_count\n";
close(OT2);
    
print OT "sh $bash\n"; 

close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -g <gtf> -t <threads> -s <rscript> -c <occurrence> -nu <number> -o <out> 
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure   
    -g|gtf      :[essential]. gtf file
    -t|threads  :[essential]. threads
    -s|rscript  :[essential]. rscript
    -c|occurrence  :[essential]. occurrence
    -nu|number     :[essential]. number
	-o|out         :[essential]. output directory

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
