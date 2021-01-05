#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $index, $gtf, $threads, $script, $tx2gen, $occur, $number, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "i|index:s"     =>  \$index,
    "g|gtf:s"       =>  \$gtf,
    "t|threads:s"   =>  \$threads,
    "s|script:s"    =>  \$script,
    "tx|tx2gen:s"   =>  \$tx2gen,
    "c|occur:s"     =>  \$occur,
    "n|number:s"    =>  \$number,            
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_align = "$out/align_kallisto/"; 
system "mkdir -p $dir_align" unless(-d $dir_align);

# script
my $dir_script = "$out/script/align_kallisto/"; 
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

my ($fq1, $fq2, $abundance);
my $abundance_list = join("/", $dir_align, "abundance_list.tsv");
open(AB, "> $abundance_list") or die "can't open $abundance_list\n";
print AB "SampleID\tPath\n";
open(OT, "> $prefix") or die "can't open $prefix\n";
foreach my $key (keys %file_name){
    $fq1 = join("", "$out/flexbar/", $key, "_1.fastq.gz");
    $fq2 = join("", "$out/flexbar/", $key, "_2.fastq.gz");
    $abundance = join("", $dir_align, $key, "/abundance.tsv"); 
    print AB "$key\t$abundance\n";  
    my $bash = join("", $dir_script, $key, ".align_kallisto.sh");
    open(OT2, "> $bash") or die "can't open $bash\n"; 
    print OT2 "kallisto quant -i $index -g $gtf -o $dir_align/$key -b 30 -t $threads $fq1 $fq2\n"; 
    close(OT2);
    
    print OT "sh $bash\n";        
}

print OT "Rscript $script -s $abundance_list -d $dir_align -g $tx2gen -t kallisto -o $dir_align\n";

close(OT);
close(AB);

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
