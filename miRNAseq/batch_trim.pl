#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
#use Cwd 'abs_path';

#my $cwd = abs_path;
my ($file, $out, $adapter, $adapter2, $help, $version);
GetOptions(
    "f|file:s"       =>  \$file,
    "a|adapter:s"    =>  \$adapter,
    "a2|adapter2:s"  =>  \$adapter2,
    "o|out:s"        =>  \$out,
    "h|help:s"       =>  \$help,
    "v|version"      =>  \$version
);
&usage if(!defined $out);
my $out_name = (split(".sh", $out))[0];
my $out_path = join("/", "./result", $out_name);
my $out_path_trim = join("/", $out_path, "cutadapt");
my $out_path_filter = join("/", $out_path, "trimmomatic");
system "mkdir -p $out_path $out_path_trim $out_path_filter" unless(-d $out_path && $out_path_trim && $out_path_filter);


open(IN1, $file) or die "can't open $file";
my %file_name;
<IN1>;
while(<IN1>){
    chomp;
    my @tmp = split("\t", $_);
    $file_name{$tmp[1]} = $tmp[2];
}
close(IN1);


open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    my $new_key = $1 if $key =~ /(\S+)\_miRNA/;
    #print "$new_key\t$file_name{$key}\n";
    my $prefix_trim = join("/", $out_path_trim, $new_key);
    my $prefix_filter = join("/", $out_path_filter, $new_key);
    if(-e $file_name{$key}){
        if($key =~ /novo$/){
            print OT "cutadapt -a $adapter2 --minimum-length=17 --maximum-length=35 -o $prefix_trim\.trimmed.fq $file_name{$key}\n";    
        }else{
            print OT "cutadapt -a $adapter --minimum-length=17 --maximum-length=35 -o $prefix_trim\.trimmed.fq $file_name{$key}\n";
        }
        print OT "trimmomatic SE -threads 4 -phred33 $prefix_trim\.trimmed.fq $prefix_filter\.clean.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:17\n";
        print OT "rm $prefix_trim\.trimmed.fq\ngzip $prefix_filter\.clean.fq\nrm $prefix_filter\.clean.fq\n";
    }
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -o <out> -a <adapter> -a2 <adapter2>
options:
	-f|file		     :[essential].
	-o|out           :[essential].
    -a|adapter       :[essential] adapter.
    -a2|adapter2     :[essential] adapter2.
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.3
    update:     20200831 - 20200831
    author:     zou_hua\@grmh-gdl.cn
VERSION
};

