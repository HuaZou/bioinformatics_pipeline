#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $f_adapter, $r_adapter, $quality, $error, $length, $stringency, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "a|f_adapter:s"    =>  \$f_adapter,
    "a2|r_adapter:s"   =>  \$r_adapter,
    "q|quality:s"      =>  \$quality,
    "e|error:s"        =>  \$error,
    "l|length:s"       =>  \$length,
    "s|stringency:s"   =>  \$stringency,            
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_trim = "$out/trim_galore/"; 
system "mkdir -p $dir_trim" unless(-d $dir_trim);

# script
my $dir_script = "$out/script/trim_galore/"; 
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
    if (${$file_name{$key}}[0] =~ /[r|R]1/){
        $fq1 = ${$file_name{$key}}[0];
    }else{
        $fq2 = ${$file_name{$key}}[0];
    }

    if (${$file_name{$key}}[1] =~ /[r|R]2/){
        $fq2 = ${$file_name{$key}}[1];
    }else{
        $fq1 = ${$file_name{$key}}[1];
    }
    my $bash = join("", $dir_script, $key, ".trim_galore.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "trim_galore --paired  $fq1 $fq2 --quality $quality -e $error --length $length --stringency $stringency -o $dir_trim -a $f_adapter -a2 $r_adapter --fastqc\n";       
    close(OT2);
    
    print OT "sh $bash\n";        
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -a <f_adapter> -a2 <r_adapter> -q <quality> -e <error> -l <length> -s <stringency> -o <out> 
options:
	-f|file       :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix     :[essential]. the name of this processure   
    -a|f_adapter  :[essential]. forward adapter sequence
    -a2|r_adapter :[essential]. reverse adapter sequence
    -q|quality    :[essential]. reads quality
    -e|error      :[essential]. error rate 
    -l|length     :[essential]. minus length
    -s|stringency :[essential]. stringency 
    -o|out        :[essential]. output directory
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
