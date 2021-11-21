#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
#use Cwd 'abs_path';

#my $cwd = abs_path;
my ($file, $out, $help, $version);
GetOptions(
    "f|file:s"       =>  \$file,   
    "o|out:s"        =>  \$out,
    "h|help:s"       =>  \$help,
    "v|version"      =>  \$version
);
&usage if(!defined $out);
my $out_name = (split(".sh", $out))[0];
my $out_path = join("/", "./result", $out_name);
system "mkdir -p $out_path" unless(-d $out_path);


open(IN1, $file) or die "can't open $file";
my %file_name;
<IN1>;
while(<IN1>){
    chomp;
    my @tmp = split("\t", $_);
    print "$tmp[1]\n";
    $file_name{$tmp[1]} = $tmp[2];
}
close(IN1);

my ($fq);
open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    if(-e $file_name{$key}){
        print OT "fastqc  -o $out_path --noextract -t 2 $file_name{$key}\n";
    } 
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -o <out> 
options:
	-f|file		     :[essential].
	-o|out           :[essential].

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

