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
my $out_path_fasta = join("/", $out_path, "fasta");
my $out_path_collapse = join("/", $out_path, "collapse");
system "mkdir -p $out_path $out_path_fasta $out_path_collapse" unless(-d $out_path && $out_path_fasta && $out_path_collapse);


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
    my $new_key = $1 if $key =~ /(\S+)\.clean/;
    #print "$new_key\t$file_name{$key}\n";
    my $prefix_fasta = join("/", $out_path_fasta, $new_key);
    my $prefix_collapse = join("/", $out_path_collapse, $new_key);
    if(-e $file_name{$key}){
        print OT "seqtk seq -a  $file_name{$key} > $prefix_fasta\.fa\n";
        print OT "collapse_reads_md.pl $prefix_fasta\.fa mmu > $prefix_collapse\.collapse\.fa\n";
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
    version:    v1.0
    update:     20200901 - 20200901
    author:     zou_hua\@grmh-gdl.cn
VERSION
};

