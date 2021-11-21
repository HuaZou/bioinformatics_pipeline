#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;


my ($file, $prefix, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_qc = "$out/quality/fastqc"; 
system "mkdir -p $dir_qc" unless(-d $dir_qc);

# script
my $dir_script = "$out/script/quality/"; 
system "mkdir -p $dir_script" unless(-d $dir_script);

my @array_name;
open(IN, $file) or die "can't open $file\n";
open(OT, "> $prefix") or die "can't open $prefix\n";
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    if(-e $tmp[2]){
        my $bash = join("", $dir_script, $tmp[1], ".fastqc.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";
        print OT2 "fastqc -o $dir_qc --noextract $tmp[2]\n";
        close(OT2);

        print OT "sh $bash\n";
    }
}
close(IN);

my $dir_mc = "$out/quality/multiqc"; 
system "mkdir -p $dir_mc" unless(-d $dir_mc);
print OT "multiqc $dir_qc --outdir $dir_mc\n";
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
