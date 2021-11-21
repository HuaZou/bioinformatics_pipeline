#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
#use Cwd 'abs_path';

#my $cwd = abs_path;
my ($file, $out, $database, $help, $version);
GetOptions(
    "f|file:s"       =>  \$file,
    "o|out:s"        =>  \$out,
    "d|database:s"   =>  \$database,
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
    $file_name{$tmp[1]} = $tmp[2];
}
close(IN1);

# database 
my($exon, $intron);
$exon   = join("/", $database, "Musculus_exon");
$intron = join("/", $database, "Musculus_intron");

open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    my $prefix = join("/", $out_path, $key);
    if(-e $file_name{$key}){
        print OT "bowtie -S -f -a -v 1 --best --strata -m 20 -p 2 --al $prefix\_for\_intron.fa --un $prefix\_no\_exon.fa -x $exon $file_name{$key} $prefix\_exon\_alignment.sam\n";
        print OT "rm $prefix\_exon\_alignment.sam\n";
        print OT "bowtie -S -f -a -v 1 --best --strata -m 20 -p 2 --al $prefix\_intron\_positive.fa -x $intron $prefix\_for\_intron.fa $prefix\_intron\_alignment.sam\n";
        print OT "rm $prefix\_intron\_alignment.sam\n";
        print OT "rm $prefix\_for\_intron.fa\n";
        print OT "cat $prefix\_no\_exon.fa $prefix\_intron\_positive.fa > $prefix\_exonintron\_filtered.fa\n";
        print OT "rm $prefix\_no\_exon.fa\n";        
        print OT "rm $prefix\_intron\_positive.fa\n";
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
	-d|database      :[essential].
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

