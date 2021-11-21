#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,     
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);


open(OT, "> $prefix") or die "can't open $prefix\n";
open(IN, "$file") or die "can't open $file\n";    
<IN>;
print OT "Gene_id\tCoverage\tFPKM\tTPM\n";
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    my $result = join("\t", $tmp[0], $tmp[6], $tmp[7], $tmp[8]);
    print OT "$result\n";
}
close(IN);
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> 
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure   

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
