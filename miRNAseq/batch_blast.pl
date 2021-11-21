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

# RNA database 
my($rRNA, $sRNA, $tRNA, $snRNA);
$rRNA = join("/", $database, "rRNA_Rfam_Mus");
$sRNA = join("/", $database, "sRNA_Rfam_Mus");
$tRNA = join("/", $database, "tRNA_Rfam_Mus");
$snRNA = join("/", $database, "snRNA_Rfam_Mus");

open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    my $prefix = join("/", $out_path, $key);
    if(-e $file_name{$key}){
        my $fa = join(".", $prefix, "remained.fa");
        print OT "blastn -task blastn-short -query $file_name{$key} -out $prefix\_rRNA.blast -db $rRNA -outfmt 6 -evalue 0.01 -num_threads 4\n";
        print OT "blastn -task blastn-short -query $file_name{$key} -out $prefix\_sRNA.blast -db $sRNA -outfmt 6 -evalue 0.01 -num_threads 4\n"; 
        print OT "blastn -task blastn-short -query $file_name{$key} -out $prefix\_tRNA.blast -db $tRNA -outfmt 6 -evalue 0.01 -num_threads 4\n";
        print OT "blastn -task blastn-short -query $file_name{$key} -out $prefix\_snRNA.blast -db $snRNA -outfmt 6 -evalue 0.01 -num_threads 4\n";
        print OT "cat $prefix\*.blast | awk \'{print \$1}\' | sort | uniq > $prefix\.annotated_RNA.txt\n";
        print OT "sed -n \'/^>/p\' $file_name{$key} | sed \'s/^>//g\' | sort > $prefix\.all_seq.txt\n";
        print OT "sort $prefix\.all_seq.txt $prefix\.annotated_RNA.txt $prefix\.annotated_RNA.txt | uniq -u > $prefix\.unannotated_RNA.txt\n";
        print OT "python get_unannotated_fasta.py -f $file_name{$key} -d $prefix\.unannotated_RNA.txt -o $prefix\.remained.fa\n";
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

