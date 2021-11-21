#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
#use Cwd 'abs_path';

#my $cwd = abs_path;
my ($file, $out, $database, $fasta, $mirbase, $help, $version);
GetOptions(
    "f|file:s"       =>  \$file,
    "o|out:s"        =>  \$out,
    "d|database:s"   =>  \$database,
    "fa|fasta:s"     =>  \$fasta,
    "m|mirbase:s"    =>  \$mirbase,
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
my $genome = join("/", $database, "Mus_musculus.GRCm38");
my $mature = join("/", $mirbase, "Mus_mature.fa");
my $hairpin = join("/", $mirbase, "Mus_hairpin.fa");


my $all_fa = join("/", $out_path, "all_samples.fa");
open(OT, "> $all_fa") or die "$! $all_fa\n";
foreach my $key (keys %file_name){
    my $prefix = join("/", $out_path, $key);
    if(-e $file_name{$key}){
        # merge all fa into onfile
        open(IN2, $file_name{$key}) or die "can't open $file_name{$key}\n";
        while(<IN2>){
            chomp;
            if($_ =~ /^>/){
                $_ =~ s/mmu/$key/;
                print OT "$_\n";
            }else{
                print OT "$_\n";
            }
        }
    }
}
close(OT);
open(OT2, "> $out") or die "can't open $out\n";
if(-e $all_fa){
    print OT2 "mapper.pl $all_fa -v -u -c -l 17 -r 10 -n -q -p $genome -t all\_aligned\.arf\n";
    print OT2 "miRDeep2.pl $all_fa $fasta all\_aligned\.arf $mature none $hairpin -t Mouse -P\n"; 
}
close(OT2);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -o <out> -d <database> -fa <primary fasta> -m <mirbase>
options:
	-f|file		     :[essential].
	-o|out           :[essential].
	-d|database		 :[essential].
	-fa|fasta        :[essential].
	-m|mirbase		 :[essential].
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.0
    update:     20200902 - 20200902
    author:     zou_hua\@grmh-gdl.cn
VERSION
};

