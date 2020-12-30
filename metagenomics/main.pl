#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use Cwd 'abs_path';

my ($file, $adapter, $out, $help, $version);
GetOptions(
    "f|file:s"  =>  \$file,
    "o|out:s"   =>  \$out,
    "a|adapter:s"   =>  \$adapter,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

my $Bin = $RealBin;
my $cwd = abs_path;

# output
my $dir = "$cwd/result/"; 
system "mkdir -p $dir" unless(-d $dir);

########## output #########################################
# bash script per step
# combine all steps in one script
my $qc        = join("", $dir, "Run.s1.qc.sh");
my $kneaddata = join("", $dir, "Run.s2.kneaddata.sh");
my $merge     = join("", $dir, "Run.s3.merge.sh");
my $humann    = join("", $dir, "Run.s4.humann.sh");
my $metaphlan = join("", $dir, "Run.s5.metaphlan.sh");


# scripts in bin
my $bin         = "$Bin/bin";
my $s_qc        = "$bin/qc.pl";
my $s_kneaddata = "$bin/kneaddata.pl";
my $s_merge     = "$bin/merge.pl";
my $s_humann    = "$bin/humann.pl";
my $s_metaphlan = "$bin/metaphlan.pl";

########## Steps in metagenomics pipeline #################
##################################################
#  step1 reads quality scan
system("perl $s_qc -f $file -d $cwd -o $qc");

##################################################
#  step2 filter and trim low quality reads;
#        remove host sequence
system("perl $s_kneaddata -f $file -d $cwd -a $adapter -o $kneaddata");

##################################################
#  step3 merge PE reads
system("perl $s_merge -f $file -d $cwd -o $merge");

##################################################
#  step4 get function profile
system("perl $s_humann -f $file -d $cwd -o $humann");

##################################################
#  step5 get taxonomy profile
system("perl $s_metaphlan -f $file -d $cwd -o $metaphlan");


open(OT, "> $out") or die "can't open $out\n";
print OT "sh $qc\nsh $kneaddata\nsh $merge\nsh $humann\nsh $metaphlan\n";
close(OT);


sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -o <out> -a <adapter>
options:
	-f|file     :[essential].
	-o|out      :[essential].
	-a|adapter  :[essential].
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.1
    update:     20201224 - 20201225
    author:     zouhua1\@outlook.com
VERSION
};
