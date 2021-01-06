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
my $profile   = join("", $dir, "Run.s6.profile.sh");

# scripts in bin
my $bin         = "$Bin/bin";
my $util        = "$Bin/util";
my $s_qc        = "$bin/qc.pl";
my $s_kneaddata = "$bin/kneaddata.pl";
my $s_merge     = "$bin/merge.pl";
my $s_humann    = "$bin/humann.pl";
my $s_metaphlan = "$bin/metaphlan.pl";
my $s_profile   = "$bin/convert2matrix.pl";
my $s_unifrac   = "$util/calculate_unifrac.R";

# mpa_v30_CHOCOPhlAn_201901_species_tree.nwk
my $species_tree = "$util/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk";


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

##################################################
#  step5 get taxonomy profile
system("perl $s_metaphlan -f $file -d $cwd -o $metaphlan");

##################################################
#  step6 get profile matrix
system("perl $s_profile -s $s_unifrac -t $species_tree -d $cwd -o $profile");

open(OT, "> $out") or die "can't open $out\n";
print OT "sh $qc\n";
print OT "sh $kneaddata\n";
print OT "sh $merge\n";
print OT "sh $humann\n";
print OT "sh $metaphlan\n";
print OT "sh $profile\n";
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
