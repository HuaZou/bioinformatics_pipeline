#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use Cwd 'abs_path';

my ($file, $config, $prefix, $type, $out, $help, $version);
GetOptions(
    "f|file:s"     =>  \$file,
    "c|config:s"   =>  \$config,
    "p|prefix:s"   =>  \$prefix,
    "t|type:s"  =>  \$type,
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

my $Bin = $RealBin;
my $cwd = abs_path;

# output
my $dir = "$cwd/$out"; 
system "mkdir -p $dir" unless(-d $dir);
# scripts in bin & util
my $bin   = "$Bin/bin";
my $util  = "$Bin/util";

# configure file 
my %CFG;
open(CFG,"$config") or die "failed to open configure file $config. $!\n";
while(<CFG>){
    chomp;
    next if $_ =~ /^#/;next if $_ eq "";
	$_ =~ /^(\S+)\s*=\s*(\S+)/;
	$CFG{$1} = $2;
}
close(CFG);

####################################################
# star + featurecounts
if($type eq "star"){
    #################################################
    # scripts and output 
    # combine all steps in one script
    my $qc             = join("/", $dir, "Run.qc.sh");
    my $trim_galore    = join("/", $dir, "Run.trim_galore.sh");
    my $align_star     = join("/", $dir, "Run.align_star.sh");
    my $featureCounts  = join("/", $dir, "Run.featureCounts.sh");

    # scripts in bin & util
    my $s_qc             = "$bin/qc.pl";
    my $s_trim_galore    = "$bin/trim_galore.pl";
    my $s_align_star     = "$bin/align_star.pl";
    my $s_featureCounts  = "$bin/featureCounts.pl";
    my $s_get_star_prf   = "$util/get_star_profile.R";
    ##################################################
    #  step1 reads quality scan
    system("perl $s_qc -f $file -p $qc -o $dir");

    ##################################################
    #  step2 filter and trim low quality reads;
    system("perl $s_trim_galore -f $file -p $trim_galore -a $CFG{\"forward_adapter\"} -a2 $CFG{\"reverse_adapter\"} -q $CFG{\"trim_quality\"} -e $CFG{\"trim_error\"} -l $CFG{\"trim_length\"} -s $CFG{\"trim_stringency\"} -o $dir");

    ##################################################
    #  step3 alignment
    system("perl $s_align_star -f $file -p $align_star -i $CFG{\"db_STAR\"} -t $CFG{\"align_cpu\"} -o $dir");

    ##################################################
    #  step4 get gene expression profile
    system("perl $s_featureCounts -f $file -p $featureCounts -g $CFG{\"db_gtf\"} -t $CFG{\"align_cpu\"} -s $s_get_star_prf -c $CFG{\"occurrence\"} -nu $CFG{\"ncount\"} -o $dir");

    ##################################################
    # output 
    open(OT, "> $prefix") or die "can't open $prefix\n";
    print OT "sh $qc\n";
    print OT "sh $trim_galore\n";
    print OT "sh $align_star\n";
    print OT "sh $featureCounts\n";
    close(OT);
}
######################################################


######################################################
# hisat2 + stringtie
if($type eq "hisat2"){
    ########## output ################################
    # bash script per step
    # combine all steps in one script
    my $qc             = join("/", $dir, "Run.qc.sh");
    my $trim_flexbar   = join("/", $dir, "Run.flexbar.sh");
    my $align_hisat2   = join("/", $dir, "Run.align_hisat2.sh");
    my $stringtie_prf   = join("/", $dir, "Run.stringtie.sh");

    # scripts in bin & util
    my $s_qc             = "$bin/qc.pl";
    my $s_flexbar        = "$bin/flexbar.pl";
    my $s_align_hisat2   = "$bin/align_hisat2.pl";  
    my $s_stringtie_prf  = "$bin/stringtie.pl";
    my $s_get_profile    = "$util/stringtie_extract_profile.pl";
    my $s_get_matrix     = "$util/combineTable.pl";
    my $s_get_counts     = "$util/preDE.py";
    ##################################################
    #  step1 reads quality scan
    system("perl $s_qc -f $file -p $qc -o $dir");

    ##################################################
    #  step2 filter and trim low quality reads;
    system("perl $s_flexbar -f $file -p $trim_flexbar -a $CFG{\"all_adapter\"} -q $CFG{\"trim_quality\"} -u $CFG{\"trim_uncall\"} -l $CFG{\"trim_length\"} -t $CFG{\"trim_cpu\"} -o $dir");

    ##################################################
    #  step3 alignment
    system("perl $s_align_hisat2 -f $file -p $align_hisat2 -i $CFG{\"db_hisat2\"} -t $CFG{\"align_cpu\"} -o $dir");

    ##################################################
    #  step4 get gene expression profile 
    system("perl $s_stringtie_prf -f $file -p $stringtie_prf -g $CFG{\"db_gtf\"} -t $CFG{\"align_cpu\"} -s $s_get_profile -s2 $s_get_matrix -s3 $s_get_counts -o $dir");

    ##################################################
    # output 
    open(OT, "> $prefix") or die "can't open $prefix\n";
    print OT "sh $qc\n";
    print OT "sh $trim_flexbar\n";
    print OT "sh $align_hisat2\n";
    print OT "sh $stringtie_prf\n";    
    close(OT);
}
#######################################################


#######################################################
# kallisto + salmon
if($type eq "kallisto"){
    ########## output #################################
    # bash script per step
    # combine all steps in one script
    my $qc             = join("/", $dir, "Run.qc.sh");
    my $trim_flexbar   = join("/", $dir, "Run.flexbar.sh");
    my $align_kallisto = join("/", $dir, "Run.align_kallisto.sh");
    my $align_salmon   = join("/", $dir, "Run.align_salmon.sh");
    my $get_cdna_prf   = join("/", $dir, "Run.get_profile.sh");

    # scripts in bin & util
    my $s_qc             = "$bin/qc.pl";
    my $s_flexbar        = "$bin/flexbar.pl";
    my $s_align_kallisto = "$bin/align_kallisto.pl";
    my $s_align_salmon   = "$bin/align_salmon.pl";
    my $s_get_cdna_prf   = "$util/get_kallisto_profile.R";
    ##################################################
    #  step1 reads quality scan
    system("perl $s_qc -f $file -p $qc -o $dir");

    ##################################################
    #  step2 filter and trim low quality reads;
    system("perl $s_flexbar -f $file -p $trim_flexbar -a $CFG{\"all_adapter\"} -q $CFG{\"trim_quality\"} -u $CFG{\"trim_uncall\"} -l $CFG{\"trim_length\"} -t $CFG{\"trim_cpu\"} -o $dir");

    ##################################################
    #  step3 alignment
    system("perl $s_align_kallisto -f $file -p $align_kallisto -i $CFG{\"db_kallisto\"} -g $CFG{\"db_gtf\"} -t $CFG{\"align_cpu\"} -s $s_get_cdna_prf -tx $CFG{\"tx2gene\"} -c $CFG{\"occurrence\"} -n $CFG{\"ncount\"} -o $dir");    
    system("perl $s_align_salmon -f $file -p $align_salmon -i $CFG{\"db_salmon\"} -t $CFG{\"align_cpu\"} -s $s_get_cdna_prf -tx $CFG{\"tx2gene\"} -c $CFG{\"occurrence\"} -n $CFG{\"ncount\"} -o $dir");

    ##################################################
    #  step4 get gene expression profile 
    

    ##################################################
    # output 
    open(OT, "> $prefix") or die "can't open $prefix\n";
    print OT "sh $qc\n";
    print OT "sh $trim_flexbar\n";
    print OT "sh $align_kallisto\n";
    print OT "sh $align_salmon\n";    
    close(OT);
}
#######################################################


sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -c <config> -p <prefix> -t <type> -o <out>
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -c|config   :[essential]. the parameters 
    -p|prefix   :[essential]. the name of all processures
    -t|type     :[essential]. the type of alignment   
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
