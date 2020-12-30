#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $real_dir, $out, $help, $version);
GetOptions(
    "f|file:s"  =>  \$file,
    "d|real_dir:s"  => \$real_dir,  
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

my $dir = "$real_dir/result/03.humann"; 
system "mkdir -p $dir" unless(-d $dir);

my $dir_log = "$real_dir/result/03.humann/log"; 
system "mkdir -p $dir_log" unless(-d $dir_log);

my $genefamilies = "$real_dir/result/03.humann/genefamilies"; 
system "mkdir -p $genefamilies" unless(-d $genefamilies);
my $pathabundance = "$real_dir/result/03.humann/pathabundance"; 
system "mkdir -p $pathabundance" unless(-d $pathabundance);
my $pathcoverage = "$real_dir/result/03.humann/pathcoverage"; 
system "mkdir -p $pathcoverage" unless(-d $pathcoverage);

my $dir_metaphlan = "$real_dir/result/03.humann/metaphlan"; 
system "mkdir -p $dir_metaphlan" unless(-d $dir_metaphlan);

# script
my $dir_script = "$real_dir/result/script/03.humann/"; 
system "mkdir -p $dir_script" unless(-d $dir_script);

open(IN, $file) or die "can't open $file";
my %file_name;
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    push (@{$file_name{$tmp[0]}}, $tmp[2]);
}
close(IN);

my ($fq);
open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    $fq = join("", "./result/02.merge/", $key, "_merge.fastq.gz");
    #if($fq){

        my $bash = join("", $dir_script, $key, ".humann.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";         
        print OT2 "humann --input $fq --output $dir --threads 10\n";
        print OT2 "mv $dir/$key\_merge_humann_temp/$key\_merge_metaphlan_bugs_list.tsv $dir_metaphlan/$key\_metaphlan.tsv\n";
        print OT2 "mv $dir/$key\_merge_humann_temp/$key\_merge.log $dir_log\n";
        print OT2 "mv $dir/$key\_merge_genefamilies.tsv $genefamilies/$key\_genefamilies.tsv\n";
        print OT2 "mv $dir/$key\_merge_pathabundance.tsv $pathabundance/$key\_pathabundance.tsv\n";
        print OT2 "mv $dir/$key\_merge_pathcoverage.tsv $pathcoverage/$key\_pathcoverage.tsv\n";
        print OT2 "rm -r $dir/$key\_merge_humann_temp/\n";
        close(OT2);

        print OT "sh $bash\n";         
    #}
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -d <real_dir> -o <out> 
options:
	-f|file		:[essential].
    -d|real_dir :[essential].    
	-o|out      :[essential].
USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.1
    update:     20201223 - 20201224
    author:     zouhua1\@outlook.com
VERSION
};
