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

my $dir = "$real_dir/result/04.metaphlan"; 
system "mkdir -p $dir" unless(-d $dir);

# script
my $dir_script = "$real_dir/result/script/04.metaphlan/"; 
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

        my $bash = join("", $dir_script, $key, ".metaphlan.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";         
        print OT2 "metaphlan $fq --bowtie2out $dir/$key\_metagenome.bowtie2.bz2 --nproc 10 --input_type fastq -o $dir/$key\_metagenome.tsv --unknown_estimation -t rel_ab_w_read_stats\n";
        close(OT2);

        print OT "sh $bash\n";         
    #}
}

print OT "merge_metaphlan_tables.py $dir/*metagenome.tsv > $dir/merge_metaphlan.tsv\n";
print OT "Rscript /data/share/database/metaphlan_databases/calculate_unifrac.R $dir/merge_metaphlan.tsv /data/share/database/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk $dir/unifrac_merged_mpa3_profiles.tsv";
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
