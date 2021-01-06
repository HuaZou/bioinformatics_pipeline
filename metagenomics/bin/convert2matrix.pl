#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($script, $tree, $real_dir, $out, $help, $version);
GetOptions(
    "s|script:s"    => \$script,
    "t|tree:s"      => \$tree, 
    "d|real_dir:s"  => \$real_dir,   
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

my $metaphlan_dir = "$real_dir/result/04.metaphlan/";
my $humann_pathabundance_dir = "$real_dir/result/03.humann/pathabundance/";
my $humann_pathcoverage_dir = "$real_dir/result/03.humann/pathcoverage/";
my $humann_genefamilies_dir = "$real_dir/result/03.humann/genefamilies/"; 
my $dir = "$real_dir/result/05.profile/"; 
system "mkdir -p $dir" unless(-d $dir);

# script
my $dir_script = "$real_dir/result/script/05.profile/"; 
system "mkdir -p $dir_script" unless(-d $dir_script);


open(OT, "> $out") or die "can't open $out\n";
my $bash = join("", $dir_script, "profile.sh");
open(OT2, "> $bash") or die "can't open $bash\n";
my $metaphlan_file = "$dir/all_merge_metaphlan.tsv";
my $humann_pathabundance_file = "$dir/all_merge_pathabundance.tsv";
my $humann_pathcoverage_file  = "$dir/all_merge_pathcoverage.tsv";
my $humann_genefamilies_file  = "$dir/all_merge_genefamilies.tsv";
print OT2 "merge_metaphlan_tables.py $metaphlan_dir/*metagenome.tsv | sed \'s/_metagenome//g\' > $metaphlan_file\n";
print OT2 "humann_join_tables --input $humann_pathabundance_dir --output $humann_pathabundance_file && sed -i \'s/_merge_Abundance//g\'  $humann_pathabundance_file\n";
print OT2 "humann_join_tables --input $humann_pathcoverage_dir --output $humann_pathcoverage_file && sed -i \'s/_merge_Coverage//g\'  $humann_pathcoverage_file\n";
print OT2 "humann_join_tables --input $humann_genefamilies_dir --output $humann_genefamilies_file && sed -i \'s/_merge_Abundance\-RPKs//g\'  $humann_genefamilies_file\n";

my $weighted_unifrac = "$dir/weighted_unifrac.tsv";
my $unweighted_unifrac = "$dir/unweighted_unifrac.tsv";
print OT2 "Rscript $script $metaphlan_file $tree $weighted_unifrac weighted\n";
print OT2 "Rscript $script $metaphlan_file $tree $unweighted_unifrac unweighted\n";

close(OT2);
print OT "sh $bash\n";
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -s <script> -t <tree> -d <real_dir> -o <out> 
options:
    -s|script   :[essential]. 
    -t|tree     :[essential]. 
    -d|real_dir :[essential].    
	-o|out      :[essential].

USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.0
    update:     20210106 - 20210106
    author:     zouhua1\@outlook.com
VERSION
};
