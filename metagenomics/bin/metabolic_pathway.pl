#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $real_dir, $out, $db_dir, $help, $version);
GetOptions(
    "f|file:s"      => \$file, 
    "d|real_dir:s"  => \$real_dir,   
    "db|db_dir:s"   =>  \$db_dir,
    "o|out:s"   =>  \$out,
    "h|help:s"  =>  \$help,
    "v|version" =>  \$version
);
&usage if(!defined $out);

# output dir
my $metabolic_dir = "$real_dir/result/06.metabolic_pathway/"; 
my $MetaCyc_Reactions_dir = "$metabolic_dir/MetaCyc_Reactions/";
my $KEGG_Orthogroups_dir = "$metabolic_dir/KEGG_Orthogroups/";
my $Pfam_domains_dir = "$metabolic_dir/genefamilies/"; 
my $Level4_enzyme_dir = "$metabolic_dir/Level4_enzyme/";
my $EggNOG_COGs_dir = "$metabolic_dir/EggNOG_COGs/";
my $Gene_Ontology_dir = "$metabolic_dir/Gene_Ontology/"; 
#my $Informative_GO_dir = "$metabolic_dir/Informative_GO/"; 

system "mkdir -p $metabolic_dir" unless(-d $metabolic_dir);
system "mkdir -p $MetaCyc_Reactions_dir" unless(-d $MetaCyc_Reactions_dir);
system "mkdir -p $KEGG_Orthogroups_dir" unless(-d $KEGG_Orthogroups_dir);
system "mkdir -p $Pfam_domains_dir" unless(-d $Pfam_domains_dir);
system "mkdir -p $Level4_enzyme_dir" unless(-d $Level4_enzyme_dir);
system "mkdir -p $EggNOG_COGs_dir" unless(-d $EggNOG_COGs_dir);
system "mkdir -p $Gene_Ontology_dir" unless(-d $Gene_Ontology_dir);
#system "mkdir -p $Informative_GO_dir" unless(-d $Informative_GO_dir);

# reference database
my $MetaCyc_Reactions_db = "$db_dir/map_ec_name.txt.gz";
my $KEGG_Orthogroups_db = "$db_dir/map_ko_uniref90.txt.gz";
my $Pfam_domains_db = "$db_dir/map_pfam_uniref90.txt.gz"; 
my $Level4_enzyme_db = "$db_dir/map_level4ec_uniref90.txt.gz";
my $EggNOG_COGs_db = "$db_dir/map_eggnog_uniref90.txt.gz";
my $Gene_Ontology_db = "$db_dir/map_go_uniref90.txt.gz"; 
#my $Informative_GO_db = "$db_dir/Informative_GO/"; 

# script
my $dir_script = "$real_dir/result/script/06.metabolic_pathway/"; 
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

my $genefamilies_file;
open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    $genefamilies_file = join("", $real_dir, "/result/03.humann/genefamilies/", $key, "_genefamilies.tsv");
    my $bash = join("", $dir_script, $key, ".metabolic_pathway.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $MetaCyc_Reactions_db -o $MetaCyc_Reactions_dir/$key\_MetaCyc.tsv\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $KEGG_Orthogroups_db -o $KEGG_Orthogroups_dir/$key\_KEGG_Orthogroups.tsv\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $Pfam_domains_db -o $Pfam_domains_dir/$key\_Pfam_domains.tsv\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $Level4_enzyme_db -o $Level4_enzyme_dir/$key\_Level4_enzyme.tsv\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $EggNOG_COGs_db -o $EggNOG_COGs_dir/$key\_EggNOG_COGs.tsv\n";
    print OT2 "humann_regroup_table -i $genefamilies_file -c $Gene_Ontology_db -o $Gene_Ontology_dir/$key\_Gene_Ontology.tsv\n";
    close(OT2);

    print OT "sh $bash\n";
}

print OT "humann_join_tables --input $MetaCyc_Reactions_dir --output $metabolic_dir/all_merge_MetaCyc_Reactions.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_MetaCyc_Reactions.tsv\n";
print OT "humann_join_tables --input $KEGG_Orthogroups_dir --output $metabolic_dir/all_merge_KEGG_Orthogroups.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_KEGG_Orthogroups.tsv\n";
print OT "humann_join_tables --input $Pfam_domains_dir --output $metabolic_dir/all_merge_Pfam_domains.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_Pfam_domains.tsv\n";
print OT "humann_join_tables --input $Level4_enzyme_dir --output $metabolic_dir/all_merge_Level4_enzyme.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_Level4_enzyme.tsv\n";
print OT "humann_join_tables --input $EggNOG_COGs_dir --output $metabolic_dir/all_merge_EggNOG_COGs.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_EggNOG_COGs.tsv\n";
print OT "humann_join_tables --input $Gene_Ontology_dir --output $metabolic_dir/all_merge_Gene_Ontology.tsv && sed -i \'s/_merge_Abundance-RPKs//g\' $metabolic_dir/all_merge_Gene_Ontology.tsv\n";

print OT "humann_renorm_table --input $metabolic_dir/all_merge_MetaCyc_Reactions.tsv --units relab --output $metabolic_dir/all_merge_MetaCyc_Reactions_relab.tsv\n";
print OT "humann_renorm_table --input $metabolic_dir/all_merge_KEGG_Orthogroups.tsv --units relab --output $metabolic_dir/all_merge_KEGG_Orthogroups_relab.tsv\n";
print OT "humann_renorm_table --input $metabolic_dir/all_merge_Pfam_domains.tsv --units relab --output $metabolic_dir/all_merge_Pfam_domains_relab.tsv\n";
print OT "humann_renorm_table --input $metabolic_dir/all_merge_Level4_enzyme.tsv --units relab --output $metabolic_dir/all_merge_Level4_enzyme_relab.tsv\n";
print OT "humann_renorm_table --input $metabolic_dir/all_merge_EggNOG_COGs.tsv --units relab --output $metabolic_dir/all_merge_EggNOG_COGs_relab.tsv\n";
print OT "humann_renorm_table --input $metabolic_dir/all_merge_Gene_Ontology.tsv --units relab --output $metabolic_dir/all_merge_Gene_Ontology_relab.tsv\n";

close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <script> -d <real_dir> -db <database> -o <out> 
options:
    -f|file     :[essential].
    -d|real_dir :[essential].  
    -db|db_dir  :[essential].
    -d|real_dir :[essential].    
    -o|out      :[essential].

USAGE
    exit;
};

sub version {
    print <<VERSION;
    version:    v1.0
    update:     20210111 - 20210111
    author:     zouhua1\@outlook.com
VERSION
};
