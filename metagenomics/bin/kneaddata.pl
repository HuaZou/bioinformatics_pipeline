#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $real_dir, $out, $adapter, $help, $version);
GetOptions(
    "f|file:s"   =>  \$file,
    "d|real_dir:s"  => \$real_dir,
    "a|adapt:s"  =>  \$adapter,  
    "o|out:s"    =>  \$out,
    "h|help:s"   =>  \$help,
    "v|version"  =>  \$version
);
&usage if(!defined $out);

my $dir = "$real_dir/result/01.kneaddata"; 
system "mkdir -p $dir" unless(-d $dir);

# script
my $dir_script = "$real_dir/result/script/01.kneaddata/"; 
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

my ($fq1, $fq2);
open(OT, "> $out") or die "can't open $out\n";
foreach my $key (keys %file_name){
    if (${$file_name{$key}}[0] =~ /R1/){
        $fq1 = ${$file_name{$key}}[0];
    }else{
        $fq2 = ${$file_name{$key}}[0];
    }

    if (${$file_name{$key}}[1] =~ /R2/){
        $fq2 = ${$file_name{$key}}[1];
    }else{
        $fq1 = ${$file_name{$key}}[1];
    }
    my $trim_opt = "ILLUMINACLIP:$adapter:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50";
    #if(-e $fq1 && -e $fq2){
        #print OT "kneaddata -i $fq1 -i $fq2 --output-prefix $key -o $dir -v -t 5 --remove-intermediate-output  --trimmomatic /data/share/anaconda3/share/trimmomatic/ --trimmomatic-options \'$trim_opt\'  --bowtie2-options \'--very-sensitive --dovetail\' -db /data/share/database/kneaddata_database/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens\n";       
        my $bash = join("", $dir_script, $key, ".kneaddata.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";
        print OT2 "kneaddata -i $fq1 -i $fq2 --output-prefix $key -o $dir -v -t 5 --remove-intermediate-output  --trimmomatic /data/share/anaconda3/share/trimmomatic/ --trimmomatic-options \'$trim_opt\'  --bowtie2-options \'--very-sensitive --dovetail\' -db /data/share/database/kneaddata_database/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens\n";       
        close(OT2);

        print OT "sh $bash\n";        
    #}
}

print OT "kneaddata_read_count_table --input $dir --output $dir/01kneaddata_sum.tsv\n";
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -d <real_dir> -o <out> -a <adapter>
options:
	-f|file		:[essential].
    -d|real_dir :[essential].  
	-o|out      :[essential].
    -a|adapt    :[essential].
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
