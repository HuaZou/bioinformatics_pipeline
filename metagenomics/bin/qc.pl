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

# output
my $dir_qc = "$real_dir/result/00.quality/fastqc"; 
system "mkdir -p $dir_qc" unless(-d $dir_qc);

# script
my $dir_script = "$real_dir/result/script/00.quality/"; 
system "mkdir -p $dir_script" unless(-d $dir_script);

my @array_name;
open(IN, $file) or die "can't open $file\n";
open(OT, "> $out") or die "can't open $out\n";
<IN>;
while(<IN>){
    chomp;
    my @tmp = split("\t", $_);
    if(-e $tmp[2]){
        #print OT "fastqc -o $dir_qc --noextract $tmp[2]\n";
        my $bash = join("", $dir_script, $tmp[1], ".fastqc.sh");
        open(OT2, "> $bash") or die "can't open $bash\n";
        print OT2 "fastqc -o $dir_qc --noextract $tmp[2]\n";
        close(OT2);

        print OT "sh $bash\n";
    }
}
close(IN);

my $dir_mc = "$real_dir/result/00.quality/multiqc"; 
system "mkdir -p $dir_mc" unless(-d $dir_mc);
print OT "multiqc $dir_qc --outdir $dir_mc\n";
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
