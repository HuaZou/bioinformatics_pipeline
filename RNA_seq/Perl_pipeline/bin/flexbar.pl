#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($file, $prefix, $adapter, $quality, $uncalled, $length, $threads, $out, $help, $version);
GetOptions(
    "f|file:s"      =>  \$file,
    "p|prefix:s"    =>  \$prefix,
    "a|adapter:s"   =>  \$adapter,
    "q|quality:s"   =>  \$quality,
    "u|uncalled:s"  =>  \$uncalled,
    "l|length:s"    =>  \$length,
    "t|threads:s"   =>  \$threads,            
    "o|out:s"       =>  \$out,
    "h|help:s"      =>  \$help,
    "v|version"     =>  \$version
);
&usage if(!defined $prefix);

# output
my $dir_trim = "$out/flexbar/"; 
system "mkdir -p $dir_trim" unless(-d $dir_trim);

# script
my $dir_script = "$out/script/flexbar/"; 
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
open(OT, "> $prefix") or die "can't open $prefix\n";
foreach my $key (keys %file_name){
    if (${$file_name{$key}}[0] =~ /[r|R]1/){
        $fq1 = ${$file_name{$key}}[0];
    }else{
        $fq2 = ${$file_name{$key}}[0];
    }

    if (${$file_name{$key}}[1] =~ /[r|R]2/){
        $fq2 = ${$file_name{$key}}[1];
    }else{
        $fq1 = ${$file_name{$key}}[1];
    }
    my $bash = join("", $dir_script, $key, ".flexbar.sh");
    open(OT2, "> $bash") or die "can't open $bash\n";
    print OT2 "flexbar  -r $fq1 -p $fq2 --qtrim-threshold $quality -ao 5 --min-read-length $length --threads $threads --target $dir_trim/$key -a $adapter  --zip-output GZ --max-uncalled $uncalled\n";       
    close(OT2);
    
    print OT "sh $bash\n";        
}
close(OT);

sub usage{
	print <<USAGE;
usage:
	perl $0 -f <file> -p <prefix> -o <out>
options:
	-f|file     :[essential]. fqpath (sampleID/laneid/fqpath)
    -p|prefix   :[essential]. the name of this processure   
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
