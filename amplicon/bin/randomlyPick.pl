#!/usr/bin/env perl
# Randomly split reads into multiple samll files
#-----------------------------------------------------------------------------
# Author : Chao Fang
# Email  : fangchao@genomics.cn
# Create : Nov 2018
#-----------------------------------------------------------------------------
# see detail below
use strict;
use Getopt::Long;

sub usage {
  my $msg = shift;
  $msg||="Randomly split reads into multiple samll files";
print <<USAGE;
$msg
  * Undetected barcodes will be randomly output.
  * Otherwize sequences will be randomly output using barcode as unit.
usage:
  $0 -i input -o out_prefix -n num -s seed
    -i  input file (fastq/fasta)
    -o  output prefix
    -n  number of split files
    -m  How many to output (default is equal to n)
    -s  seed of srand(), default is randomly
    -t  threads num for compress, default 1
    -v  verbose
    -h  show help info
USAGE
}
&usage && exit unless @ARGV;

my ($inf,$out,$num,$outNum,$seed,$thread,$verbose,$help);
GetOptions(
  "i=s" => \$inf,
  "o=s" => \$out,
  "n=i" => \$num,
  "m=i" => \$outNum,
  "s:i" => \$seed,
  "t:i" => \$thread,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;
&usage("[fatal] Essential input is missing") && exit unless defined $inf && $num && $out;

$thread||=1;
$outNum||=$num;
my $fmt = $1 if $inf =~ /.(f(|ast)[aq](|.gz))$/;
&verbose("[log] input format detected: $fmt\n");
if($fmt =~ /.gz$/){
  open INF, "pigz -p $thread -dc $inf|" or die $!;
}else{
  open INF, "< $inf" or die $!;
}
my $fmtR = ($fmt =~ /f(|ast)q/)?"fq":"fa";
my %FH;
for(my $i=0;$i<$outNum;$i++){
  my $oFile = sprintf("%s_%02d.%s",$out,$i,$fmt);
  if($fmt =~ /.gz$/){
    open $FH{$i}, "|pigz -p $thread > $oFile" or die $!;
  }else{
    open $FH{$i}, ">$oFile" or die $!;
  }
  &verbose("[log] open filehandle: $oFile\n");
}
# Main start
$seed = ($seed)?srand($seed):srand();
&verbose("[log] Use seed: $seed\n");

&verbose("\n[log] Start IO ...\n");
my ($chunk,$idSnap,$seedi,$readCount);
while(<INF>){
  my @bc = &getIdBarcode($_);
  $chunk = ($fmtR eq "fq")?$_.<INF>.<INF>.<INF>:$_.<INF>;

  if($bc[1] ne $idSnap || $bc[1]=~/0000/){
    $seedi = int(rand($num));
    $idSnap = $bc[1];
  }
  next if $seedi >= $outNum;
  $FH{$seedi} -> print($chunk);
  $readCount ++;
  if($readCount % 1000000 == 0){
    &verbose("[log] Working on $bc[1] ...\n");
  }
}
close INF;
close OUT;
# Main end

&verbose("done!\n");

exit;

sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}

sub getIdBarcode {
  my $id = shift;
  my @bc;
  if($id =~ /[@\/](\d+)_(\d+)_(\d+)\//){
    my @bcode= ($1,$2,$3);
    @bc = ($bcode[0],"$bcode[0]_$bcode[1]_$bcode[2]");
  }else{
    my $unknown="0000";
    @bc = ($unknown,"$unknown\_$unknown\_$unknown");
  }
  return(@bc);
}
