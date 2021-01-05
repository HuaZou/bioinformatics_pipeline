#!/usr/bin/perl -w
use strict;
use warnings;
################################################################################
unless(3<=@ARGV) {
    &usage;
    exit;
}
################################################################################
my($list_f,$row,$out,$alpha,$beta) = @ARGV;
my(@order,%list,@info,$i,%class,%cover,%Count,%sum_abun);
################################################################################
open IN,"<$list_f" || die "read $list_f $!\n";
while(<IN>) {
    chomp;
    @info=split /\t|\s+/;
    my $name =$1;
    $list{$info[0]}=$info[1];
	push(@order,$info[0]);
}
close IN;
################################################################################
for($i=0;$i<@order;++$i) {
	next if not defined $list{$order[$i]};
	my $openMethod = ($list{$order[$i]} =~ /gz$/)?"gzip -dc $list{$order[$i]} |":"$list{$order[$i]}";
	open IN,$openMethod  or die "$!\n";
    while(<IN>) {
        chomp;
        @info=split /\t/;
	if ($.==1) {
		next unless &checkNumber($info[$row]);
	}
	if($beta){
		$class{$info[0]}{$i} = "\t".$info[$row];
	}else{
	        $class{$info[0]} .= "\t".$info[$row];
	}
	$cover{$info[0]} ||= 0;
	$sum_abun{$info[0]} ||=0;
	$Count{$i} ||= 0;
	if ($info[$row] > 0){ $cover{$info[0]} ++ ; $Count{$i} ++ ;$sum_abun{$info[0]} += $info[$row] };
    }
    close IN;
}
################################################################################
open OT,">$out.tsv" or die "write $out $!\n";
open CrT,">$out.cover" || die $!;
open CtT,">$out.count" || die $!;
for($i=0;$i<@order;++$i) {
	next if not defined $list{$order[$i]};
    print OT "\t",$order[$i];
	print CtT "$order[$i]\t$Count{$i}\n";
}
print OT "\n";
my @index;
if($alpha){@index = sort keys %class}else{@index = sort {$a<=>$b} keys %class}
foreach my $index(@index){
#foreach my $index (sort {$a<=>$b} keys %class) {
#	if ( not defined $cover{$index} || $cover{$index}==0 ){
#		$cover{$index}=1 unless $cover{$index} > 0;
#	}
	my $avg = ($sum_abun{$index} >0 )?($sum_abun{$index} / $cover{$index}):$sum_abun{$index};
	if($beta){
		print OT $index;
		for($i=0;$i<@order;++$i) {
			$class{$index}{$i} ||= "\t0";
			print OT $class{$index}{$i};
		}
		print OT "\n";
	}else{
		print OT $index,$class{$index},"\n";
	}
	print CrT "$index\t$cover{$index}\t$avg\n";
}
close OT;
close CrT;
################################################################################
sub checkNumber{    
    #return shift =~ /^[+\-]?([1-9]\d*|0)(\.\d+)?([eE][+\-]?([0-9]\d*|0)(\.\d+)?)?$/;
	return shift =~ /^[+\-]/;
}
sub usage {
    print STDERR<<USAGE;
Description
	This programme is to combine profiling table
Usage:
	perl $0 [file.list] [row] [outfile prefix]
		file.list needs 2 columns with one of the sample ID and following the path of the abundance file
		Row must be in 1(pairs),2(base_abundance),3(reads_abundance),4(depth_abundance)
USAGE
}
