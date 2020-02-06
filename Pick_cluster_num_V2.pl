#!/usr/local/bin/perl -w
# 
# Copyright (c) 
# Writer:         wangxy <wxyang1988@126.com>
# Program Date:   2011.
# Modifier:       wangxy <wxyang1988@126.com>
# Last Modified:  2011.
##########################################################

my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#########################################################

my %opts;
GetOptions(\%opts,"id=s","o=s","h" ); 

#&help()if(defined $opts{h});
if(!defined($opts{id}) || !defined($opts{o}) ||defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-id     Input folder                    must be given

		-o      Output file                  must be given

		-h    Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
my $programe_dir=basename($0);
my $ind=$opts{id};
my $out=$opts{o};

#S	0	573	*	*	*	*	*	T12_45085 G5JZO0C01BPBTT orig_bc=ACGAGTGCGT new_bc=ACGAGTGCGT bc_diffs=0 1..573	*

my @file=glob"$ind/*.uc";
open OUT,">$out" || die "Can't create $out,$!\n";
for (my $j=0;$j<@file ;$j++) 
{
	my %cluster;
	my $file_name=(split/\//,$file[$j])[-1];
	my $seq_name=(split/\.uc/,$file_name)[0];
	open IN,"$file[$j]" || die "Can't open $file[$j],$!\n" ;
	while (<IN>) 
	{
		chomp;
		my $cluster_id;
		next if (/^\#/);
		my @a=split/\t/,$_;
		if ($a[0] eq "S") 
		{
			$cluster_id=(split/\s+/,$a[8])[0];
			$cluster{$a[1]}=$cluster_id;
		}
		elsif ($a[0] eq "H") 
		{
#			for (my $i=8;$i<@a ;$i++) 
#			{
				$cluster_id=(split/\s+/,$a[8])[0];
#				next if ($cluster_id eq $cluster{$a[1]});
				$cluster{$a[1]}.="\t$cluster_id";
#			}
		}
	}
	close (IN);
	foreach my $id (sort {$a<=>$b} keys %cluster) 
	{
		print OUT $seq_name,"_",$id,"\t",$cluster{$id},"\n";
	}
}

close (OUT);
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
