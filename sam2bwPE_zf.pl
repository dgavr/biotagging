#!/usr/local/bin/perl

use Pod::Usage;
use Getopt::Long;
use strict;

system ("module load ucsctools");



=head1 NAME

 sam2bwPE.pl 

=head1 SYNOPSIS

 sam2bwPE.pl -sam data.sam -build hg18 -name lane1

=head1 EXAMPLE

 sam2bwPE.pl -sam lane1.sam -build hg18 -name lane1 -desc GSM449526_lane1 -path /hts/data0/public/JBloggs/ -strict n
 
=head1 DESCRIPTION

 Takes sam file as input and produces a windowed Density track across the genome,
 ready to load into UCSC as a binary wig file.  Dependant on the UCSC Script wigToBigWig
 being installed on the system.  This script requires a file of chromosome sizes for each
 build to be used.  HG18 , HG19, Gal4 , MM9 and yeast  are installed at present.
 Both paired end and single end Sam files can be analysed and the stript will automatically
 recognise each as such.
 The output file is moved by the script to a public access folder and chmoded.
 A path to this public folder can be supplied using -path, or if left blank the script will
 attempt to predict the path from your user name.

=head1 OPTIONS

   -sam 	input file (SAM format)
   -build 	genome build used during alignment  (default: mm9, alternatively use hg18, hg19, or yeast)
   -window 	window size  (default: 300)
   -inc 	window increment  (default: 30)
   -name 	track name  (default: track_name)
   -desc 	track description  (default: track_desc)
   -path 	Path to public folder to move bwig file to.
   -strict	Only use reads mapped as a proper pair, if set to n , the reads that are mapped but not in a pair will be included. (default: y)
   -clip	unless set to 'n', wigToBigWig will issue warning messages rather than dying if wig file contains items off end of chromosome.

=head1 AUTHOR

 Written by Jim Hughes 2010

=cut



my $help=0;
my $man=0;	


my $window_size;
my $window_incr;
my $track_name;
my $track_desc;
my $build;
my $strict;
my %chrlength;
my $clip;

my $wigToBigWig = '/package/wigtool/default/bin/wigToBigWig';

my $username;

my $cmd = "/usr/bin/whoami";   
open (PIPE, "$cmd |") || die "$!";
while (<PIPE>)
{
	$username = $_;  
	chomp $username; 
}
close PIPE;

my $path = "/public/$username";


&GetOptions (
	      "sam=s"=>\my $samfile,
	     "window=i"=> \$window_size,
	     "name=s"=>\$track_name,
	     "desc=s"=>\$track_desc,
	     "inc=i" => \$window_incr,
	     "build=s" => \$build,
	     "strict=s" => \$strict,
	     "path=s" => \ $path,
	     "clip=s" => \ $clip,
	     'h|help'=>\$help,
	     'man'=>\$man);


pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($samfile); 
    
# set some defaults:
unless($window_size) {$window_size = 300;}
unless($window_incr) {$window_incr = 30;}
unless($track_name) {$track_name = 'track_name';}
unless($track_desc) {$track_desc = ' track_desc';}
unless($build) {$build = 'danRer10';}
unless($strict) {$strict = 'y';}
unless($clip) {$clip = 'y';}

my $chr_lengths_file = '/hts/data1/dariag/genomes/danRer10/Zv10.chrom.sizes';

if ($build eq 'danRer10') {$chr_lengths_file = '/hts/data1/dariag/genomes/danRer10/Zv10.chrom.sizes';}
elsif ($build eq 'hg18') {$chr_lengths_file = '/hts/data0/config/bigwig/hg18_sizes.txt';}
elsif ($build eq 'yeast') {$chr_lengths_file = '/hts/data0/config/bigwig/yeast_sizes.txt';}
elsif ($build eq 'hg19') {$chr_lengths_file = '/hts/data0/config/bigwig/hg19_sizes.txt';}
elsif ($build eq 'mm10') {$chr_lengths_file = '/hts/data0/config/bigwig/mm10_sizes.txt';}

open (SIZES, $chr_lengths_file);
	
while (<SIZES>)
{
	chomp;
	my ($gchr, $gsize)= split(/\t+/);
	print "$gchr $gsize\n";

	$gchr =~ s/chr//gi;
	unless ($gchr =~ /M/gi)  { $chrlength{$gchr}=$gsize; }		
}	
close SIZES;
	





my $pos_fix = int($window_size/2);




open (OUTPUT, ">$track_name.wig");
open (OUTPUT1, ">$track_name.url.txt");


#reads chr lengths from sam file

#open(INFO, $samfile) || die $!;

#while (<INFO>) {
	#if (/\@SQ(\s+)SN:(chr)?(\S+)(\s+)(\S+)?(\s+)?LN:(\d+)/){
	#	&lengths;
	#} 
#}
#close INFO;



open(INFO1, $samfile) || die $!;

my %bins;
my %read_posns;
my $bintotal;

while (<INFO1>) {
	&HASHIT;
}

close INFO1;

my $binsnumbers;

foreach my $dChr(sort by_number keys %bins) {
	print OUTPUT "variableStep  chrom=chr$dChr span=$window_incr\n";
	foreach my $dcoor(sort by_number keys %{$bins{$dChr}}) {
		
		my $dcount = $bins{$dChr}{$dcoor};
		my $chromolength = $chrlength{$dChr};

		my $adjustcoor = $dcoor + $pos_fix;   ####move coordinate from start of window to middle of window

		unless ($adjustcoor > $chromolength){
			print OUTPUT "$adjustcoor\t$dcount\n";
			$binsnumbers++;
			$bintotal += $dcount;
			}	
}}


print OUTPUT1 "track type=bigWig name=\"$track_name\" description=\"$track_desc\" bigDataUrl=http://sara.molbiol.ox.ac.uk/public/$username/$track_name.bw\n";	


#run wigToBigWig (edited by SMcG, 22Oct2013):
my $chr_size_file;
if ($build =~ /danRer10/i) {$chr_size_file = '/hts/data1/dariag/genomes/danRer10/Zv10.chrom.sizes';}
elsif ($build =~ /hg18/i) {$chr_size_file = '/hts/data0/config/bigwig/hg18_sizes.txt';}
elsif ($build =~ /hg19/i) {$chr_size_file = '/hts/data0/config/bigwig/hg19_sizes.txt';}
elsif ($build =~ /mm9/i) {$chr_size_file = '/hts/data0/config/bigwig/mm9_sizes.txt';}


my $wig_to_bigwig_cmd = "$wigToBigWig -clip $track_name.wig $chr_size_file $track_name.bw";
if ($clip =~ /n/gi) {$wig_to_bigwig_cmd = "$wigToBigWig $track_name.wig $chr_size_file $track_name.bw";}

system ($wig_to_bigwig_cmd) == 0 or die "couldn't bigwig $build files (using $wig_to_bigwig_cmd): $!\n";


#if (-e $path){
#	system ("mv $track_name.bw $path/") == 0 or die "couldn't move files\n";		
#	system ("chmod 755 $path/$track_name.bw") == 0 or die "couldn't chmod files\n";
#} else {
#	print "Path to public folder not valid, manually move bigwig file to a public folder!\n";
#}




my $avwindowSig = $bintotal/$binsnumbers;

print "Average window signal is $avwindowSig\n";


##########subs#########


#sub lengths {
    
	#if (/\@SQ(\s+)SN:(chr)?(\S+)(\s+)(\S+)?(\s+)?LN:(\d+)/){
    
	#  my $snchr = $3;
	#  my $snlenght = $7;
	 # print "$snchr $snlenght\n";
    
    
	#unless (($snchr =~ /M/gi) || ($snchr =~ /_/g)){   
	#	 $chrlength{$snchr}=$snlenght;
	#}
#}

#if ((/NM:/g)  && (/MD:/g)){
#last;
#}
#}





sub HASHIT {

	if (($_ =~ /^@/) || ($_ =~ /^\#/)){ next; }	

	my ($readsam, $bitwisesam, $chrsam, $startsamm)=split(/\s+/);
	if ($chrsam =~ /_/g){next;}
	$chrsam =~ s/chr//;
	
	unless(exists $chrlength{$chrsam}){next;}


unless ($bitwisesam & 0x0001)  {										#####     unless the read is paired end
	if (($bitwisesam == 0) || ($bitwisesam == 16)){								#####     Bitwise for forward or reverse mapping for single end
		
		my $int = int($startsamm / $window_size);
		my $start_bin = ($int * $window_size);
		my $diff = $startsamm - $start_bin;
		my $incr = (int(($window_size - $diff) / $window_incr) * $window_incr);
		$start_bin -= $incr;
	
		for (my $bin=$start_bin; $bin<($start_bin+$window_size); $bin+=$window_incr){
			unless (($chrsam =~ /M|m/)||($bin <= 0)) {
			$bins{$chrsam}{$bin} ++;
		}

		}
	}
}


if ($bitwisesam & 0x0001)  {                                                                                    #####     the read is paired in sequencing, (but not necesarily mapped in a pair)
	
	if ($strict eq 'y') {
		unless ($bitwisesam & 0x0004){										#####     the query sequence itself is unmapped
			if ($bitwisesam & 0x0002){									#####	  read is mapped in a proper pair
		
				my $int = int($startsamm / $window_size);
				my $start_bin = ($int * $window_size);
				my $diff = $startsamm - $start_bin;
				my $incr = (int(($window_size - $diff) / $window_incr) * $window_incr);
				$start_bin -= $incr;
	
				for (my $bin=$start_bin; $bin<($start_bin+$window_size); $bin+=$window_incr){
					unless (($chrsam =~ /M|m/)||($bin <= 0)) {
					$bins{$chrsam}{$bin} ++;
				}

			}	
		
		
	}	
	}
	}else {
		
		unless ($bitwisesam & 0x0004){	
		
			my $int = int($startsamm / $window_size);
			my $start_bin = ($int * $window_size);
			my $diff = $startsamm - $start_bin;
			my $incr = (int(($window_size - $diff) / $window_incr) * $window_incr);
			$start_bin -= $incr;
	
			for (my $bin=$start_bin; $bin<($start_bin+$window_size); $bin+=$window_incr){
			
				unless (($chrsam =~ /M|m/)||($bin <= 0)) {
				$bins{$chrsam}{$bin} ++;
			}

			}	
		}
	}
}
}

sub by_number {
	($a <=> $b);
	}
