#!/usr/bin/env perl
# Modified by August Woerner
# This was originally written by Amy Williams and
# released as part of ibis.
#
# now modified to add cM positions to a VCF

# The HapMap genetic map is available from
#   ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/
#
# NOTE: the names of the chromosomes in the PLINK bim/map file and genetic map
# need to match. the script also handles differences in the presence/absence of
# a 'chr' prefix, but will fail otherwise if no match is found.
#
# NOTE: it is imperative that the physical positions of the two files are from
# the same genome build. For example, both could have physical positions from
# build 36 or both could have physical positions from build 37, but avoid
# mixing physical positions from different builds.
#
# Prints the updated bim/map file to STDOUT

use strict;
use warnings;
use Getopt::Long;

my $use_cM = 1;
my $no_header = 0;
my $beagleMap=0;

GetOptions( "-cm" => \$use_cM,
	    "-beagle" => \$beagleMap,
	    "-noheader" => \$no_header);

if (@ARGV < 1) {
  &print_usage();
}

# an ordered list of 2-long arrays. Two values are physical and genetic position
my %map;

foreach my $map_file (@ARGV) {
  print STDERR "Reading map file $map_file... ";
  open MAP, "$map_file" or die "\nCouldn't open $map_file: $!\n";

  if (!$no_header) {
    # omit header
    <MAP>;
  }

  my $last_chr = '';
  my $last_pos = -1;

  while (my $line = <MAP>) {
    $line =~ s/^\s+//;
    my @fields = split /\s+/, $line;
    if (@fields != 4) {
      my $num_cols = @fields;
      die "\nRead line with $num_cols columns from map but should always have 4 columns\n";
    }

    my $chr = $fields[0];
    my $phys_pos = $fields[1];
    my $genet_pos = $fields[3];

    if ($beagleMap) {
	$phys_pos = $fields[3];
	$genet_pos = $fields[2];
    }
    
    my $chr_key = $chr;
    $chr_key =~ s/^chr//; # remove 'chr' prefix, if present

    if ($chr ne $last_chr) {
      if (exists $map{ $chr_key }) {
	die "\nERROR: ensure only one file contains data for any chromosome and that\nchromosomes are not interspersed in a file";
      }
      $last_chr = $chr;
      $last_pos = -1;
      $map{ $chr_key } = [];
    }

    if ($phys_pos < $last_pos) {
      die "\nERROR: physical positions not in ascending order on chromosome $chr\nposition $phys_pos is after $last_pos\n";
    }
    $last_pos = $phys_pos;

    if ($use_cM) {
      # $genet_pos already in cM
      push @{ $map{ $chr_key } }, [ $phys_pos, $genet_pos ];
    }
    else {
      # convert to Morgans:
      push @{ $map{ $chr_key } }, [ $phys_pos, $genet_pos / 100 ];
    }
  }

  close MAP;
  print STDERR "done\n";
}


my $cur_chr = '';
my $cur_index = 0;
my $last_pos = -1;
my $have_warned = 0;

while ($_ = <STDIN>) {

    if (m/^#/) {

	if (m/^##/) {
	    print $_;
	    next;
	}
	print '##INFO=<ID=GPOS,Number=A,Type=Float,Description="cM position in genetic map">' , "\n";
	print $_;
	next;
    }
    
  my @fields = split /\s+/;

    
    my $chr = $fields[0];
    $chr =~ s/^chr//; # remove 'chr' prefix, if present
    my $phys_pos = $fields[1];

  if ($chr ne $cur_chr) {
    if (not exists $map{ $chr }) {
      #die "\nERROR: no chromosome $chr in genetic map, but is present in bim/map file";
    }
    $cur_chr = $chr;
    $cur_index = 0;
    $last_pos = 0;
    print $_;
    next;
  }

  if ($phys_pos < $last_pos) {
    die "\nError: physical positions not in ascending order in bim/map file\nAt SNP $fields[1] (chromosome $chr)\n";
  }
  $last_pos = $phys_pos;

  my $genet_pos;
  if ($phys_pos == 0) {
    $genet_pos = 0.0;
  }
  else {
    # find the index immediately after $phys_pos:
    for( ; $cur_index < scalar @{ $map{ $cur_chr } } &&
	     $map{ $cur_chr }[ $cur_index ][0] < $phys_pos; $cur_index++) { }

    if ($cur_index >= scalar @{ $map{ $cur_chr } }) {
      if (!$have_warned) {
	print STDERR "\n";
	$have_warned = 1;
      }
      #print STDERR "Warning: position $phys_pos after chrom $cur_chr end; map and physical pos set to 0\n";
      $genet_pos = 0.0;
      #$fields[3] = 0;
    }
    elsif ($map{ $cur_chr }[ $cur_index ][0] == $phys_pos) {
      $genet_pos = $map{ $cur_chr }[ $cur_index ][1];
    }
    else {
      if ($cur_index == 0) {
	if (!$have_warned) {
	  print STDERR "\n";
	  $have_warned = 1;
	}
	print STDERR "Warning: position $phys_pos before chrom $cur_chr start; genetic map pos set to 0\n";
	$genet_pos = 0.0;
#	$fields[3] = 0;
      }
      else {
	# linear interpolation:
	my $rel_distance = ($phys_pos - $map{ $cur_chr }[ $cur_index - 1 ][0]) /
		    ($map{ $cur_chr }[ $cur_index ][0] -
					$map{ $cur_chr }[ $cur_index - 1 ][0]);
        die "\nLinear interploation bug" if ($rel_distance > 1.0 || $rel_distance < 0.0);
	$genet_pos = $map{ $cur_chr }[ $cur_index - 1][1] +
		    $rel_distance *
		    ($map{ $cur_chr }[ $cur_index ][1] -
					$map{ $cur_chr }[ $cur_index - 1 ][1]);
      }
    }
  }

  if (length($fields[7])== 0 || $fields[7] eq '.') {
      $fields[7]='GPOS='.$genet_pos;
  } else {
      $fields[7].=';GPOS='.$genet_pos;
  }

    print join("\t", @fields) , "\n";
 # print "$fields[0]\t$fields[1]\t$genet_pos\t$fields[3]";
  #if (@fields == 6) {
   # print "\t$fields[4] $fields[5]";
  #}
  #print "\n";
}



sub print_usage() {
  print STDERR "Usage:\n";
  print STDERR "cat vcfFile |  $0 genetic map file(s)\n";
  print STDERR "   Accepts any number of genetic map files, so long as data for one chromosome\n";
  print STDERR "   is in only one file\n";
  print STDERR "\n";
  print STDERR "Example usage:\n";
  print STDERR "bcftools view myBcf.bcf | $0 /path/to/map/genetic_map_GRCh37_chr{1..22}.txt > new.vcf\n";
  exit;
}
