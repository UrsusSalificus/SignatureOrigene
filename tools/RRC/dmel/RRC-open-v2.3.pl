#!/usr/bin/perl

##############################################################################
# The RRC tool is a free program: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright © 2009, 2010, 2011, 2012, 2013 Anna-Sophie Fiston-Lavier 
##############################################################################



use strict;
use warnings;
use Getopt::Long;

# Coordinates of Rate no NULL (in Mb):
# X  => [1.22-21.21]
# 2L => [0.53-18.87]
# 2R => [1.87-20.86]
# 3L => [0.75-19.02]
# 3R => [0.75-19.02]

my %dicoX;
my %dico2L;
my %dico2R;
my %dico3L;
my %dico3R;


&startup;

sub startup { 
    my $RELEASE="5.36";
    my $mapfile;
    my $help;
    my @order_chromosomes = ("2L","2R","3L","3R","X", "4");
    
    
    my %bp_limits_3 = (
	"2L" => "22,217,931",
    	"2R" => "20,302,755",
    	"3L" => "23,352,213",
    	"3R" => "27,890,790",
    	"X" => "21,780,003",
    	"4" => "1,237,870",
	);
    my %bp_limits_4 = (
	"2L" => "22,407,834",
    	"2R" => "20,736,785",
    	"3L" => "23,771,897",
    	"3R" => "27,905,053",
    	"X" => "22,224,390",
    	"4" => "1,281,640",
	);
    my %bp_limits_5 = (
	"2L" => "23,011,544",
    	"2R" => "21,146,708",
    	"3L" => "24,543,557",
    	"3R" => "27,905,053",
    	"X" => "22,422,827",
    	"4" => "1,351,857",
	);
    
    my %avg_rates_3 = (
	"2L" => 2.27,
	"2R" => 2.51,
	"3L" => 1.95,
	"3R" => 2.03,
	"X" => 3.06,
	"4" => 0.00,
	"Whole genome" => 2.31,
	);
    my %max_rates_3 = (
	"2L" => 4.04,
	"2R" => 3.81,
	"3L" => 3.40,
	"3R" => 3.28,
	"X" => 4.24,
	"4" => 0.00,
	"Whole genome" => 4.24,
	);
    my %avg_rates_4 = (
	"2L" => 2.22,
	"2R" => 2.43,
	"3L" => 1.87,
	"3R" => 1.99,
	"X" => 3.02,
	"4" => 0.00,
	"Whole genome" => 2.26,
	);
    my %max_rates_4 = (
	"2L" => 3.98,
	"2R" => 3.76,
	"3L" => 3.33,
	"3R" => 3.21,
	"X" => 4.25,
	"4" => 0.00,
	"Whole genome" => 4.25,
	);
    my %avg_rates_5 = (
	"2L" => 2.14,
	"2R" => 2.77,
	"3L" => 2.20,
	"3R" => 1.97,
	"X" => 3.09,
	"4" => 0.00,
	"Whole genome" => 2.46,
	);
    my %max_rates_5 = (
	"2L" => 4.01,
	"2R" => 3.78, 
	"3L" => 3.45,
	"3R" => 3.21,
	"X" => 4.22, 
	"4" => 0.00,
	"Whole genome" => 4.22,
	);
   my %avg_rates_5comeron = (
	"2L" => 2.39,
	"2R" => 2.66,
	"3L" => 1.79,
	"3R" => 1.96,
	"X" => 2.95,
	"4" => 0.00,
       "Whole genome" => 2.32,
       );
    my %max_rates_5comeron = (
	"2L" => 10.20,
	"2R" => 8.89, 
	"3L" => 7.82,
	"3R" => 14.80,
	"X" => 14.47, 
	"4" => 0.00,
	"Whole genome" => 14.80,
	);
    

    GetOptions ('R=s' => \$RELEASE,
		'M=s' => \$mapfile,
		'h|help' => \$help);
     
     if ($help) {
	&help();
	die "\n";
    }

    if($RELEASE) {
	print "Release: $RELEASE\n";
    }
    unless ($mapfile) {
	&help();
	die " *** Must specifiy a map file as input ***\n\n";
    }
    
    my ($pos,$rec);
    open(OX, "Comeron_tables/Comeron_100kb_chrX.txt"); 
    while(<OX>){
	chomp;
	($pos,$rec)=split();
	$dicoX{$pos}=$rec;
    }
    close OX;
    
    open(O2L, "Comeron_tables\/Comeron_100kb_chr2L.txt"); 
    while(<O2L>){
	chomp;
	($pos,$rec)=split();
	$dico2L{$pos}=$rec;
    }
    close O2L;
    
    open(O2R, "Comeron_tables\/Comeron_100kb_chr2R.txt"); 
    while(<O2R>){
	chomp;
	($pos,$rec)=split();
	$dico2R{$pos}=$rec;
    }
    close O2R;
    
    open(O3L, "Comeron_tables\/Comeron_100kb_chr3L.txt"); 
    while(<O3L>){
	chomp;
	($pos,$rec)=split();
	$dico3L{$pos}=$rec;
    }
    close O3L;
    
    open(O3R, "Comeron_tables\/Comeron_100kb_chr3R.txt"); 
    while(<O3R>){
	chomp;
	($pos,$rec)=split();
	$dico3R{$pos}=$rec;
    }
    close O3R;
  

    open(MAP, "$mapfile");
    open(RRC, ">${mapfile}.rrc");
    while(<MAP>){
	chomp;
	my @results;
	my $limit_string;
	my $limit;
	my ($start_rate, $end_rate, $center_rate, $center, $start_end);
        my ($start_comeron_rate, $end_comeron_rate, $center_comeron_rate);
	my @location = split(/\:/, $_);
	my $chromosome=$location[0];
	my @location_coord = split(/\../,$location[1]);
	my $start=$location_coord[0];
	my $end=$location_coord[1];
	
	if (!($chromosome =~ /^(2R|2L|3R|3L|4|X)/)) {
	    print "$chromosome Incorrect locus specification\n";
	} 
	else { 
	    if ($RELEASE eq "3") {
		$limit_string = $bp_limits_3{$chromosome};
	    } 
	    elsif ($RELEASE eq "4") {
		$limit_string = $bp_limits_4{$chromosome};
	    } 
	    elsif ($RELEASE eq "5.19") {
		$limit_string = $bp_limits_5{$chromosome};
	    }
	    elsif ($RELEASE eq "5.36") {
		$limit_string = $bp_limits_5{$chromosome};
	    }
	    $limit = $limit_string;
	    $limit =~ s/,//g;
	    if (!(($start >= 0) && ($start <= $limit) && ($end >= 0) && ($end <= $limit))) {
		print "Locus out of bounds on $chromosome: $limit_string; $limit \n";
	    } 
	    elsif ($start > $end) {
		print "Start coordinate greater than end coordinate \n"; 
	    } 
	    else {
		$center = ($start+$end)/2.0;
		print "Estimates on the chromosome $chromosome in the window [ $start - $end ]\t";
		$start_rate = get_rate($chromosome, $start, $RELEASE);
		$center_rate = get_rate($chromosome, $center, $RELEASE);
		$end_rate = get_rate($chromosome, $end, $RELEASE);
		
		$start_comeron_rate = comeron($chromosome, $start);
		$end_comeron_rate = comeron($chromosome, $end);
		$center_comeron_rate = comeron($chromosome, $center);

		printf ("\t\t-\t\t$_ \t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ($start_rate, $center_rate, $end_rate, $start_comeron_rate, $center_comeron_rate, $end_comeron_rate));
		printf RRC ("$_ \t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ($start_rate, $center_rate, $end_rate, $start_comeron_rate, $center_comeron_rate, $end_comeron_rate));
	    }
	}
    }
}

# Comeron recombination estimates 
sub comeron
{
    my ($chromosome,$pos)=@_;
    my $rate;
    
    if ($chromosome eq "X") {
	for my $key ( keys %dicoX ) {
	    my $nkey= $key + 100000;
	    if ($key<=$pos && $pos<=$nkey){
		$rate = $dicoX{$key};
		last;
	    }
	}
    }
    elsif ($chromosome eq "2L") {
	for my $key ( keys %dico2L ) {
	    my $nkey= $key + 100000;
	    if ($key<=$pos && $pos<=$nkey){
		$rate = $dico2L{$key};
		last;
	    }
	}
    }
    elsif ($chromosome eq "2R") {
	for my $key ( keys %dico2R ) {
	    my $nkey= $key + 100000;
	    if ($key<=$pos && $pos<=$nkey){
		$rate = $dico2R{$key};
		last;
	    }
	}
    }
    elsif ($chromosome eq "3L") {
	for my $key ( keys %dico3L ) {
	    my $nkey= $key + 100000;
	    if ($key<=$pos && $pos<=$nkey){
		$rate = $dico3L{$key};
		last;
	    }
	}
    }
    elsif ($chromosome eq "3R") {
	for my $key ( keys %dico3R ) {
	    my $nkey= $key + 100000;
	    if ($key<=$pos && $pos<=$nkey){
		$rate = $dico3R{$key};
		last;
	    }
	}
    }
    else{
	$rate = 0;
    }
    
    return $rate;
}


sub get_rate{
    my $chromosome=$_[0];
    my $pos=$_[1];
    my $RELEASE=$_[2];
    my $rate;
    my $coord = $pos / 1000000;
    
    if ($RELEASE eq "3") {
	if ($chromosome eq "2L") {
	    $rate = 2.4313 + 0.4474 * $coord - 0.0312 * $coord**2;
	} elsif ($chromosome eq "2R") {
	    $rate = -1.3784 + 0.7318 * $coord - 0.0258 * $coord**2;
	} elsif ($chromosome eq "3L") {
	    $rate = 2.8995 + 0.185 * $coord- 0.0171 * $coord**2;
	} elsif ($chromosome eq "3R") {
	    $rate = -1.5529 + 0.4632 * $coord - 0.0111 * $coord**2;
	} elsif ($chromosome eq "X") {
	    $rate = 1.2102 + 0.5966 * $coord - 0.0294 * $coord**2;
	} elsif ($chromosome eq "4") {
	    $rate = 0;
	}
    } 
    
    elsif ($RELEASE eq "4") {
	if ($chromosome eq "2L") {
	    $rate = 2.5183 + 0.4174 * $coord - 0.0297 * $coord**2;
	} elsif ($chromosome eq "2R") {
	    $rate = -1.5224 + 0.7254 * $coord - 0.0249 * $coord**2;
	} elsif ($chromosome eq "3L") {
	    $rate = 2.6663 + 0.2182 * $coord- 0.018 * $coord**2;
	} elsif ($chromosome eq "3R") {
	    $rate = -1.5567 + 0.466 * $coord - 0.0114 * $coord**2;
	} elsif ($chromosome eq "X") {
	    $rate = 1.1626 + 0.5992 * $coord - 0.0291 * $coord**2;
	} elsif ($chromosome eq "4") {
	    $rate = 0;
	}	    
    } 
	
    elsif ($RELEASE eq "5.19") { 
	if ($chromosome eq "2L") {
	    if ($coord>=0.3 and $coord<=17.84){
		$rate = 1.69634 + 0.65386 * $coord - 0.04197 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "2R") {
	    if ($coord>=2.05 and $coord<=20.66){
		$rate = -1.278363 + 0.66925 * $coord - 0.022218 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "3L") {
	    if ($coord>=0.86 and $coord<=19.45){
		$rate = 2.638370 + 0.244726 * $coord - 0.019554 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "3R") {
	    if ($coord>=4.53 and $coord<=25.74){
		$rate = -1.781510 + 0.502048 * $coord - 0.012588 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "X") {
	    if ($coord>=1.53 and $coord<=22.32){
		$rate = 1.132521 + 0.602636 * $coord - 0.029259 * $coord**2;
	    }
	    else{$rate = 0;}   
	} elsif ($chromosome eq "4") {
	    $rate = 0;
	}
    }
    
    elsif ($RELEASE eq "5.36") { 
	if ($chromosome eq "2L") {
	    if ($coord>=0.53 and $coord<=18.87){
		$rate = 2.58909 + 0.40558 * $coord - 0.02886 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "2R") {
	    if ($coord>=1.87 and $coord<=20.86){
		$rate = -1.435345 + 0.698356 * $coord - 0.023364 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "3L") {
	    if ($coord>=0.75 and $coord<=19.02){
		$rate = 2.940315 + 0.19103 * $coord - 0.017793 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "3R") {
	    if ($coord>=2.58 and $coord<=27.44){
		$rate = -1.629112 + 0.484498 * $coord - 0.012117 * $coord**2;
	    }
	    else{$rate = 0;}
	} elsif ($chromosome eq "X") {
	    if ($coord>=1.22 and $coord<=21.21){
		$rate = 1.151769 + 0.588228 * $coord - 0.029142 * $coord**2;
	    }
	    else{$rate = 0;}   
	} elsif ($chromosome eq "4") {
	    $rate = 0;
	}
    }
    
    if ($rate < 0) {
	$rate = 0;
	}
    
    return $rate;
}

sub help {
    print "RRC version 2.2\n\n";
    print "\tUsage:  RRC_2_2.pl -M map file -R genome sequence release (Release 5.36 by default)\n\n"; 
    print "\t\tThe map file should contain the coordinates of genomic windows such as:\n";
    print "\t\t2L:7500..9400\n\t\t2L:210918..250151\n\t\t...\n\n";
    print "\t\t NOTE: the recombination rate estimates are returned in a .rrc file!\n\n";
}
