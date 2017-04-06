#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# -----------------------------------------------------------------------------
# GLOBALS

use vars qw ($help $file $zip);

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions (
"i=s"       => \$file,
"z"         => \$zip,
"help"      => \$help,
"h"         => \$help);
usage() if ($help || !$file);

# -----------------------------------------------------------------------------
# MAIN

printHeader();

if(`gzip -t $file 2>&1 1>/dev/null` !~ /not in gzip format/){
    open(FILE, "gunzip -c $file | ") || die "cannot open $file\n";
}else{
    open(FILE, "<$file") || die "cannot open $file\n";
}
if($zip){
    open(OUT1, " | gzip -c > $file\_1.gz") || die "cannot open $file\_1.gz\n";
    open(OUT2, " | gzip -c > $file\_2.gz") || die "cannot open $file\_2.gz\n";
}else{
    open(OUT1, ">$file\_1") || die "cannot open $file\_1\n";
    open(OUT2, ">$file\_2") || die "cannot open $file\_2\n";
}
while(<FILE>){
    chomp;
    my @line = split(/\s+/,$_); my $length = 0;
    foreach(@line){if(/length/){$_=~s/length\=//;$length=$_/2;}}
    $_ =~ s/length\=\d+/length\=$length/;
    print OUT1 "$_\n";
    print OUT2 "$_\n";
    my $newline = <FILE>; chomp($newline);
    print OUT1 substr($newline, 0, length($newline)/2)."\n";
    print OUT2 substr($newline, length($newline)/2, length($newline)/2)."\n";
    $newline = <FILE>; chomp($newline);
    $newline =~ s/length\=\d+/length\=$length/;
    print OUT1 "$newline\n";
    print OUT2 "$newline\n";
    $newline = <FILE>; chomp($newline);
    print OUT1 substr($newline, 0, length($newline)/2)."\n";
    print OUT2 substr($newline, length($newline)/2, length($newline)/2)."\n";
}
close(OUT1);close(OUT2);close(FILE);

# -----------------------------------------------------------------------------
# FUNCTIONS

sub usage {
    print STDERR "usage: splitPairedEndReads.pl -i <file>\n";
    print STDERR "\n";
    print STDERR "[INPUT]\n";
    print STDERR " -i <file>    fastq file (.fastq or .fastq.gz)\n";
    print STDERR " -z           zip output files\n";
    print STDERR " -h           this (usefull) help message\n";
    print STDERR "[VERSION]\n";
    print STDERR " 03-06-2012\n";
    print STDERR "[BUGS]\n";
    print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
    print STDERR "\n";
    exit(-1);
}

sub prettyTime{
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, 
    $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    return "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}

sub printHeader{
    print "\n";
    print "splitPairedEndReads.pl startet at ".prettyTime()."\n";
    print "\n";
    print "input:  $file\n";
}
