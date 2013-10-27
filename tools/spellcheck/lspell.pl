#!/bin/perl
# From http://www.kegel.com/kerspell/

# C comment spell checker
# For each given source file, print the filename, a colon, and the number
# of misspelled words, then a list of misspelled words.
# Words contained in the file stopwords.txt are not considered spelling errors.
# Copyright 2003, Dan Kegel.  Licensed under GPL.  See the file ../COPYING for details.

sub check_content($) {
	my $content = shift;
	$content =~ tr/*/ /;
	print POUT "$content\n";
}

$TEMPFILE="/tmp/spell.tmp";
$STOPFILE=shift(@ARGV);

open(STOPFILE, $STOPFILE) || die "can't open stopword file $STOPFILE";
while (<STOPFILE>) {
	chomp;
	$stopped{$_}++;
}
close(STOPFILE);

foreach $file (@ARGV) {
	open (FI, $file) || die $file;
	$content = join ("", <FI>);
	close (FI);

	open(POUT, "> $TEMPFILE") || die;
	$content =~ s!//(.+)$!check_content($1)!egm;
	$content =~ s!/\*(.+?)\*/!check_content($1)!egs;
	close(POUT);

  `sed -i 's/[{}]/XXX/' $TEMPFILE`;

  # open(PIN, "aspell -a -e -l en < $TEMPFILE | sort -uf |") || die;
  open(PIN, "aspell -a -e -l en < $TEMPFILE 2> /dev/null | sort | cut -f 2 -d ' ' | uniq |") || die;
	undef @badwords;
	while (<PIN>) {
		chomp;
    if ($stopped{$_} == 0) 
    {
			push(@badwords, $_);
    }
	}
	close(PIN) || die;

  if (@badwords) 
  {
    foreach (@badwords)
    {
      print "$file: $_\n";
    }
	}
  # else
  # {
  #   print "no bad for $file\n";
  # }
  # print "===================================\n\n";
}
