#!/usr/bin/perl

# transform the numbers to the strings according to the table file
# basically, this script is for untransform the numbers transformed by 
# the transnum.pl

$ARGC = @ARGV;
if ( $ARGC < 1 ){
      # error routine
  printf ("untransnum.pl: output-table-file [separator] < input-file > output-file\n");
  exit (1);
}
    # initialization
open ( TABLEFILE, "<$ARGV[0]" );
@table = <TABLEFILE>;
%numbers = ();
$c=0;
$sep=" ";
if ( $ARGC >1 ){ $sep = $ARGV[1]; }  # separator

    # read transform-table
foreach $trans( @table ) {
    @eles = split(" ", $trans);
    if ( $#eles == 0 ){ $numbers{$c} = $eles[0]; }
    else { $numbers{$eles[0]} = $eles[1]; }
    $c++;
}

    # read transform-file
while (<STDIN>){
  chomp;
  $_ =~ s/$sep$sep/$sep/g;
  @eles = split($sep, $_);
  $all = @eles;
  $c = 0;
  foreach $item( @eles ){
	if ( $item < 0 ){
      print "*";
    } elsif (!exists $numbers{$item}) { 
      print "$item";
    } else { print "$numbers{$item}"; }
    $c++;
    if ($c<$all){ print " ";}
  }
  print "\n"
}
