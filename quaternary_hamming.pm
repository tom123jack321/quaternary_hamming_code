#!/usr/bin/env perl
##!/Volumes/820g/Mac_app/homebrew/bin/perl
##!/usr/bin/perl
# naming convention: 1. A: array; 2. H: hash; 3. R: reference; 4. Sr: subroutine; 5.
#mark a#
package DNAstorage;

use v5.24.1;
use utf8;
use strict;
use warnings;
# use YAML;   # say Dump
use POSIX qw/ ceil floor /;
use Data::Dump qw/ dump /;  #dump @AoA;
use Exporter qw/ import /;
use Storable qw/ store retrieve /;
use List::Util qw/ sum /;
#mark b#

our @ISA = qw/ Exporter /;
our @EXPORT_OK = qw / 
  /;
#mark c#

my $eps = 1.e-9;
# my $a_in = "0123";      # $c : 0102123
my $a_in = "0123" . "2311" . "010";   # $c : 320012302311010
print '$a_in   '; say $a_in;
my $c = & encode_q_hamming ($a_in, );
print '$c   '; say $c;  

$c = "320012302311011";   # the last digit is changed to 1
$c = "320012302311030";   # the -2nd digit is changed to 3
$c = "120012302311010";   # the 1st digit is changed to 1
# $c = "0102120";   # the last digit changed to 0
my $a_corrected = & decode_q_hamming ($c, );
print '$a_corrected   '; say $a_corrected;

#mark z#

sub encode_q_hamming (){
  # $a_in : the input number sequence
  my ($a_in, ) = @_;
  my @a = split //, $a_in;
  my $n = scalar @a;
  # $m : the length of parity checksum series
  my $m = 0;
  ++$m until 2**$m > $n + $m;
  # print '$m   '; say $m;
  my (@c, );
  # initialise @c
  for my $i (1 .. $m){
    for my $j ( (2**($i-1)+1) .. (2**$i-1) ){
      $c[$j-1] = $a[$j-$i-1] if $j <= $n+$m; 
    }
  }
  # print '@a  '; map {print "|$_"; } @a; say '|';
  print 'dump @c   '; say dump @c;
  for my $i (1 .. $m){
    my $p = 0;          # $p : checksum
    my $k_max = 0;
    ++$k_max until 2**($i-1) *(2*$k_max+1) > $n + $m;
    --$k_max;
    for my $k (0 .. $k_max ){
      my $j_sta = 2**($i-1) * (2*$k+1);
      my $j_end = 2**($i-1) * (2*$k+2) -1;
      for my $j ( $j_sta .. $j_end ){
        die '($i, $k, $j_sta, $j_end, $j, ) ' . "($i, $k, $j_sta, $j_end, $j, )"
          if $j > 2**($i-1) and not defined $c[$j-1];
        $p += $c[$j-1] if $j > 2**($i-1) and $j <= $n+$m;
      }
    }
    $p %= 4;
    $c[2**($i-1)-1] = $p;
  }
  return join '', @c;
}

sub decode_q_hamming (){
  # $c_in : codeword. 
  # The codeword might differ from the correct one in one number at most
  my ($c_in, ) = @_;
  my @c = split //, $c_in;
  my $len_c = scalar @c;
  my $m = ceil( log($len_c)/log(2) );
  my $n = $len_c - $m;
  my $v = 0;
  # attention: to use indices later, minus 1
  my @set_index_p = 1 .. $m;
  my ($set_index_each_p_R, @w, );
  for my $i (@set_index_p){
    my $p = 0;          # $p : checksum
    my $k_max = 0;
    ++$k_max until 2**($i-1) *(2*$k_max+1) > $n + $m;
    --$k_max;
    for my $k (0 .. $k_max ){
      my $j_sta = 2**($i-1) * (2*$k+1);
      my $j_end = 2**($i-1) * (2*$k+2) -1;
      for my $j ( $j_sta .. $j_end ){
        die '($i, $k, $j_sta, $j_end, $j, ) ' . "($i, $k, $j_sta, $j_end, $j, )"
          if $j > 2**($i-1) and not defined $c[$j-1];
        if ( $j > 2**($i-1) and $j <= $n+$m ){
          $p += $c[$j-1];
          push @{ $set_index_each_p_R->[$i-1] }, $j;
        }
      }
    }
    $p %= 4;
    $w[$i-1] = $p;
    if ( abs($c[2**($i-1)-1] - $p)<$eps ){ 
      $v += 2**($i-1) * 0;
    } else { 
      $v += 2**($i-1) * 1;
    }
  }
  # print 'dump $set_index_each_p_R     '; say dump $set_index_each_p_R  ;
  # print '$v   '; say $v;
  # if no errors or only one of the checksum digits is wrong
  if ( abs($v)<$eps or grep { abs($v-$_)<$eps } @set_index_p ){
    # $a_crt : corrected a
    my $a_crt = & extract_data_qh (\@c, $m);
    return $a_crt;
  }
  # if one element of the data sequence is wrong
  for my $i (@set_index_p){
    if ( grep { abs($v-$_)<$eps } @{$set_index_each_p_R->[$i-1]} ){
      my $p_i = $c[2**($i-1)-1];
      my $c_v = $c[$v-1];         # the wrong element
      $c_v = ($p_i - ($w[$i-1]-$c_v) )%4;     # the correct one
      $c[$v-1] = $c_v;
      my $a_crt = & extract_data_qh (\@c, $m);
      return $a_crt;
    }
  }
}

#mark d#

# my $c_R = [ split //, $c ];
# my $len_c = length $c;
# my $m = ceil( log($len_c)/log(2) );
# my $a_crt = & extract_data_qh ($c_R, $m);
# print '@$a_crt  '; map {print "|$_"; } @$a_crt; say '|';
sub extract_data_qh (){
  # $c_R : reference to the array of codeword
  # $m : the number of checksum digits
  my ($c_R, $m, ) = @_; 
  my @c = @{$c_R};
  for my $i (1 .. $m){
    undef $c[2**($i-1)-1];
  }
  my @a;
  for (@c){
    push @a, $_ if defined $_ 
  }
  return join '', @a;
}
