#!/usr/local/bin/perl -w

use strict;

my $regionfile=$ARGV[0];
my $goterm=$ARGV[1];

my $GO_obofile="gene_ontology_ext2015.obo";
my $GOSLIMfile="ATH_GO_GOSLIM2015.txt";
my $gfffile="TAIR10_GFF3_genes.gff";

my %parents;
my %children;
my %descriptions;

#First parse and read in GO
open(IN,"<$GO_obofile");
while(<IN>){
   my $line=$_;
   chomp($line);
   if($line=~/^\[Term\]/){
      my $curid;
      my $curname;      
       
      #Fetch id
      do{
        $_=<IN>;
        $line=$_;
        chomp($line);
      }while($line!~/^id: (GO:\w+)/);
      $curid=$1;
      
      #Fetch name
      do{
        $_=<IN>;
        $line=$_;
        chomp($line);
      }while($line!~/^name: (.+)/);
      $curname=$1;
      $descriptions{$curid}=$curname;
      
      #Fetch parents
      while(length($line)>0){
        if($line=~/^is_a: (GO:\w+).*/){
          my $curparent=$1;
          push(@{$parents{$curid}},$curparent);
          push(@{$children{$curparent}},$curid);   
        }
        if($line=~/^relationship: part_of (GO:\w+).*/){
          my $curparent=$1;
          push(@{$parents{$curid}},$curparent);
          push(@{$children{$curparent}},$curid);   
        }
        if($line=~/^relationship: regulates (GO:\w+).*/){
          my $curparent=$1;
          push(@{$parents{$curid}},$curparent);
          push(@{$children{$curparent}},$curid);   
        }
        if($line=~/^relationship: positively_regulates (GO:\w+).*/){
          my $curparent=$1;
          push(@{$parents{$curid}},$curparent);
          push(@{$children{$curparent}},$curid);   
        }
        if($line=~/^relationship: negatively_regulates (GO:\w+).*/){
          my $curparent=$1;
          push(@{$parents{$curid}},$curparent);
          push(@{$children{$curparent}},$curid);   
        }
        $_=<IN>;
        $line=$_;
        chomp($line);
      }  
   }
}
close IN;

#Next find all the GO ids that contain the term we're looking for
my %endlist=();
foreach my $key(keys %descriptions){
  my $curname=$descriptions{$key};
  if($curname=~/.*$goterm.*/){
    $endlist{$key}=1;
  }
}

#Now find all the genes that have a label that is in the list
open(IN,"<$GOSLIMfile");
my %genelist;
while(<IN>){
  my $line=$_;
  chomp($line);
  my @vals=split/\t/,$line;
  my $gene=$vals[0];
  my $golabel=$vals[5];
  #my $godesc=$vals[3]." ".$vals[4];
  if(defined($endlist{$golabel})){
    push(@{$genelist{$gene}},$golabel);
  }  
}
close IN;

my @golist=keys %genelist;

#Now read and analyze the regions
#First parse GO
my %GO;
my %annotations;
my %gstart;
my %gend;
my %gchrom;
my %genes;

open(IN,"<$GOSLIMfile");
while(<IN>){
  my $line=$_;
  chomp($line);
  my @vals=split/\t/,$line;
  my $gene=$vals[0];
  my $golabel=$vals[5];
  my $godesc=$vals[3]." ".$vals[4];
  $GO{$golabel}=$godesc;
  push(@{$annotations{$gene}},$golabel);
}
close IN;

#Next parse GFF to get gene locations
my %chrom;
open(IN,"<$gfffile");
my $numgenes=0;
while(<IN>){
  my $line=$_;
  chomp($line);
  my @vals=split/\t/,$line;
  if($vals[2] eq "gene"){
     #Extract gene name
     $line=~/.*ID=([^;]+);.*/;
     my $atcode=$1;
     $genes{$atcode}=1;
     $gstart{$atcode}=$vals[3];
     $gend{$atcode}=$vals[4];
     my $chr=$vals[0];
     $chr=~/^Chr(.+)/;
     $gchrom{$atcode}=$1;
     push(@{$chrom{$gchrom{$atcode}}},$atcode);
     $numgenes++;
  }
}
close IN;

#Now parse the region file
open(IN,"<$regionfile");
while(<IN>){
  my $line=$_;
  chomp($line);
  my @vals=split/\t/,$line;
  if(scalar(@vals)==9){
    print "$line\n";  
    my @genespresent=();  
    my $curchrom=$vals[1];
    my $curstart=$vals[6];
    my $curend=$vals[7];
    
    #Read SNP positions and logp
    my @snppositions=();
    my @snplogps=();
    my $pos=tell;
    while(<IN>){
      $line=$_;
      chomp($line);
      @vals=split/\t/,$line;
      if(scalar(@vals)==3){
        push(@snppositions,$vals[1]);
        push(@snplogps,$vals[2]);        
        #print "".($vals[1])."\t".($vals[2])."\n";
      }
      else{
        #Rewind to previous position
        seek(IN,$pos,0);
        last;        
      }
      $pos=tell;
    }
    
    #Order genes on the current chromosome
    my @sortedgenes=sort{$gstart{$a}<=>$gstart{$b}} @{$chrom{$curchrom}};
    
    #Check the genes that are contained within this region
    foreach my $gene(@sortedgenes){
      if( ($gstart{$gene} > $curstart) and
           ($gend{$gene} < $curend)){
        push(@genespresent,$gene);
        #print "$gene ".($gchrom{$gene})."(".($gstart{$gene})."-".($gend{$gene})."):\n";
        #foreach my $label(uniq(@{$annotations{$gene}})){
        #   print "$label\t".($GO{$label})."\n";
        #}
      }
    }
    #Check how many of the present genes can be associated to the GO label
    my %genespresenth=map{$_ =>1} @genespresent;
    my %golisth=map{$_=>1} @golist;
    my @geneswithgo = grep( $genespresenth{$_}, @golist );
    my $tgenes=scalar(@genespresent);
    my $fgenes=scalar(@geneswithgo);
    print "$tgenes genes found in region, $fgenes of these associated to GO label\n";
    if($tgenes >0){
      foreach my $gene(@genespresent){
        print "\t$gene ".($gchrom{$gene})."(".($gstart{$gene})."-".($gend{$gene})."):\n";
        #See if we find SNPS in the gene or 1000bp upstream of the start (promoter)
        for(my $s=0;$s<scalar(@snppositions);$s++){
          if( ($gstart{$gene}<=$snppositions[$s]) and ($snppositions[$s]<=$gend{$gene})){
             print "\t\tSNP\t".($snppositions[$s])."\t".($snplogps[$s])."\tGENE\n";   
          }
          elsif( ($gstart{$gene}-1000<=$snppositions[$s]) and ($snppositions[$s]<$gstart{$gene})){
             print "\t\tSNP\t".($snppositions[$s])."\t".($snplogps[$s])."\tPROMOTER\n";    
          }
        }
        #Print GO terms
          foreach my $label(uniq(@{$annotations{$gene}})){
             print "\t\t$label\t".($GO{$label})."\n";
          }
       }
    }         
    
    #if($fgenes >0){
    #   my $prob=combi($fgenes,scalar(@golist))*combi($tgenes-$fgenes,scalar(keys %genes)-scalar(@golist))/combi($tgenes,scalar(keys %genes));   
    #   print "Prob=$prob\n";
    #}
  }
}
close IN;

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub fac {
  $_[0]>1?$_[0]*fac($_[0]-1):1;
}

sub combi{
  (my $a, my $b)=@_;
  return(fac($b)/(fac($a)*fac($b-$a)));  
}


