#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::SimpleAlign;
use Bio::AlignIO;

#usage = ./pair_calc.pl TABLE FASTA outfile

my $table = shift; #the comma-separated table of  strain allelism from MLSTEZ
my $seqs = shift; #the fasta formated file of all alleles from MLSTEZ
my $outfile = shift; #name of the table to be written

my $fasdb = Bio::DB::Fasta->new($seqs);
my @ids = $fasdb->get_all_primary_ids;

#this section determines which strains have multiple alleles for each gene
#creates a hash, with strains as key, values are another hash with gene as key and percent similarity as value

open my $fh, "<", $table or die "$table: $!";
my $csv = Text::CSV->new ({
        binary    => 1, # allow special character. Always set this 
        auto_diag => 1, # report irregularities immediately
    });

my $header = $csv->getline($fh); 
my %hash;
my %outhead;
$csv->column_names (@$header); # use header

while (my $row = $csv->getline_hr ($fh)) {
    foreach my $col (@$header) {
         if ($row->{$col} eq "Yes"){
             my $strain = $row->{Strain} ;
             my $psim = pssub($col, $strain); #send to subroutine to align and calculate attributes
             $outhead{$col} = 1;
             $hash{$strain}{$col} = $psim; #add values to %hash
         }
         else {next;}
     }
}
close $fh;

open (OUT, '>', $outfile); #begin writing to output file

my @keys = sort keys %hash; #sorted names of allelic isolates

print OUT join("\t", @$header),"\n"; #writes tab delimited file header

foreach my $sample_key (sort keys %hash) { #cycles through each sample name
        my @out_line = (); #empty array that is populated with percentages in loop below
        foreach my $gene (sort keys %outhead) { #this loop queries genes contained in %outhead
            if (exists $hash{$sample_key}{$gene}) { 
            push (@out_line, $hash{$sample_key}{$gene}); #get percentage 
            } else { 
                    push (@out_line, ''); #if the sample does not have percentage, push an empty value
                    }
        }
        print OUT $sample_key,"\t",join("\t", @out_line),"\n";
    }

#subroutine here receives $col (the gene) and $strain (the isolate)
#retrives both alles for a given strain and performs percent identity calculation
#returns percent identity
#creates folders for each gene and saves an alignment of allels for each isolate

sub pssub {
    my $gene = shift;
    my $isolate = shift;
    my @seq_array;
    my $seq_array_ref =\@seq_array;
    my @isolate_array = grep {/^${isolate}_${gene}_*/} @ids;
    foreach my $id (@isolate_array) {
        my $seq = $fasdb->get_Seq_by_id($id);
            push (@seq_array, $seq);
        }
    #make the alignment with MAFFT. Params default.
    my $factory = Bio::Tools::Run::Alignment::MAFFT->new;
    my $aln = $factory->align($seq_array_ref);
    my $psim = $aln->percentage_identity();
    
    #makes folders and writes allele alignments
    if (!-d $gene) {
        mkdir $gene;
    }
    my $out = Bio::AlignIO->new(-file => ">$gene/$isolate.nex", -format => 'nexus'); {
        $out->write_aln($aln);
    }
    return ($psim);
}
__END__
