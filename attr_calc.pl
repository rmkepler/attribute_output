#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::SimpleAlign;
use Bio::AlignIO;

#usage: ./attr_calc.pl TABLE FASTA outfile
#no alignments saved or folders created

my $table = shift; #The comma-separated table of  strain allelism from MLSTEZ
my $seqs = shift; #The fasta formated file of all alleles from MLSTEZ
my $outfile = shift; #Name of the table to be written

my $fasdb = Bio::DB::Fasta->new($seqs) or die "$seqs: $!";
my @ids = $fasdb->get_all_primary_ids;
my @attr = qw(psim glen uglen pindel);
my %hash;

#This section determines which strains have multiple alleles for each gene
#Creates a multidimensional hash; isolate key to a hash of genes for each isolate, the value of which is a hash of the attributes

open my $fh, "<", $table or die "$table: $!";
my $csv = Text::CSV->new ({
        binary    => 1, # Allow special character. Always set this 
        auto_diag => 1, # Report irregularities immediately
    });

my $header = $csv->getline($fh);
$csv->column_names (@$header); # use header

while (my $row = $csv->getline_hr ($fh)) {
    foreach my $col (@$header) {
        if ($col eq "Strain") {
            next
        }
         elsif ($row->{$col} eq "Yes"){
             my $strain = $row->{Strain} ;
             $hash{$strain}{$col} = align($col, $strain, @attr);
         }
     }
}
close $fh;

open (OUT, '>', $outfile); #begin writing to output file

print OUT "strain","\t","gene","\t",join("\t", @attr),"\n";

foreach my $indiv (sort keys %hash) { #cycles through sample names
    foreach my $gene (sort keys $hash{$indiv}) {
        my @attributes = ();
        foreach my $attr (@attr) {
                push (@attributes, $hash{$indiv}{$gene}{$attr})
            } 
            print OUT $indiv,"\t",$gene,"\t",join("\t", @attributes),"\n";
        }
    }
print OUT "\n\npsim is the percent similarity of non-gap alignment columns\n";
print OUT "glen is the length of the alignment with gaps\n";
print OUT "uglen is the length of the alingment with columns containins gaps removed\n";
print OUT "pindel is the percent of glen attributed to gaps, assumed to be indels\n";

close (OUT);

#subroutine here receives $col (the gene), $strain (the isolate) and an array of attributes to obtain
#retrives both alleles for a given strain and performs percent identity calculation
#returns a hash of attribute names (key) and their values 

sub align {
    my ($gene, $isolate, @attr) = @_;
    my @seq_array;
    my @isolate_array = grep {/^${isolate}_${gene}_*/} @ids;
    foreach my $id (@isolate_array) {
        my $seq = $fasdb->get_Seq_by_id($id);
            push (@seq_array, $seq);
        }
    #make the alignment with MAFFT. Params default.
    my $factory = Bio::Tools::Run::Alignment::MAFFT->new;
    my $aln = $factory->align(\@seq_array);
    my $ugaln = $aln->remove_gaps;
    #make the hash of alignment attributes
    my %attr_hash;
    foreach my $attr (@attr) {
        if ($attr eq "psim") {
	    $attr_hash{$attr} = $aln->percentage_identity();
		}
		elsif ($attr eq "glen") {
		    $attr_hash{$attr} = $aln->num_residues;
		}
		elsif ($attr eq "uglen") {
		    $attr_hash{$attr} = $ugaln->num_residues;
		}
		elsif ($attr eq "pindel") {
		    $attr_hash{$attr} = ((($attr_hash{"glen"}-$attr_hash{"uglen"})/$attr_hash{"glen"})*100);
		}
    }
    return \%attr_hash;
}
__END__
