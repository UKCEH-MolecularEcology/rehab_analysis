#!/usr/bin/perl -w

## Help text in perldoc format
=head1 SYNOPSIS

get_all_contigs.pl - returns a tab list with name, filename, and size of sequences 
                      in all .fa or .fasta files in a folder, concatenated to the standard output,
                      and some stats on the standard error.

=head1 USAGE

get_all_contigs.pl [-f format] folder_path

=head1 AUTHOR

Lionel Guy (lionel.guy@ebc.uu.se)

=cut

## load perl packages
use strict;
use Getopt::Std;
use Bio::SeqIO;
use File::Basename;  # Added for file basename

## handle input parameters
our $opt_f = 'fasta';
getopts('f:');

usage() unless @ARGV;

## subroutines
## Display helptext and exit program
sub usage {
    system("perldoc $0");
    exit;
}

## instatiate variables
my ($count, $sum) = (0, 0);
my ($min, $max);
my $folder_path = shift @ARGV;

## print header
print "bin\t0_Contig\tlength\n";

## Find all .fa and .fasta files in the folder
opendir(my $dh, $folder_path) or die "Cannot open directory $folder_path: $!";
my @files = grep {/\.fa$|\.fasta$/i && -f "$folder_path/$_"} readdir($dh);
closedir($dh);

## Iterate over found files
foreach my $file (@files) {
    my $filename = fileparse($file, qr/\.[^.]*$/);  # Extract basename without extension
    my $full_path = "$folder_path/$file";
    
    ## use Bio::SeqIO package to read input file
    my $seq_io = Bio::SeqIO->new(-file => $full_path, -format => $opt_f);

    ## Iterate over sequences in input file to 
    ## generate individual sequence and summary 
    ## statistics.
    while (my $seq = $seq_io->next_seq) { ## selects one sequence at a time
        ## set variables for THIS sequence
        my $id = $seq->display_id;
        my $len = $seq->length;
        
        ## print output for THIS sequence 
        print "$filename\t$id\t$len\n";
        
        ## set stats for ALL sequences
        $count++; ## seq count
        $sum += $len; ## add length of this sequence to all other sequences
        $min = $len if (!$min || $len < $min); ## find shortest sequence
        $max = $len if (!$max || $len > $max); ## find longest sequence
    }
}

my $average = int($sum / $count + 0.5); ## Compute average

