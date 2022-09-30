#!/usr/bin/env perl
=head1 Generate phylogenetic tree based on multiple sequence alignment in fasta format.

    p3x-retrieve-family-sequences-for-genomes.pl [options] families genomes

This script runs a multiple sequence alignment program for all pgfam sequences for a set of genomes.
Protein families and genomes can be specified as either a comma-separated list or by a file with one p[gl]fam (or SOG) or genome_id per line.
Output is one fasta file per protein family.
Can optionally align with mafft, yielding an aligned fasta (.afa) file per family.

=head2 Parameters

The two positional parameters are the protein families (comma-separated or file of ids) and genomes (comma-separated or file of ids).

The command-line options are as follows.

=over 4

=item program

mafft or muscle (case insensitive) (default mafft)

=item alphabet

DNA or protein (default DNA)

=item output_dir

Directory where alignment files and log file will be written. Defaults to current directory.

=item verbose

Output extra comments and leave intermediate files

=back

=cut

use strict;
use Cwd qw(getcwd);
use P3DataAPI;
use P3Utils;
#use File::Copy ('copy', 'move');
#use File::Spec;

#$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts(['families', 'genomes'],
                ['align', 'Align sequences with mafft'],
                ['threads|t', 'Threads for mafft'],
                ['alphabet|a=s', 'DNA or protein (default DNA)', { default => 'DNA'}],
                ['output_dir|o=s', 'Output directory.' ],
                ['verbose|debug|v', 'Write status messages to STDERR'],
        );
# Check the parameters.
my ($families, $genomes) = @ARGV;
if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}
if (! $families and $genomes) {
    print "Fewer than needed two inputs specified.";
    exit(1);
}
if ($opt->output_dir) {
    die "output dir " .$opt->output_dir . " does not exist. You must create it first." unless -d $opt->output_dir;
}

my @families; # the list of protein families
if ( -f $families ) {
    print STDERR "looking for family ids in file $families\n" if $opt->verbose;
    open my $f, $families or die "cannot open $families";
    while (<$f>) {
        if (/PGF\S+/ or /PLF\S+/ or /SOG\S+/) {
            push @families, $&;
        }
    }
}
else { # try splitting on commas
    print "Try splitting families out of $families\n";
    my @candidates = split(",", $families);
    for my $candidate (@candidates) {
        push @families, $candidate if $candidate =~ /PGF\S+|PLF\S+|SOG\S+/;
    }
}

my $num_families = scalar @families;
print "Got $num_families families.\n", join(",", @families), "\n\n" if $opt->verbose;

my @genome_list; # the list of genomes
if ( -f $genomes ) {
    print STDERR "looking for genome_ids in file $genomes\n" if $opt->verbose;
    open my $f, $genomes or die "cannot open $genomes";
    while (<$f>) {
        if (/\d+\.\d+/) {
            push @genome_list, $&;
        }
    }
}
else { # try splitting on commas
    print "Try splitting genomes out of $genomes\n";
    my @candidates = split(",", $genomes);
    for my $candidate (@candidates) {
        push @genome_list, $candidate if $candidate =~ /\d+\.\d+/;
    }
}
my $num_genomes = scalar @genome_list;
print "Got $num_genomes genomes.\n". join(",", @genome_list). "\n\n" if $opt->verbose;

my $api = P3DataAPI->new;
my $md5_type = 'na_sequence_md5';
$md5_type = 'aa_sequence_md5' if $opt->alphabet =~ /protein/;
my $threads = 1;
$threads = $opt->threads if $opt->threads;
print "alphabet=", $opt->alphabet, ", md5 type = $md5_type, threads = $threads\n";

my $original_dir = getcwd();
if ($opt->output_dir) {
    chdir($opt->output_dir);
}

for my $fam (@families) {
    print "Try getting features for $fam from genomes.\n";
    my $feat_aref = get_features_for_one_family($api, $fam, \@genome_list, $md5_type);
    print "result = $feat_aref\n";
    print "num items in result = ", scalar @$feat_aref, "\n";
    next if @$feat_aref == 0;
    my %md5;
    for my $item (@$feat_aref) {
        $md5{$item->{$md5_type}} = 1;
        print("$item->{patric_id} : $item->{$md5_type}\n");
    }
    print "Now get sequences for ", scalar(keys %md5), " uniq md5s.\n";
    my @md5s = keys %md5;
    my $seq_hash = $api->lookup_sequence_data_hash(\@md5s);
    my $unaligned_output_file = "${fam}_unaligned.fa"; 
    print("Writing sequenes to $unaligned_output_file\n");
    open my $f, ">", $unaligned_output_file;
    for my $item (@$feat_aref) {
        print $f ">$item->{patric_id}\n$seq_hash->{$item->{$md5_type}}\n";
    }
    close $f;
    if ($opt->align) {
        my $aligned_output_file = "${fam}_aligned.afa";
        print "aligning to $aligned_output_file\n";
        my $command = ("mafft --auto --threads $threads $unaligned_output_file > $aligned_output_file"); #using a string instead of array because redirection is needed - better option is to use open and captuer output and write to file -- TODO
        #print "alignment command = ". join(" ", @command), "\n";
        print "alignment command = $command\n";
        system($command);
    }
}
chdir($original_dir);

sub get_features_for_one_family {
    my ($p3api, $family, $genomes_aref, $md5_type) = @_;
    my $family_type = '';
    print "Test type of $family\n";
    if ($family =~ /PGF_\d+/) {
        $family_type = 'pgfam_id';
    }
    elsif ($family =~ /PLF_\d/) {
        $family_type = 'plfam_id';
    }
    elsif ($family =~ /SOG_\d+/) {
        $family_type = 'sog_id';
    }
    else { die "could not figure out family type for $family\n"}
    print "search for $family_type = $family in genomes: ", join(', ', @$genomes_aref), "\n";
    my @features;
    $p3api->query_cb(
        "genome_feature",
        sub {
            my ($data) = @_;
            push @features, @$data;
            return 1;
        },
        [ "eq",     $family_type, $family ],
        [ "in",     "genome_id",   "(".join(",", @$genomes_aref).")" ],
        [ "select", "patric_id,$md5_type" ]    );
    return \@features;
}

