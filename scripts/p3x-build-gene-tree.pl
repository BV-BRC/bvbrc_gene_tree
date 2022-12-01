#!/usr/bin/env perl
=head1 Generate phylogenetic tree based on multiple sequence alignment in fasta format.

    p3x-build-gene-tree.pl [options] msaFile

This script runs a phylogenetic tree-building program on the input alignment.
Programs include RAxML, PhyML, and FastTree
Output is a newick formatted tree as well as a log file.

=head2 Parameters

The one positional parameter is the path of the fasta-formatted MSA file.

The command-line options are as follows.

=over 4

=item program

raxml, phyml, or fasttree (case insensitive)

=item model

defaults to GTR for DNA and LG for protein, must be one recognized by program

=item alphabet

DNA or protein (inferred from sequences if missing)

=item output_dir

Directory where newick tree file and log file will be written. Defaults to current directory.

=item verbose

Output extra comments and leave intermediate files

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use File::Copy ('copy', 'move');
use File::Spec;
use Tree_Builder;

#$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('msaFile',
                ['program|p=s', 'Program to use to build tree (raxml, phyml, or fasttree) (default raxml)', { default => 'raxml' }],
                ['model|m=s', 'Substitution model (or use GTR for DNA, LG for protein)'],
                ['alphabet|a=s', 'DNA or protein (infer if not specified)'],
                ['bootstrap|b=n', 'Perform n bootstrap pseudoreplicates for branch support values', { default => 0 }],
                ['output_dir|o=s', 'Output directory. (default "./")', { default => "./" }],
                ['verbose|debug|v', 'Write status messages to STDERR'],
        );
# Check the parameters.
my $msaFile = $ARGV[0];
if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}
if (! $msaFile) {
    print "No input file specified.";
    exit(1);
}
if (! -f $msaFile) {
    print "Cannot find alignment file $msaFile";
    exit(1);
}
my $output_dir = File::Spec->rel2abs($opt->output_dir);
# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    Tree_Builder::set_debug($debug)
}

my $alphabet;
if ($opt->alphabet) {
    $alphabet = $opt->alphabet; 
}
my $tree_builder = new Tree_Builder($msaFile, $alphabet);

if ($opt->model) {
   $tree_builder->set_model($opt->model);
}
$tree_builder->set_output_dir($output_dir);
my $treeFile;
if ($opt->program eq 'raxml') {
    $treeFile = $tree_builder->build_raxml_tree($opt->bootstrap);
} elsif ($opt->program eq 'phyml') {
    $treeFile = $tree_builder->build_phyml_tree($opt->bootstrap);
} elsif ($opt->program eq 'fasttree') {
    $treeFile = $tree_builder->build_fasttree($opt->bootstrap);
} else {
    die "Unrecognized program: $opt->program \n";
}
my $logFile = $tree_builder->get_log_file();
print "Output directory is $output_dir\n";
print "Tree file is $treeFile\nLog file is $logFile\n";

