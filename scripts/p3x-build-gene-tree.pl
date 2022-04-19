#!/usr/bin/env perl
=head1 Generate phylogenetic tree based on multiple sequence alignment in fasta format.

    p3x-infer-gene-tree.pl [options] msaFile

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

=item verbose

Output extra comments and leave intermediate files

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use File::Copy ('copy', 'move');
use Sequence_Alignment;

#$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('msaFile',
                ['program|p=s', 'Program to use to build tree (raxml, phyml, or fasttree) (default raxml)', { default => 'raxml' }],
                ['model|m=s', 'Substitution model (or use GTR for DNA, LG for protein)'],
                ['alphabet|a=s', 'DNA or protein (infer if not specified)'],
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
# Get the debug flag.
my $debug = $opt->verbose;

my $alphabet;
if ($opt->alphabet) {
    $alphabet = $opt->alphabet; 
}
else {
    open F, $msaFile or die "cannot open $msaFile";
    my $dnaCount = 0;
    my $proteinCount = 0;
    while (<F>) {
        next if /^>/;
        $dnaCount += tr/ACGTacgt/ACGTacgt/;
        $proteinCount += tr/DEFHIKLMNOPQRSVWYdefhiklmnopqrsvwy/DEFHIKLMNOPQRSVWYdefhiklmnopqrsvwy/;
    }
    print "Autodetecting alphabet: DNA = $dnaCount, protein = $proteinCount\n" if $debug;
    $alphabet = ('DNA', 'protein')[$proteinCount > $dnaCount];
}
my $model;
if ($opt->model) {
    $model = $opt->model; # default for protein
}
elsif ($alphabet =~ /DNA/i) {
        $model = "GTR"
}
elsif ($alphabet =~ /protein/i) {
        $model = "LG"
}
my $program = "raxml"; #default
if ($opt->program) {
    $program = $opt->program;
    $program = lc($program);
}
my $output_name = $msaFile;
$output_name =~ s/\..{2,4}$//;

print STDERR "About to call tree program $program\n" if $debug;
my @output_files;
if ($program eq 'raxml') {
    @output_files = run_raxml($msaFile, $alphabet, $model, $output_name);
} elsif ($program eq 'phyml') {
    @output_files = run_phyml($msaFile, $alphabet, $model, $output_name);
} elsif (lc($program) eq 'fasttree') {
    #$alignment->write_fasta($msaFile);
    @output_files = run_fasttree($msaFile, $alphabet, $model, $output_name);
} else {
    die "Unrecognized program: $program \n";
}
print STDERR "Tree file is $output_files[0]\nLog file is $output_files[1]\n";

sub run_raxml {
    my ($alignment_file, $alphabet, $model, $output_name) = @_;
    print STDERR "In run_raxml (with RELL support), alignment = $alignment_file\n";
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    $model = uc($model); 
    if ($alphabet eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = 'LG' if $model !~ /DAYHOFF|DCMUT|JTT|MTREV|WAG|RTREV|CPREV|VT|BLOSUM62|MTMAM|LG|MTART|MTZOA|PMB|HIVB|HIVW|JTTDCMUT|FLU|STMTREV|DUMMY|DUMMY2|AUTO|LG4M|LG4X|PROT_FILE|GTR_UNLINKED|GTR/i;
        $model = "PROTCAT". $model;
    }

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $parallel);
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $alignment_file);
    push @cmd, ("-n", $output_name);
    push @cmd, ("-f", "D"); # generate RELL support values
    
    my $comment = "command = ". join(" ", @cmd);
    print STDERR "$comment\n\n";
   
    system(@cmd);

    # now map RELL bootstrap replicates support numbers onto ML tree
    @cmd = ("raxmlHPC-PTHREADS-SSE3", "-f", "b", "-m", "GTRCAT"); #to map RELL bootstrap support onto ML tree
    push @cmd, ("-t", "RAxML_bestTree.".$output_name);
    push @cmd, ("-z", "RAxML_rellBootstrap.". $output_name);
    push @cmd, ("-n", $output_name . "_rell");
    my $comment = "command = ". join(" ", @cmd);
    print STDERR "$comment\n\n";
    system(@cmd);
    
    # now rename tree file 
    my $treeFile = $output_name . "_raxml_tree.nwk";
    move("RAxML_bipartitions.".$output_name. "_rell", $treeFile);
    # now concatenate the two info file to the desired log file name
    my $logFile = $output_name . "_raxml_log.txt";
    move("RAxML_info.".$output_name, $logFile);
    system("cat RAxML_info.${output_name}_rell >> $logFile");
    # delete other unwanted files unless $debug
    system("rm RAxML*$output_name") unless $debug;
    unlink("$alignment_file.reduced") if -e "$alignment_file.reduced" and !$debug;
    return ($treeFile, $logFile);
}

sub run_phyml {
    my ($msaFile, $alphabet, $model, $output_name) = @_;
    import Sequence_Alignment;
    my $alignment = new Sequence_Alignment($msaFile);
    my $file_base = $msaFile;
    $file_base =~ s/\..{2,4}$//;
    my $phylip_file = $file_base . ".phy";
    $alignment->write_phylip($phylip_file);

    my $datatype = 'aa';
    if ($alphabet =~ /DNA/i) {
        $datatype = 'nt';
        $model = 'GTR' if $model !~ /HKY85|JC69|K80|F81|F84|TN93|GTR/;
    }
    else {
        $model = 'LG' if $model !~ /WAG|JTT|MtREV|Dayhoff|DCMut|RtREV|CpREV|VT|AB|Blosum62|MtMam|MtArt|HIVw|HIVb/;
    }

    my @cmd = ("phyml");
    push @cmd, ("-i", $phylip_file);
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    push @cmd, ("-b", '-3'); # -1 gives approximate Likelihood Ratio Tests
    # -2 give Chi2-based parametric branch supports
    # -3 gives Shimodaira-Hasegawa (SH) support values
    
    my $comment = "command = ". join(" ", @cmd);
    print STDERR "$comment\n\n";
   
    system(@cmd);
    
    # rename tree file
    my $treeFile = $output_name."_phyml_tree.nwk";
    move($phylip_file."_phyml_tree.txt", $treeFile);
    #rename log file
    my $logFile = $output_name."_phyml_log.txt";
    move($phylip_file."_phyml_stats.txt", $logFile);
    #remove unwanted files unless $debug
    unlink("$phylip_file") unless $debug;
    return ($treeFile, $logFile);
}

sub run_fasttree {
    my ($alignment_file, $alphabet, $model, $output_name) = @_;

    my $treeFile = $output_name."_fasttree_tree.nwk";
    my @cmd = ("FastTree", "-out", $treeFile);
    if ($alphabet =~ /DNA/i) {
        push @cmd, "-nt", "-gtr";
    }
    else {
        $model = 'lg' if $model !~ /LG|WAG/i;
        $model = lc($model);
        push @cmd, "-$model";
    }
    
    push @cmd, $alignment_file;
    my $logFile = $output_name."_fasttree_log.txt";
    push @cmd, "2>", $logFile;
    my $command = join(" ", @cmd);
    print STDERR "command = $command\n\n";
    system($command); # need to run as a string instead of array to get the STDERR redirect to work   
    return ($treeFile, $logFile);
}
