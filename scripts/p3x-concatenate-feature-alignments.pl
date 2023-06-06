use strict;
use P3Utils;

my($opt, $usage) = P3Utils::script_opts('alignments',
                ['trimsoftmask|t', 'Trim any sequence in lowercase'],
                ['verbose|debug|v', 'Write status messages to STDERR.'],
            );
if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}

my @families = @ARGV;

my %genomes;
my %family_genome_count;
my %family_genome_features;
my %feature_seq;

for my $family_file (@families) {
    # read aligned fasta file and update genomes and family_genome_count
    my $family = $family_file;
    $family =~ s/\..{1,5}//;
    print "reading alignment for $family from $family_file\n";
    open my $f, $family_file;
    my $feature;
    while (<$f>) {
        if (/^>(\S+)/) {
            $feature = $1;
            my $genome_id = $1;
            if ($feature =~ /fig\|(\d+\.\d+)/ or $feature =~ /^(\d+\.\d+)$/ or /::(\d+\.\d+)/)
            {
                $genome_id = $1;
            }
            else {
                />\S+\s+(\S+)/;
                $feature = $genome_id . "_$1"; 
                print "cannot parse genome_id, using $genome_id : $feature";
                
            }
            $genomes{$genome_id}++;
            $family_genome_count{$family}{$genome_id}++;
            push @{$family_genome_features{$family}{$genome_id}}, $feature;
            $feature_seq{$feature} = "";
        }
        else {
            chomp;
            if ($opt->trimsoftmask) {
                s/[a-z]//g;
            }
            $feature_seq{$feature} .= $_;
        }
    }
}

my %family_alignment_length;
for my $family (keys %family_genome_features) {
    print "Study $family\n";
    for my $genome (keys %{$family_genome_features{$family}}) {
        for my $feature (@{$family_genome_features{$family}{$genome}}) {
            my $len = length($feature_seq{$feature});
            $family_alignment_length{$family} = $len unless exists $family_alignment_length{$family};
            die "aligned sequence $feature of length $len, mismatching expectation of $family_alignment_length{$family} for $family" unless $len == $family_alignment_length{$family};
        }
    }
}
my $outfile = "concatenated_alignments_" . scalar(keys %genomes). "_genomes_" . scalar(keys %family_genome_features). "_genes.afa";
print "outfile = $outfile\n";

open my $f, ">", $outfile;
my @genomes = sort keys %genomes;
for my $genome (@genomes) {
    print "Study $genome\n";
    print $f ">$genome\n";
    for my $family (keys %family_genome_features) {
        if (exists $family_genome_features{$family}{$genome}) {
            my $feature = @{$family_genome_features{$family}{$genome}}[0]; # just get the first one -- better option is to see which is best fit to the rest of the alignment
            print $f $feature_seq{$feature};
        }
        else {
            my $gapseq = '-' * $family_alignment_length{$family};
            print $f $gapseq;
        }
    }
    print $f "\n";
}
