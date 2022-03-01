#
# The GeneTree application.
use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use JSON;
use File::Copy ('copy', 'move');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use URI::Escape;
use Sequence_Alignment; # should be in lib directory
use Phylo_Tree; # should be in lib directory

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;
my $max_genome_length = 100_000; #most/all single-sequence viruses are less than this
my @default_genome_metadata_fields = (
        "species", "strain", "geographic_group", "isolation_country", "host_group", "host_common_name", "collection_year", "subtype", "lineage", "clade");
my @default_feature_metadata_fields = ("product", "accession");

our $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
if ($debug) {
    print STDERR "debug = $debug\n" if $debug;
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}
our @analysis_step => ();# collect info on sequence of analysis steps
our @step_stack => (); # for nesting of child steps within parent steps
my @original_sequence_ids; # list of items requested, can be different from those actually obtained

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#$data_url = "https://patricbrc.org/api" if $debug;
print STDERR "data_url=\n$data_url\n" if $debug;

my $script = Bio::KBase::AppService::AppScript->new(\&build_tree, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";
    my $pf = {
	cpu => 8,
	memory => "128G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub start_step {
    my $name = shift;
    my %step_info = {};
    push @analysis_step, \%step_info;
    $step_info{name} = $name;
    $step_info{comments} = ();
    $step_info{start_time} = time();
    print STDERR "start_step($step_info{name})\n";
    my $stack_depth = scalar @step_stack;
    my $current_step = $step_stack[$stack_depth-1];
    $step_info{parent_step} = $current_step;
    push @step_stack, $name;
    return \@{$step_info{comments}}, \%step_info
}

sub end_step {
    my $name = shift;
    print STDERR "end_step($name)\n";
    my $current_step = pop @step_stack;
    unless ($name eq $current_step) {
        print STDERR "Problem! at end_step name is wrong: $name should be $current_step";
    }
    my $step_index = scalar @analysis_step - 1;
    my $step_info = $analysis_step[$step_index];
    $step_info->{end_time} = time();
}

sub write_report {
    my ($output_file, $tree_graphic_file) = @_;
    print STDERR "write_report()\n";
    open F, ">$output_file";
    print F "<HTML>\n<h1>Gene Tree Output</h1>\n";
    if (-e $tree_graphic_file) {
        my $svg_text = read_file($tree_graphic_file);
        print F $svg_text, "\n\n";
    }
    print F "<h2>Analysis Steps</h2>\n";
    my $start_time = $analysis_step[0]->{start_time};
    my $time_string = localtime($start_time);
    print F "<p>Start time $time_string\n";
    for my $step (@analysis_step) {
        print F "<h3>$step->{name}</h3>\n";
        my $duration = $step->{end_time} - $step->{start_time};
        if (scalar @{$step->{comments}}) {
            #print F "<p>Comments: ", scalar @{$step->{comments}}, "\n<ul>\n";
            print F "<ul>\n";
            for my $comment (@{$step->{comments}}) {
                print F "<li>$comment\n";
            }
            print F "</ul>\n";
        }
        if (exists $step->{details} and $step->{details} =~ /\S/) {
            my $element_name = "$step->{name}_details";
            $element_name =~ tr/ /_/;
            print F "<script>\nfunction toggle_$element_name() {
              var x = document.getElementById(\"$element_name\");
                if (x.style.display == \"none\") {
                      x.style.display = \"block\";
                  } else {
                      x.style.display = \"none\";
                  }
              }\n</script>\n";
            print F "Details: <button onclick=\"toggle_$element_name()\">Show/Hide</button>\n";
            print F "<div id=\"$element_name\" style=\"display:none; background:#f0f0f0\" \n";
            print F "    onclick=\"toggle_$element_name()\">\n";
            print F "<pre>\n", $step->{details}, "\n</pre></div>\n";
        }
        if ($duration > 10) { 
            print F "<p>Duration ", $duration, " seconds.\n";
        }
    }
    print F "</HTML>\n";
}

sub retrieve_sequence_data {
    my ($app, $params, $api,  $feature_metadata_fields, $genome_metadata_fields) = @_;
    my ($step_comments, $step_info) = start_step("Gather Sequence Data");
    my @seq_items; # array of hash refs, one for each dataset
    my $comment;
    my ($aligned_state, $any_in_memory) = (0, 0);
    $aligned_state = (scalar(@{$params->{sequences}}) == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i); 
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    #push @{$step_comments}, "number of input sequence sets = " . scalar @{$params->{sequences}};
    my $total_seqs = 0;
    for my $sequence_item (@{$params->{sequences}}) {
        print STDERR "data item: $sequence_item->{type}, $sequence_item->{filename}\n";
        push @{$step_comments}, "reading $sequence_item->{type} $sequence_item->{filename}";
        my %seq_data; # data about this set of sequences
        push @seq_items, \%seq_data;
        my $num_seqs = 0;
        $seq_data{"type"}=$sequence_item->{type};
        $seq_data{"is_aligned"}=0;
        if ($sequence_item->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            my $local_file = $sequence_item->{filename};
            $local_file =~ s/.*\///;
            print STDERR "About to copy file $sequence_item->{filename} to $local_file\n";
            $app->workspace->download_file($sequence_item->{filename}, $local_file, 1, $global_token);
            open F, $local_file;
            while (<F>) {
                $num_seqs++ if /^>(\S+)/;
                push @original_sequence_ids, $1;
            }
            close F;
            $seq_data{"name"}=$local_file;
            $seq_data{"type"}="local_file"; # makes it easier to recognize which ones are already written to files
            $seq_data{"num_seqs"}=$num_seqs;
            $seq_data{"is_aligned"} = $sequence_item->{type} =~ /ALIGNED/i; # could be aligned_dna_fasta or aligned_protein_fasta
        }
        elsif ($sequence_item->{type} eq "feature_group") {
            # need to get feature sequences from database 
            my $feature_group = $sequence_item->{filename};
            $seq_data{"name"} = $feature_group;
            $seq_data{"type"}="feature_group"; 
            $comment = "retrieving sequences for feature group $feature_group";
            print STDERR $comment, "\n";
            push @{$step_comments}, $comment;
            if ($debug) {
                my $feature_ids = $api->retrieve_patricids_from_feature_group($feature_group);
                print STDERR "\nfeature_ids = @$feature_ids\n";
                print STDERR "number feature_ids is ", scalar @$feature_ids, "\n";
            }

            my @copy_of_feature_metadata_fields = (@$feature_metadata_fields); # use copy because api methods will over-write it
            my $seq_list;
            my $na_or_aa = ('aa', 'na')[$params->{alphabet} eq 'DNA']; 
            $seq_list = $api->retrieve_sequences_from_feature_group($feature_group, $na_or_aa, \@copy_of_feature_metadata_fields);
            if ($debug) {
                print STDERR "retreive_sequences_from_feature_group return is $seq_list\n"; 
                print STDERR "length is ", scalar @$seq_list, "\n"; 
                print STDERR "First element = $seq_list->[0]\n";
                for my $key (sort keys %{$seq_list->[0]}) {
                    print STDERR "    \t$key\t$seq_list->[0]->{$key}\n";
                }
            }
            $seq_data{"seq_list"} = $seq_list;
            $num_seqs = scalar @$seq_list;
            $seq_data{"type"}='feature_group';
            $seq_data{"num_seqs"}= $num_seqs;
            for my $record (@$seq_list) {
                push @original_sequence_ids, $record->{patric_id};
            }
            $comment = "number of elements = $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
        }
        elsif ($sequence_item->{type} eq "genome_group") {
            if (scalar @{$params->{sequences}} > 1) {
                print STDERR "Genome group $sequence_item->{filename} combined with other sequence inputs. This case is not handled yet. Exiting.\n";
                `cd ..`; # to allow temp directory to be deleted
                exit(1);
            }
            $feature_metadata_fields = \(); # we will not be using these
            # need to check genomes are:
            #   1) single-sequence (not multiple contigs)
            #   2) short enough (below threshold length specified in json parameters?)
            $comment = "retrieving sequences for genome group $sequence_item->{filename}\n";
            print STDERR "$comment\n";
            push @{$step_comments}, $comment;
            my $genome_ids = $api->retrieve_patric_ids_from_genome_group($sequence_item->{filename});
            for my $id (@$genome_ids) {
                print STDERR "$id\n";
            }
            push @$genome_metadata_fields, "genome_accessions" unless grep(/genome_accession/, @$genome_metadata_fields);
            my @copy_of_genome_metadata_fields = (@$genome_metadata_fields);
            push @copy_of_genome_metadata_fields, ('genome_id', 'contigs', 'superkingdom', 'genome_length');
            my @genome_metadata = $api->retrieve_genome_metadata($genome_ids, \@copy_of_genome_metadata_fields);
            print STDERR "retrieve_genome_metadata to test for single-sequence virus under $max_genome_length:\n";
            print STDERR join("\t", ('genome_id', 'contigs', 'superkingdom', 'genome_length')), "\n";
            for my $info (@genome_metadata) {
                for my $key ('genome_id', 'contigs', 'superkingdom', 'genome_length') {
                    #for my $key (sort(keys %$info)) {
                    print STDERR "$info->{$key}\t";
                }
                print STDERR "\n";
                if ($info->{'genome_length'} > $max_genome_length) {
                    print STDERR "Problem: length of genome $info->{'genome_id'} exceeds $max_genome_length ($info->{'genome_length'}).\nExiting.\n";
                    `cd ..`;
                    exit(1);
                }
                if ($info->{'superkingdom'} ne 'Viruses') {
                    print STDERR "Problem: superkingdom of genome $info->{'genome_id'} is not 'Viruses' ($info->{'superkingdom'}).\nExiting.\n";
                    `cd ..`;
                    exit(1);
                }
                if ($info->{'contigs'} > 1) {
                    print STDERR "Problem: genome $info->{'genome_id'} has multiple contigs ($info->{'contigs'}).\nExiting.\n";
                    `cd ..`;
                    exit(1);
                }
            }
            print STDERR "All genomes are viruses, all have a single sequence, all are under $max_genome_length bases.\n";
            push @original_sequence_ids, @$genome_ids;
            $num_seqs = scalar @$genome_ids;
            my $fasta_file = $sequence_item->{filename};
            $fasta_file =~ s/.*\///; # remove everything upto and including last slash
            $fasta_file =~ tr/ /_/;  # replace any blank spaces
            $fasta_file .= ".fasta";
            my $genome_seqs = $api->retrieve_contigs_in_genomes_to_temp($genome_ids);
            #print STDERR "retrieve_contigs_in_genomes_to_temp returned $genome_seqs\n";
            print STDERR "saving sequences to $fasta_file\n";
            open OUT, ">$fasta_file";
            open IN, $genome_seqs;
            while (<IN>) {
                s/>accn\|(\d+\.\d+).con.*/>$1/; # replace definition line with simple genome_id
                print OUT $_;
            }
            unlink($genome_seqs);
            close OUT;
            
            $seq_data{"name"}=$fasta_file;
            $seq_data{"type"}="local_file"; # makes it easier to recognize which ones are already written to files
            $seq_data{"num_seqs"}=$num_seqs;
            $seq_data{"genome_metadata"} = \@genome_metadata;
            $comment = "number of genomes = $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
        }
        elsif ($sequence_item->{type} eq "feature_ids") {
            # need to get feature sequences from database 
            $comment = "retrieving sequences for feature ids";
            $comment .= ", this functionality not implemented yet";
            push @{$step_comments}, $comment;
            my $feature_ids = $sequence_item->{sequences};
            exit(1);
            #print STDERR "\tfeature_ids = ", join(", ", @$feature_ids), "\n" if $debug;
            #if ($params->{alphabet} eq 'DNA') {
            #    $sequence_item->{sequences} = $api->retrieve_nucleotide_feature_sequence($feature_ids);
            #}
            #else {
            #    $sequence_item->{sequences} = $api->retrieve_protein_feature_sequence($feature_ids);
            #}
            #$num_seqs = scalar keys %{$sequence_item->{sequences}};
        }
        $comment = "number of sequences retrieved: $num_seqs";
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
        $total_seqs += $num_seqs;
    }
    $comment = $aligned_state ? "sequences are aligned" : "sequences need aligning";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    if ($total_seqs < 4) {
        $comment = "total sequences: $total_seqs, too few to build a tree (minimum 4)";
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
    }
    end_step("Gather Sequence Data");
    return \@seq_items;
}

sub write_sequences_to_file {
    my ($seq_set, $unaligned_file) = @_;
    open F, ">$unaligned_file";
    my $id_field = "";
    if ($seq_set->{type} eq "feature_group") {
        $id_field = "patric_id";
    }
    elsif ($seq_set->{type} eq "genome_group") {
        $id_field = "genome_id";
    }
    else {
        print STDERR "warning: write_sequences_to_file called on a seq_set that is neither feature_group nor genome_group\n";
        print STDERR "name = $seq_set->{name}, type = $seq_set->{type}\n";
        exit(1);
    }
    for my $seq_item (@{$seq_set->{seq_list}}) {
        # try series of possible identifiers 
        my $seq_id = $seq_item->{$id_field};
        $seq_id = $seq_item->{accession} unless $seq_id;
        $seq_id = $seq_item->{feature_id} unless $seq_id;
        if ($seq_id) {
            print F ">$seq_id\n$seq_item->{sequence}\n";
        }
        else {
            print STDERR "Cannot find appropriate sequence id for sequence record:\n";
            for my $key (keys %$seq_item) {
                print "$key\t$seq_item->{$key}\n";
            }
            print STDERR "\n";
        }
    }
    close F;
}

sub merge_sequence_sets_to_file {
    my ($seq_items, $output_file) = @_;
    my ($step_comments, $step_info) = start_step("Merge Sequence Sets");
    my $comment = "merge_sequence_sets_to_file: $seq_items, $output_file";
    print STDERR $comment, "\n" if $debug;
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    # goal is to end up with all sequences written to a single unaligned FASTA file
    # only needed when multiple input sequence sets are passed
    my @seq_list;
    for my $seq_set (@$seq_items) {
        if ($seq_set->{type} eq 'local_file') {
            $comment = "read sequences from $seq_set->{name}";
            #push @{$step_comments}, $comment;
            print STDERR "$comment\n";
            open INDATA, $seq_set->{name};
            my $seqid = '';
            my $sequence = '';
            my $num_seqs = 0;
            while (<INDATA>) {
                if (/^>(\S+)/) {
                    if ($seqid) {
                        push @seq_list, [$seqid, $sequence];
                        $num_seqs++;
                    }
                    $seqid = $1;
                    $sequence = '';
                }
                else {
                    chomp;
                    tr/-//d; # strip out gap characters;
                    $sequence .= $_;
                }
            }
            if ($seqid) {
                push @seq_list, [$seqid, $sequence];
                $num_seqs++;
            }
            $comment = "found $num_seqs sequences in $seq_set->{name}.";
            #push @{$step_comments}, $comment;
            print STDERR "$comment\n";
        }
        elsif ($seq_set->{type} eq 'feature_group') { #handle case of feature_group (sequences in memory)
            my $num_seqs = 0;
            for my $seq_item (@{$seq_set->{seq_list}}) {
                #print "keys to seq_list item as hash: " . join("\n", keys %$seq_item);
                $num_seqs++;
                push @seq_list, [$seq_item->{patric_id}, $seq_item->{sequence}];
            }
            $comment = "found $num_seqs sequences in $seq_set->{name}.";
            #push @{$step_comments}, $comment;
            print STDERR "$comment\n";
        }
    }

    open OUTDATA, ">$output_file";
    my %seq_hash;
    my $numseqs = 0;
    for my $seq_item (@seq_list) {
        my ($seqid, $sequence) = @$seq_item;
        if (exists $seq_hash{$seqid}) { # fix up non-unique seq identifier
            $comment = "got repeated sequence ID $seqid";
            if ($sequence eq $seq_hash{$seqid}) {
                $comment .= " and sequences are identical: dropping second occurrence.";
                continue;
            }
            else {
                my $suffix = 2;
                while (exists $seq_item->{sequences}->{$seqid . "_$suffix"}) {
                    $suffix++
                }
                $seqid = $seqid . "_$suffix";
                $comment .= " but sequences differ, renaming to $seqid.";
            }
            push @{$step_comments}, $comment;
        }
        print OUTDATA ">$seqid\n$sequence\n";
        $seq_hash{$seqid} = $sequence;
        $numseqs++;
    }
    close OUTDATA;
    $comment = "$numseqs sequences written to $output_file";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    end_step("Merge Sequence Sets");
}

sub build_tree {
    my ($app, $app_def, $raw_params, $params) = @_;

    my $time1 = `date`;
    print "Proc GeneTree build_tree ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my @outputs; # array of tuples of (filename, filetype)
    my $api = P3DataAPI->new;
    my $tmpdir = File::Temp->newdir( "/tmp/GeneTree_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
    my $original_wd = getcwd();
    chdir($tmpdir); # do all work in temporary directory
   
    my @feature_metadata_fields = @default_feature_metadata_fields;
    if (exists $params->{feature_metadata_fields}) {
        @feature_metadata_fields = @{$params->{feature_metadata_fields}};
    }
    #ensure that patric_id and genome_id are retrieved
    push @feature_metadata_fields, "patric_id" unless grep(/patric_id/, @feature_metadata_fields);
    push @feature_metadata_fields, "genome_id" unless grep(/genome_id/, @feature_metadata_fields);

    my @genome_metadata_fields = @default_genome_metadata_fields;
    if (exists $params->{genome_metadata_fields}) {
        @genome_metadata_fields = @{$params->{genome_metadata_fields}};
    }

    my $seq_items = retrieve_sequence_data($app, $params, $api, \@feature_metadata_fields, \@genome_metadata_fields);
    my $unaligned_file = undef;
    my $aligned_file = undef;
    my $seqids_are_genome_ids = 0; # indicate whether seq identifiers are links to BVBRC genomes in database
    if (scalar @$seq_items == 1) {
        if ($seq_items->[0]->{type} eq 'local_file') {
            if ($seq_items->[0]->{is_aligned}) {
                $aligned_file = $seq_items->[0]->{name};
            }
            else {
                $unaligned_file = $seq_items->[0]->{name};
            }
            $seqids_are_genome_ids = 1 if exists $seq_items->[0]->{genome_metadata}; 
            print STDERR "seqids_are_genome_ids = $seqids_are_genome_ids\n";
        }
        elsif ($seq_items->[0]->{type} eq 'feature_group' | $seq_items->[0]->{type} eq 'genome_group') {
            $unaligned_file = "$params->{output_file}_unaligned.fa";
            write_sequences_to_file($seq_items->[0], $unaligned_file);
        }
    }
    else {
        $unaligned_file = "$params->{output_file}_unaligned.fa";
        merge_sequence_sets_to_file($seq_items, $unaligned_file);
    }
    my $num_seqs = 0;
    for my $seq_set (@$seq_items) {
        $num_seqs += $seq_set->{num_seqs};
    }
    if ($num_seqs < 4) { #need at least 4 seuqences to build a tree
        print STDERR "After retrieval, number of sequences is $num_seqs, less than 4. Cannot build a tree.\n";
        #die "Too few sequences."
        exit(1);
    }
    my $alignment_modified = 0; # flag whether alignment is different from input data 
    if ($unaligned_file) {
        $aligned_file = $unaligned_file;
        $aligned_file =~ s/\..{1,6}$//;
        $aligned_file =~ s/unaligned/aligned/;
        $aligned_file .= ".afa";
        #run_muscle($unaligned_file, $aligned_file);
        run_mafft($unaligned_file, $aligned_file);
        $alignment_modified = 1;
    }
    if ($params->{trim_threshold} or $params->{gap_threshold})
    {
        my $trimmed_aligned_file = $aligned_file;
        $trimmed_aligned_file =~ s/.afa//;
        $trimmed_aligned_file .= "_trimmed.afa";
        my $retval = trim_alignment($aligned_file, $trimmed_aligned_file, $params->{trim_threshold}, $params->{gap_threshold});
        if ($retval) {
            $alignment_modified = 1;
            print STDERR "trimmed aligned file written to $trimmed_aligned_file\n" if $debug;
            $aligned_file = $trimmed_aligned_file;
        }
    }
    my $alphabet = $params->{alphabet};
    push @outputs, [$aligned_file, "aligned_${alphabet}_fasta"] if $alignment_modified;
    open F, $aligned_file;
    $num_seqs = 0;
    my %altered_name;
    my $alignment_text;
    my @seqid_list;
    my $seqids_are_patric_ids = 0;
    while (<F>) { # count sequences and replace any illegal characters in sequence IDs (not allowed by newick standard)
        if (/^>(\S+)/) {
            $num_seqs++;
            my $seq_id = $1;
            if ($seq_id =~ /[[]():]/) {
                my $orig = $seq_id;
                $seq_id =~ tr/:()[]/_____/; # replace any bad characters with underscores
                $altered_name{$seq_id} = $orig;
                print STDERR "replacing identifier $orig with $seq_id\n";
            }
            $seqids_are_patric_ids = 1 if $seq_id =~ /fig\|(\d+)/; # note that seqids may be valid for retrieving feature metadata
            push @seqid_list, $seq_id;
            $alignment_text .= ">$seq_id\n";
        }
        else {
            $alignment_text .= $_;
        }
    }
    if ($num_seqs < 4) { #need at least 4 seuqences to build a tree
        print STDERR "After looking at aligned seqs, number of sequences is $num_seqs, less than 4. Cannot build a tree.\n";
        exit(1);
    }
    if (scalar keys %altered_name) {  # if sequence IDs needed to be altered, write altered alignment
        my $altered_aligned_file = $aligned_file;
        $altered_aligned_file =~ s/\..{3,6}$//; #remove extension
        $altered_aligned_file .= "_altered_ids.afa";
        $aligned_file = $altered_aligned_file;
        write_file($aligned_file, $alignment_text);
        #open F, ">$aligned_file";
        #print F $alignment_text;
        #close F;
    }
    run("echo $tmpdir && ls -ltr $tmpdir") if $debug;

    my $model = "LG"; # default for protein
    if ($params->{substitution_model}) {
        $model = $params->{substitution_model}
    }
    elsif ($params->{alphabet} =~ /DNA/i) {
        $model = "GTR"
    }
    my $output_name = $params->{output_file};
    my $recipe = "raxml"; #default
    if (defined $params->{recipe} and $params->{recipe}) {
        $recipe = lc($params->{recipe})
    }
    print STDERR "About to call tree program $recipe\n";
    my @tree_outputs;
    if ($recipe eq 'raxml') {
        @tree_outputs = run_raxml($aligned_file, $alphabet, $model, $output_name);
    } elsif ($recipe eq 'phyml') {
        my $alignment = new Sequence_Alignment($aligned_file);
        my $phylip_file = $aligned_file;
        $phylip_file =~ s/\..{2,4}$//;
        $phylip_file .= ".phy";
        $alignment->write_phylip($phylip_file);
        @tree_outputs = run_phyml($phylip_file, $alphabet, $model, $output_name);
    } elsif (lc($recipe) eq 'fasttree') {
        #$alignment->write_fasta($aligned_file);
        @tree_outputs = run_fasttree($aligned_file, $alphabet, $model, $output_name);
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    push @outputs, @tree_outputs;
    # optionally re-write newick file with original sequence IDs if any were altered
    # not sure we want to do that as it could screw up downstream use of newick data i
    # (ideally, bad characters should be escaped, but does raxml know about escaped characters?)
    # generate tree graphic using figtree
    my $tree_file = undef;
    my $tree_graphic = undef;
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        if ($type eq 'nwk') {
            $tree_file = $ofile;
            $tree_graphic = generate_tree_graphic($tree_file, $num_seqs);
            push @outputs, $tree_graphic;
        }
    }
    print STDERR "tree_file $tree_file\nare_patric_ids = $seqids_are_patric_ids\nare_genome_ids = $seqids_are_genome_ids\n" if $debug;
    if ($tree_file and ($seqids_are_patric_ids or $seqids_are_genome_ids)) # avoid looking for metadata if ids don't link to database
    {  
        my $metadata = gather_metadata($seq_items, \@seqid_list, \%altered_name, \@feature_metadata_fields, \@genome_metadata_fields); 
        print STDERR "gather_metadata: hash size = " . scalar(keys %$metadata). "\n" if $debug;
        my $tree = new Phylo_Tree($tree_file);
        if ($metadata and scalar keys %$metadata) {
            my $metadata_file = "$params->{output_file}_metadata.txt";
            open F, ">$metadata_file";
            my @header_fields = ();
            for my $field (keys %$metadata) {
                push @header_fields, $field unless $field eq 'patric_id'; # dont put patric_id in header fields
            }
            print F join("\t", "SeqId", @header_fields) . "\n";
            if ($debug) {
                print STDERR "seqid_list: ", join(' ', @seqid_list), "\n";
                for my $field (@header_fields) {
                    print STDERR "For md $field: seq keys = ", join(" ", sort keys %{$metadata->{$field}}), "\n";
                }
            }
            for my $id (@seqid_list) {
                print F $id;
                for my $field (@header_fields) {
                    my $val = "na";
                    $val = $metadata->{$field}{$id} if exists $metadata->{$field}{$id};
                    print F "\t$val";
                }
                print F "\n";
            }
            close F;
            push @outputs, [$metadata_file, "tsv"];
            
            for my $field (@feature_metadata_fields, @genome_metadata_fields) {
                print STDERR "add metadata field: $field\n";
                $tree->add_tip_phyloxml_properties($metadata->{$field}, $field, "BVBRC");
            }
        }
        my ($step_comments, $step_info) = start_step("Write PhyloXML");
        my $phyloxml_data = $tree->write_phyloXML();
        my $phyloxml_file = "$params->{output_file}.xml";
        write_file($phyloxml_file, $phyloxml_data);
        push @outputs, [$phyloxml_file, "phyloxml"];
        end_step("Write PhyloXML");

        if (exists $params->{relabel_tree_fields}) {
            my @relabel_fields;
            for my $field (@{$params->{relabel_tree_fields}}) { #"species", "product") {
                push @relabel_fields, $field if (exists $metadata->{$field});
            }
            if (scalar @relabel_fields) {
                print STDERR "Now write a tree with ids substituted with metadata fields: ", join(", ", @relabel_fields), "\n";
                my $relabeled_newick_file = label_tree_with_metadata($tree_file, $metadata, \@relabel_fields);
                push @outputs, [$relabeled_newick_file, 'nwk'];
            }
        }
    }
    
    my $html_file = "$params->{output_file}_gene_tree_report.html";
    write_report($html_file, $tree_graphic);
    push @outputs, [$html_file, "html"];

    print STDERR '\@outputs = '. Dumper(\@outputs);
    my $output_folder = $app->result_folder();
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n" if $debug;
        if (0) { # for some reason this doesn't work
           $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
               (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
               $global_token);
        }
        else { # fall back to calling CLI
            my $ext = $1 if $ofile =~ /.*\.(\S+)$/;
            my @cmd = ("p3-cp", "-f", "-m", "${ext}=$type", $ofile, "ws:" . $app->result_folder);
            print STDERR "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        }
    }
    chdir($original_wd); # change back to the starting working directory
    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");
}

sub gather_metadata {
    my ($seq_items, $seqid_ar, $altered_names_hr, $feature_fields, $genome_fields) = @_; 
    my ($step_comments, $step_info) = start_step("Gather Metadata");
    my $comment = "feature ids: ". join(", ", @{$seqid_ar});
    push @{$step_comments}, $comment;
    print STDERR "$comment\n"; 
    $comment = "feature metadata fields: ". join(", ", @{$feature_fields});
    push @{$step_comments}, $comment;
    print STDERR "$comment\n"; 
    my $comment = "in gather_metadata: feature ids= $seqid_ar\n";
    print STDERR $comment;
    if (exists $seq_items->[0]->{genome_metadata}) {
        # case of aligning whole viral genomes, we already have the genome metadata, no features as such
        print STDERR "Genome metadata in hand, reformat and return.\n" if $debug;
        print STDERR "Genome fields: @$genome_fields\n";
        my %genome_metadata;
        for my $field (@$genome_fields) {
            $genome_metadata{$field} = ();
            print STDERR "\ngenmetrec[0] for $field: $seq_items->[0]->{genome_metadata}->[0]->{$field}\n" if $debug;
            for my $record (@{$seq_items->[0]->{genome_metadata}}) {
                my $genome_id = $record->{genome_id};
                $genome_metadata{$field}{$genome_id} = $record->{$field};
                #print STDERR "$field\t$genome_id\t$genome_metadata{$field}{$genome_id}\n" if $debug;
            }
        }
        end_step("Gather Metadata");
        return \%genome_metadata;
    }
    my $feature_metadata = get_feature_metadata($seq_items, $seqid_ar, $altered_names_hr, $feature_fields);
    #get_feature_metadata($feature_metadata, $seq_ids_ar, $feature_fields);
    $comment = "number of feature metadata fields: " . scalar keys %{$feature_metadata};
    push @{$step_comments}, $comment;
    print STDERR "$comment\n"; 
    if ($feature_metadata and scalar keys %{$feature_metadata}) {
        my %genome_ids = ();
        for my $feature_id (keys %{$feature_metadata->{genome_id}}) {
            if (exists $feature_metadata->{genome_id}{$feature_id}) {
                my $genome_id = $feature_metadata->{genome_id}{$feature_id};
                $genome_ids{$genome_id} = 1;
            }
        }
        my @genome_ids = keys %genome_ids;
        $comment = "genome ids: ". join(", ", @genome_ids);
        push @{$step_comments}, $comment;
        print STDERR "$comment\n"; 
        $comment = "genome metadata fields: ". join(", ", @{$genome_fields});
        push @{$step_comments}, $comment;
        print STDERR "$comment\n"; 
        print STDERR "before merging genome_metadata, feature_metadata keys = \n", join(", ", keys %$feature_metadata), "\n";
        my @tmp = keys %$feature_metadata;
        print STDERR "tmp[0] = $tmp[0]\n";
        print STDERR "for $tmp[0]: ", join(", ", keys %{$feature_metadata->{$tmp[0]}}), "\n";

        if (scalar @genome_ids) {
            for my $genome_field (@$genome_fields) {
                $feature_metadata->{$genome_field} = ();
            }
            my $genome_metadata = get_genome_metadata(\@genome_ids, $genome_fields);
            for my $feature_id (keys %{$feature_metadata->{genome_id}}) {
                my $genome_id = $feature_metadata->{genome_id}{$feature_id};
                for my $genome_field (@$genome_fields) {
                    my $value = $genome_metadata->{$genome_id}{$genome_field};
                    $feature_metadata->{$genome_field}{$feature_id} = $value;
                    print STDERR "gm $feature_id $genome_id $genome_field $value\n" if $debug;
                }
            }
        }
        #my @fields = @{$feature_fields};
        #for my $genome_field (@{$genome_fields}) {
        #    push(@fields, $genome_field) unless grep(/^$genome_field\$/, @fields);
        #}
    }
    end_step("Gather Metadata");
    return $feature_metadata;
}

sub trim_alignment {
    my ($aligned_file, $trimmed_aligned_file, $trim_threshold, $gap_threshold) = @_;
    my ($trim_comments, $trim_info) = start_step("Trim Alignment");
    my $comment = "performing trimming on alignment: trim_threshod=$trim_threshold, gap_threhold=$gap_threshold";
    print STDERR "$comment\n";
    push @{$trim_comments}, $comment;
    my $alignment = new Sequence_Alignment($aligned_file);
    print STDERR "trim_alignment: alignment=$alignment\n" if $debug;
    $trim_info->{details} = "Before trimming:\n" .  $alignment->write_stats();
    my $alignment_changed = 0;
    if ($trim_threshold > 0) {
        $comment = "trim ends of alignment to density $trim_threshold";
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        my ($left_trim_cols, $right_trim_cols) = $alignment->end_trim($trim_threshold);
        if ($left_trim_cols or $right_trim_cols) {
            $comment = "ends trimmed: $left_trim_cols columns on left, $right_trim_cols columns on right.";
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment_changed = 1;
        }
        else {
            $comment = "no end-columns needed to be trimmed";
            push @{$trim_comments}, $comment;
            print STDERR "$comment\n";
        }
    }
    if ($gap_threshold > 0) {
        $comment = "delete any sequences with gap proportion greater than $gap_threshold";
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        my $deleted_seqs = $alignment->delete_gappy_seqs($gap_threshold);
        if (scalar @$deleted_seqs) {
            $comment = scalar(@$deleted_seqs) . " gappy sequences deleted: ". join(", ", @$deleted_seqs);
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment_changed = 1;
        }
        else {
            $comment = "no sequences needed to be deleted";
            push @{$trim_comments}, $comment;
            print STDERR "$comment\n";
        }
    }
    if ($alignment_changed) {
        $trim_info->{details} .= "After trimming:\n" . $alignment->write_stats();
    }
    my $total_seqs = $alignment->get_ntaxa();
    if ($total_seqs < 4) {
        $comment = "total sequences: $total_seqs, too few to build a tree (minimum 4)";
        push @{$trim_comments}, $comment;
        print STDERR "$comment\n";
    }
    if ($alignment_changed) {
        $alignment->write_fasta($trimmed_aligned_file); 
        $comment = "writing trimmed alignment to $trimmed_aligned_file\n";
        push @{$trim_comments}, $comment;
        print STDERR $comment;
    }
    end_step("Trim Alignment");
    return $alignment_changed;
}

sub run_raxml {
    my ($alignment_file, $alphabet, $model, $output_name) = @_;
    my ($step_comments, $step_info) = start_step("Phylogenetic Inference with RAxML");
    print STDERR "In run_raxml (with RELL support), alignment = $alignment_file\n";
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
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
    push @cmd, ("-s", basename($alignment_file));
    push @cmd, ("-n", $output_name);
    push @cmd, ("-f", "D"); # generate RELL support values
    
    my $comment = "command = ". join(" ", @cmd);
    push @{$step_comments}, $comment;
    print STDERR "$comment\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    $step_info->{details} = $out;

    # now map RELL bootstrap replicates support numbers onto ML tree
    @cmd = ("raxmlHPC-PTHREADS-SSE3", "-f", "b", "-m", "GTRCAT"); #to map RELL bootstrap support onto ML tree
    push @cmd, ("-t", "RAxML_bestTree.".$output_name);
    push @cmd, ("-z", "RAxML_rellBootstrap.". $output_name);
    push @cmd, ("-n", $output_name . "_rell");
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    $step_info->{details} .= $out;
    
    my @outputs;
    my $treeFile = $output_name . "_RAxML_tree_rell.nwk";
    move("RAxML_bipartitions.".$output_name. "_rell", $treeFile);
    my $logFile = $output_name . "_RAxML_log.txt";
    move("RAxML_info.".$output_name, $logFile);
    #my $details = read_file($logFile);
    #$step_info->{details} .= $details;
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");
    end_step("Phylogenetic Inference with RAxML");
    return @outputs;
}

sub run_phyml {
    my ($alignment_file, $alphabet, $model, $output_name) = @_;
    my ($step_comments, $step_info) = start_step("Phylogenetic Inference with Phyml");

    my $cwd = getcwd();
    
    my $datatype = 'aa';
    if ($alphabet =~ /DNA/i) {
        $datatype = 'nt';
        $model = 'GTR' if $model !~ /HKY85|JC69|K80|F81|F84|TN93|GTR/;
    }
    else {
        $model = 'LG' if $model !~ /WAG|JTT|MtREV|Dayhoff|DCMut|RtREV|CpREV|VT|AB|Blosum62|MtMam|MtArt|HIVw|HIVb/;
    }

    my @cmd = ("phyml");
    push @cmd, ("-i", $alignment_file);
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    
    my $comment = "command = ". join(" ", @cmd);
    push @{$step_comments}, $comment;
    print STDERR "$comment\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    $step_info->{details} = $out;
    #$step_info->{stderr} = $err;
    #my $comment = "return code = $rc\n";
    #push @{$step_comments}, $comment;
    
    my @outputs;
    my $treeFile = $output_name."_phyml_tree.nwk";
    move($alignment_file."_phyml_tree.txt", $treeFile);
    my $logFile = $output_name."_phyml_log.txt";
    move($alignment_file."_phyml_stats.txt", $logFile);
    my $details = read_file($logFile);
    $step_info->{details} .= $details;
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");

    end_step("Phylogenetic Inference with Phyml");
    return @outputs;
}

sub run_fasttree {
    my ($alignment_file, $alphabet, $model, $output_name) = @_;
    my ($step_comments, $step_info) = start_step("Phylogenetic Inference with FastTree");

    my $treeFile = $output_name."_fasttree.nwk";
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
    my $comment = join(" ", @cmd);
    push @{$step_comments}, $comment; 
    print STDERR $comment . "\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";

    $step_info->{details} = $err;
    
    my @outputs;
    my $logFile = $output_name."_fasttree_log.txt";
    write_file($logFile, $err);
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");

    end_step("Phylogenetic Inference with FastTree");
    return @outputs;
}

sub get_feature_metadata {
    # add/supplement any data in feature_metadata hashref from patric feature id to hash of field values
    my ($seq_items, $seqid_ar, $altered_names_hr, $feature_fields) = @_;
    $feature_fields = ['genome_id','product'] unless $feature_fields;
    my %feature_metadata;
    for my $field (@$feature_fields) {
        $feature_metadata{$field} = (); # column major order
    }
    my %seqs_with_data;
    for my $seq_set (@$seq_items) {
        if ($seq_set->{type} eq 'feature_group') {
            for my $feature_data (@{$seq_set->{seq_list}}) {
                my $seq_id = $feature_data->{patric_id};
                next unless $seq_id;
                for my $field (@$feature_fields) { # get just specified feature fields
                    $feature_metadata{$field}{$seq_id} = $feature_data->{$field} if exists $feature_data->{$field};
                }
                $seqs_with_data{$seq_id} = 1;
            }
        }
    }

    my @seqs_needing_data;
    my %alt_name_rev;
    for my $seq_id (@$seqid_ar) {
        next if exists $seqs_with_data{$seq_id};
        if (exists $altered_names_hr->{$seq_id}) {
            my $alt_id = $altered_names_hr->{$seq_id}; 
            $alt_name_rev{$alt_id} = $seq_id;
            $seq_id = $alt_id;
        }
        push @seqs_needing_data, $seq_id;
    }
    if (scalar @seqs_needing_data) {
        my $select_string = "select(" . join(",", @$feature_fields) . ")"; 
        my $escaped_ids = join(",", map { uri_escape $_ } @seqs_needing_data);
        my $url = "$data_url/genome_feature/?in(patric_id,($escaped_ids))&$select_string";
        print STDERR "query=$url\n";
        my $resp = curl_json($url);
        print STDERR "response=$resp\n";
        if ( $resp =~ /Error/) {
            warn "Problem! query did not return properly: $url\n";
            return undef;
        }
        print STDERR "length of resp = ", scalar(@$resp), "\n";
        #***** need to fix this part
        for my $member (@$resp) {
            #my $feature_id = $member->{patric_id};
            #$feature_metadata{$feature_id} = $member;
            print STDERR " next member = $member\n";
            my $seq_id = $member->{patric_id};
            for my $field (@$feature_fields) {
                $feature_metadata{$field}{$seq_id} = $member->{$field}
            }
        }
    }
    return \%feature_metadata;
}

sub get_genome_metadata {
    my ($genome_ids, $fields) = @_;
    print STDERR "in get_genome_metadata: genome ids = ", join(", ", @$genome_ids), "\n" if $debug;
    return \() unless scalar(@$genome_ids);
    $fields = ['species'] unless $fields;
    my $select_string = "select(" . join(",", 'genome_id', @$fields) . ")"; 
    my $escaped_ids = join(",", @$genome_ids);
    my $url = "$data_url/genome/?in(genome_id,($escaped_ids))&$select_string";
    print STDERR "query=$url\n";
    my $resp = curl_json($url);
    my %genome_metadata = ();
    for my $member (@$resp) {
        print STDERR "member: ", join(", ", sort(keys %$member)), "\n";
        my $genome_id = $member->{genome_id};
        for my $field (@$fields) {
            $genome_metadata{$genome_id}{$field} = $member->{$field} unless $field eq 'genome_id';
        }
    }
    return \%genome_metadata
}

sub label_tree_with_metadata {
    my ($input_newick, $metadata, $label_fields) = @_;
    my ($step_comments, $step_info) = start_step("Label Tree With Metadata");
    my $comment = "label_tree_with_metadata called for fields: " . join(", ", @$label_fields);
    for my $field (@$label_fields) {
        print STDERR "metadata does not include $field\n" unless exists $metadata->{$field};
    }
    print STDERR $comment, "\n";
    my $output_newick = $input_newick;
    $output_newick =~ s/\..{2,6}$//; #remove extension
    $output_newick .= "_relabeled_" . join("_", @$label_fields) . ".nwk";
    my $newick = read_file($input_newick);
    print STDERR "Got newick string: ", substr($newick, 0, 50), "\n";
    my $start_field = $label_fields->[0];
    for my $seq_id (keys %{$metadata->{$start_field}}) {
        my @values = ();
        for my $field (@$label_fields) {
            push @values, $metadata->{$field}{$seq_id};
        }
        my $sub = join("|", @values);
        $seq_id =~ s/\|/\\\|/g;
        $sub =~ tr/[](),://d; #remove characters that break newick structure
        $sub =~ tr/ /_/; #remove characters that break newick structure
        print STDERR "substitute $seq_id with $sub\n";
        $newick =~ s/$seq_id/$sub/;
    }
    print STDERR "Modified newick string: ", substr($newick, 0, 50), "\n";
    write_file($output_newick, $newick);
    end_step("Label Tree With Metadata");
    return $output_newick;
}

sub generate_tree_graphic {
    my ($input_newick, $num_tips) = @_;
    my ($step_comments, $step_info) = start_step("Generate Tree Graphic");
    my $file_base = basename($input_newick);
    $file_base =~ s/\..{2,6}//;
    my $tree_graphic_file = "$file_base.svg";
    my $nexus_file = "$file_base.nex";
    my $comment = "run figtree input = $input_newick, output = $tree_graphic_file";
    #push @{$step_comments}, $comment;
    print STDERR "$comment\n";

    open F, ">$nexus_file";
    print F "#NEXUS\nbegin trees;\n";
    my $tree_text = read_file($input_newick);
    print F "tree one = [&U] $tree_text\nend;\n\n";
    print F "begin figtree;\n";
    print F "set appearance.branchLineWidth=3.0;\n";
    print F "set tipLabels.fontSize=14;\n";
    print F "set tipLabels.fontName=\"sansserif\";\n";
    print F 'set trees.order=true;';
    print F 'set trees.orderType="increasing";';
    print F 'set trees.rooting=true;';
    print F 'set trees.rootingType="Midpoint";';
    print F "end;\n";
    close F;

    my @cmd = ("figtree", "-graphic", 'SVG');
    
    if ($num_tips > 40) {
        my $height = 600 + 15 * ($num_tips - 40); # this is an empirical correction factor to avoid taxon name overlap
        push @cmd, '-height', $height;
    }
    push @cmd, $nexus_file, $tree_graphic_file;
    $comment = join(" ", @cmd);
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";

    my ($stdout, $stderr) =  run_cmd(\@cmd);
    $step_info->{stdout} = $stdout;
    $step_info->{stderr} = $stderr;
    end_step("Generate Tree Graphic");
    return [$tree_graphic_file, "svg"];
}

sub run_muscle {
    my ($unaligned, $aligned) = @_;
    print STDERR "run_muscle($unaligned, $aligned)\n";
    my ($step_comments, $step_info) = start_step("Align with muscle");
    my $cmd = ["muscle", "-in", $unaligned, "-out", $aligned];
    my $comment = "command = " . join(" ", @$cmd);
    push @{$step_comments}, $comment;
    my ($stdout, $stderr) =  run_cmd($cmd);
    $step_info->{details} = $stderr;
    end_step("Align with muscle");
}

sub run_mafft {
    my ($unaligned, $aligned) = @_;
    print STDERR "run_mafft($unaligned, $aligned)\n";
    my ($step_comments, $step_info) = start_step("Align with mafft");
    my $cmd = ["mafft", "--auto", $unaligned];
    my $comment = "command = " . join("\t", @$cmd);
    print STDERR "$comment\n";
    my $out;
    open $out, ">$aligned";
    run($cmd, ">", $out);
    end_step("Align with mafft");
}

sub curl_text {
    my ($url) = @_;
    my @headers = ("-H", "Accept:text/tsv");
    if ($global_token) { push @headers, "-H", "Authorization:$global_token"}
    my @cmd = ("curl", $url, @headers);
    print STDERR join(" ", @cmd), "\n" if $debug;
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my @headers = ("-H", "Accept:application/json");
    if ($global_token) { push @headers, "-H", "Authorization:$global_token"}
    my @cmd = ("curl", $url, @headers);
    print STDERR join(" ", @cmd), "\n" if $debug;
    my ($out) = run_cmd(\@cmd);
    print STDERR $out, "\n" if $debug;
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = $global_token;
    push(@opts, "-H", "Authorization:$token");
    push(@opts, "-H", "Accept:text/tsv");
    #push(@opts, "-H", "Content-Type: multipart/form-data");
    return @opts;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}

sub get_ws {
    return $global_ws;
}

sub get_token {
    return $global_token;
}

sub get_ws_file {
    my ($id) = @_;
    # return $id; # testing
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = $base;
    # return $file; # testing
    print STDERR "get_ws_file:  base=$base, file=$file, ws=$ws, token=$token\n";
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $base $id\n";
    system("ls -lrta ");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
        die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -lrta ");

    return $file;
}

sub write_output {
    my ($string, $ofile) = @_;
    open(F, ">$ofile") or die "Could not open $ofile";
    print F $string;
    close(F);
}

sub verify_cmd {
    my ($cmd) = @_;
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}
