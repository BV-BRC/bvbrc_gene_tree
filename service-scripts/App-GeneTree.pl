#
# The GeneTree application.
#

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

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

my $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
print STDERR "args = ", join("\n", @ARGV), "\n" if $debug;
my %details = ("step" => [], "step_stack"= []); # collect info on sequence of genes

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#my $data_url = "http://www.alpha.patricbrc.org/api";

my $script = Bio::KBase::AppService::AppScript->new(\&build_tree, \&preflight);
my $rc = $script->run(\@ARGV);

my $input_sequences = undef;
my $tmpdir = undef; 

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    $tmpdir = File::Temp->newdir( "/tmp/GeneTree_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "$tmpdir\n";
    my ($all_sequences, $aligned) = retrieve_sequence_data($app, $params);

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
    my %step_details = {name => $name, comments => ()}
    push @{$details->{step}}, \%step_details;
    $step_details->{start_time} = time();
    my $stack_depth = scalar @{$details->{step_stack}};
    my $current_step = $details->{step_stack}[$stack_depth-1];
    $step_details->{parent_step} = $current_step;
    push @{$details->{step_stack}}, $name;
    return \$step_details
}

sub end_step {
    my $name = shift;
    my $current_step = pop @{$details->{step_stack}};
    assert($name eq $current_step);
    my $step_index = scalar @{$details->{step}};
    my $step_details = $details->{step}[$step_index];
    $step_details->{end_time} = time();
}

sub retrieve_sequence_data {
    my ($app, $params) = @_;
    my $api = P3DataAPI->new;
    my $step_details = start_step("retrieve_sequence_data");
    #print STDERR "params = $params\n";
    #$param_sequences = @{$param_sequences};
    my $aligned = (scalar(@{$params->{sequences}}) == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i); 
    my %all_sequences = {};
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    push @{$step_details->{comments}}, "number of input sequence sets = " . scalar @{$params->{sequences}};
    for my $sequence_item (@{$params->{sequences}}) {
        print STDERR "data item: $sequence_item->{type}, $sequence_item->{filename}\n";
        push @{$step_details->{comments}}, "reading $sequence_item->{type} $sequence_item->{filename}";
        my %item_sequences = {};
        if ($sequence_item->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            $app->workspace->download_file($sequence_item->{filename}, "$tmpdir/temp_seq_file.fa", 1, $global_token);
            open INDATA, "$tmpdir/temp_seq_file.fa";
            my $seqid = '';
            while (<INDATA>) {
                if (/^>(\S+)/) {
                    $seqid = $1;
                    if (exists $item_sequences{$seqid}) {
                        my $comment = "got non-unique sequence ID $seqid, appending suffix";
                        push @{$step_details->{comments}}, $comment;
                        my $suffix = 2;
                        while (exists $item_sequences{$seqid . "_$suffix"}) {
                            $suffix++
                        }
                        $seqid = $seqid . "_$suffix"
                    }
                }
                else {
                    chomp;
                    $item_sequences->{$seqid} .= $_
                }
            }

        }
        elsif ($sequence_item->{type} eq "feature_group" or $sequence_item->{type} eq "feature_ids") {
                # need to get feature sequences from database 
            my $comment = "retrieving sequences for feature group $sequence_item->{filename}";
            push @{$step_details->{comments}}, $comment;
            my $feature_ids;
            if ($sequence_item->{type} eq 'feature_group') {
                print STDERR "retrieving ids for feature_group: $sequence_item->{filename}\n" if $debug;
                $feature_ids = $api->retrieve_patricids_from_feature_group($sequence_item->{filename});
            }
            else {
                $feature_ids = $sequence_item->{sequences};
            }
            print STDERR "feature_ids = ", join(", ", @$feature_ids), "\n" if $debug;
            if ($params->{alphabet} eq 'DNA') {
                %item_sequences = %$api->retrieve_nucleotide_feature_sequence($feature_ids);
            }
            else {
                %item_sequences = %$api->retrieve_protein_feature_sequence($feature_ids);
            }
        }
    my $comment = "Number of sequences retrieved = " . scalar keys %item_sequences;
    push @{$step_details->{comments}}, $comment;
	print STDERR "$comment\n";
    for my $seqid (keys %item_sequences) {
        my $orig_seqid = $seqid;
        if (exists $item_sequences{$seqid}) {
            my $comment = "got non-unique sequence ID $seqid, appending suffix";
            push @{$step_details->{comments}}, $comment;
            my $suffix = 2;
            while (exists $item_sequences{$seqid . "_$suffix"}) {
                $suffix++
            }
            $seqid = $seqid . "_$suffix"
        }
        $all_sequences{$seqid} = $itme_sequences{$orig_seqid}
    }
    end_step("retrieve_sequence_data");
    return($all_sequences, $aligned);
}

sub build_tree {
    my ($app, $app_def, $raw_params, $params) = @_;

    my $step_details = start_step("build_tree");
    print "Proc GeneTree build_tree ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my $time1 = `date`;
    my @outputs; # array of tuples of (filename, filetype)

    #$params = localize_params($tmpdir, $params);
    #print "after localize_params:\n", Dumper($params);
    #
    if ($debug) {
        # Write job description.
        my $json = JSON::XS->new->pretty(1);
        my $jdesc = "$tmpdir/jobdesc.json";
        write_file($jdesc, $json->encode($params));
        push @outputs, [$jdesc, 'txt']
    }
    
    print STDERR "copy data to temp dir\n";
    #print STDERR "number of data inputs is " . scalar(@{$params->{sequences}}) . "\n";
    print STDERR "params->{sequences} = $params->{sequences}\n";
    print STDERR "Aligned = $aligned, num sequences = " . scalar(keys %$all_sequences) . "\n";
    my $alignment_file_name = "$tmpdir/$params->{output_file}.msa";
    if ($aligned)
    {
        print STDERR "Write aligned seqs to $alignment_file_name\n";
        open ALIGNED_OUT, ">$alignment_file_name";
        for my $seqid (sort keys %{$all_sequences}) {
            print ALIGNED_OUT ">$seqid\n$all_sequences->{$seqid}\n"
        }
    }
    else
    { # write unaligned seqs to file, run aligner (muscle or mafft) 
        my $unaligned_file_name = "$tmpdir/$params->{output_file}_unaligned.fasta";
        print STDERR "Write unaligned seqs to $unaligned_file_name\n";
        open(OUT, ">$unaligned_file_name");
        for my $seqid (sort keys %$all_sequences) {
            $all_sequences->{$seqid} =~ tr/-//d;
            print OUT ">$seqid\n$all_sequences->{$seqid}\n"
        }
        push @outputs, [$unaligned_file_name, "feature_" . lc($params->{alphabet}) . "_fasta"]; #file type should be unaligned nt or aa
        print STDERR "Run muscle on $unaligned_file_name, ouptut to $alignment_file_name\n";
        run_muscle($unaligned_file_name, $alignment_file_name);
        push @outputs, [$alignment_file_name, "aligned_" . lc($params->{alphabet}) . "_fasta"];
    }
    if ($params->{trim_threshold} or $params->{gap_threshold})
    {
        print STDERR "performing trimming on alignment\n";
        open my $ALIGNED, $alignment_file_name or die "could not open aligned fasta";
        my $alignment = new Sequence_Alignment($ALIGNED);
        $alignment->write_stats(*STDERR);
        if (exists $params->{trim_threshold} and $params->{trim_threshold} > 0) {
            print STDERR "performing end-trimming on alignment\n";
            $alignment->end_trim($params->{trim_threshold});
            $alignment->write_stats(*STDERR);
        }
        if (exists $params->{gap_threshold} and $params->{gap_threshold} > 0) {
            print STDERR "deleting sequences with too many gaps\n";
            $alignment->delete_gappy_seqs($params->{gap_threshold});
            $alignment->write_stats(*STDERR);
        }
        my $trimmed_alignment_file_name = $alignment_file_name;
        $trimmed_alignment_file_name =~ s/\.msa/_trimmed.msa/;
        $alignment->write_fasta($trimmed_alignment_file_name);
        push @outputs, [$trimmed_alignment_file_name, "aligned_" . lc($params->{alphabet}) . "_fasta"];
        $alignment_file_name = $trimmed_alignment_file_name;
    }
    run("echo $tmpdir && ls -ltr $tmpdir");

    my $model = "AUTO"; # default for protein
    if ($params->{substitution_model}) {
        $model = $params->{substitution_model}
    }
    elsif ($params->{alphabet} =~ /DNA/i) {
        $model = "GTR"
    }
    my $output_name = $params->{output_file};
    my $alphabet = $params->{alphabet};
    my $recipe = "raxml"; #default
    if (defined $params->{recipe} and $params->{recipe}) {
        $recipe = lc($params->{recipe})
    }
    print STDERR "About to call tree program $recipe\n";
    if ($recipe eq 'raxml') {
        my @tree_outputs = run_raxml($alignment_file_name, $alphabet, $model, $output_name, $tmpdir);
        push @outputs, @tree_outputs;
    } elsif ($recipe eq 'phyml') {
        open my $ALIGNED, $alignment_file_name or die "could not open aligned fasta";
        my $alignment = new Sequence_Alignment($ALIGNED);
        $alignment_file_name =~ s/\.msa/.phy/;
        $alignment->write_phylip($alignment_file_name);
        my @tree_outputs = run_phyml($alignment_file_name, $alphabet, $model, $output_name, $tmpdir);
        push @outputs, @tree_outputs;
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    
    print STDERR '\@outputs = '. Dumper(\@outputs);
    
    my $output_folder = $app->result_folder();
    # my $output_base   = $params->{output_file};
    #$app->workspace->create( { objects => [[$path, 'folder']] } );
    
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n";
        $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
                       (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
                       $global_token);
    }
    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");
    end_step("build_tree");
}

sub run_raxml {
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;
    print STDERR "In run_raxml, alignment = $alignment_file\n";
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
    
    if ($alphabet eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = "PROTCAT". uc($model)
    }

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $parallel);
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $alignment_file);
    push @cmd, ("-n", $output_name);
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    my @outputs;
    my $bestTreeFile = $output_name . "_RAxML_bestTree.nwk";
    move("RAxML_bestTree.".$output_name, $bestTreeFile);
    push @outputs, ["$tmpdir/$bestTreeFile", 'nwk'];
    push @outputs, ["$tmpdir/RAxML_info.".$output_name, 'txt'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");

    return @outputs;
}

sub run_phyml {
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;

    my $cwd = getcwd();
    
    my $datatype = 'aa';
    if ($alphabet =~ /DNA/i) {
        $model = 'GTR';
        $datatype = 'nt'
    }

    my @cmd = ("phyml");
    push @cmd, ("-i", $alignment_file);
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    my @outputs;
    my $treeFile = $alignment_file."_phyml_tree.nwk";
    move($alignment_file."_phyml_tree.txt", $treeFile);
    my $statsFile = $alignment_file."_phyml_stats.txt";
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$statsFile, 'txt'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");

    return @outputs;
}

sub get_md5_for_feature_group {
    #" return hashref from md5 to array of patric feature id(s) - one-to-many
    my ($group, $seq_type) = @_;
    my $seq_field_md5 = (lc($seq_type) eq "protein") ? 'aa_sequence_md5' : 'na_sequence_md5';
    my $escaped = uri_escape($group);
    my $url = "$data_url/genome_feature/?in(feature_id,FeatureGroup($escaped))&select(patric_id,$seq_field_md5)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    my %md5_to_patric_id;
    for my $member (@$resp) {
        my $md5 = $member->{$seq_field_md5};
        $md5_to_patric_id{$md5} = () unless exists($md5_to_patric_id{$md5});
        push @{$md5_to_patric_id{$md5}}, $member->{'patric_id'}
    }
    if ($debug)
    {
        print "Number of patric IDs from group: ", scalar %md5_to_patric_id, "\n";
        for my $i (0..5) {
            print "member $i: $resp->[$i]{'patric_id'}, $resp->[$i]{aa_sequence_md5}\n"; 
        }
    }
    return \%md5_to_patric_id;
}

sub get_md5_for_feature_ids {
    #" return hashref from md5 to patric feature id
    my ($ids, $seq_type) = @_;
    my $seq_field_md5 = (lc($seq_type) eq "protein") ? 'aa_sequence_md5' : 'na_sequence_md5';
    my $escaped = join(",", map { uri_escape $_ } @$ids);
    my $url = "$data_url/genome_feature/?in(patric_id,($escaped))&select(patric_id,$seq_field_md5)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    my %md5_to_patric_id;
    for my $member (@$resp) {
        my $md5 = $member->{$seq_field_md5};
        $md5_to_patric_id{$md5} = () unless exists($md5_to_patric_id{$md5});
        push @{$md5_to_patric_id{$md5}}, $member->{'patric_id'}
    }
    if ($debug)
    {
        print "Number of patric IDs from group: ", scalar %md5_to_patric_id, "\n";
        for my $i (0..5) {
            print "member $i: $resp->[$i]{'patric_id'}, $resp->[$i]{aa_sequence_md5}\n"; 
        }
    }
    return \%md5_to_patric_id;
}

sub get_feature_sequences_by_md5 {
    my $md5_to_patric_id = shift;
    my $list = join(",", keys(%$md5_to_patric_id));
    my $url = "$data_url/feature_sequence/?in(md5,($list))&select(md5,sequence)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    print "response = \n\n$resp\n";
    my %patric_id_to_sequence;
    my $i = 0;
    for my $member (@$resp) {
        my $md5 = $member->{md5};
        my $sequence = $member->{'sequence'};
        for my $patric_id (@{$md5_to_patric_id->{$md5}})
        { # handle potential one-to-many relationship 
            $patric_id_to_sequence{$patric_id} = $sequence;
        }
        if ($debug and $i < 5) {
            print STDERR "$md5\t", join(",", @{$md5_to_patric_id->{$md5}}), "\t$sequence\n";
            $i++;
        }
    }
    return \%patric_id_to_sequence;
}

sub run_muscle {
    my ($unaligned, $aligned) = @_;
    my $cmd = ["muscle", "-in", $unaligned, "-out", $aligned];
    return run_cmd($cmd);
}

sub curl_text {
    my ($url) = @_;
    my @cmd = ("curl", $url, curl_options());
    my $cmd = join(" ", @cmd);
    #$cmd =~ s/sig=[a-z0-9]*/sig=XXXX/g;
    print STDERR "$cmd\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my $out = curl_text($url);
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = $global_token;
    push(@opts, "-H", "Authorization:$token");
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
    my ($tmpdir, $id) = @_;
    # return $id; # testing
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = "$tmpdir/$base";
    # return $file; # testing
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $tmpdir $base $id\n";
    system("ls -la $tmpdir");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
	die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -la $tmpdir");

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
