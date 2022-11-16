#!/usr/bin/env perl
=head1 Generate SVG version of a tree 

    p3x-tree-to-svg.pl [options] newickFile

This script reads a newick tree file and writes the tree in SVG format.

=head2 Parameters

The one positional parameter is the path of the newick file.

The command-line options are as follows.

=over 4

=item annotationtsv

A tab-delimited file where the first column contains the tip IDs on the tree.
The column headings are the names of the fields.
The values are annotations for each field (column) for each tip (row).
An optional second header row with the first field 'provenance' will be used to prefix the phy:Property fields in the PhyloXML structure (othersize the value of annotationProvenance is used).

=item databaselink 

Identifies how the tree tip identifiers are to be interpreted when linking to the BVBRC database.
Options are 'patric_id', 'genome_id', 'feature_id'.

=item featurefields

Comma separated list of feature-level fields requested for annotating tips.

=item genomefields

Comma separated list of genome-level fields requested for annotating tips.

=back

=cut

use strict;
use IPC::Run3;
use P3DataAPI;
use P3Utils;
use URI::Escape;
use Cwd;
use File::Temp;
use Phylo_Tree; # should be in lib directory

#$| = 1;
# Get the command-line options.

my($opt, $usage) = P3Utils::script_opts('xxx',
                ['infile|i=s', 'Name of the newick tree file to be reformatted (if not read from STDIN).'],
                ['annotationtsv|tsv|a=s', 'Name of a TSV file containing annotation for tips on the tree'],
                ['databaselink|link|l=s', 'Name of database field that tree identifiers map to (feature_id, genome_id, patric_id).'],
                ['genomefields|g=s', 'Comma-separated list of genome fields to annotate each tree tip.', {default => 'genome_name,family,order'}],
                ['taxon_id=s', 'NCBI taxon id.'],
                ['taxon_name=s', 'NCBI taxon name.'],
                ['taxon_rank=s', 'NCBI taxon rank.'],
                ['midpoint|m', 'Tree will be rooted in the middle of the longest path.'],
                ['quartet|q=s', 'Four* tip labels, comma-separated. Tree will be rooted below first node subtending any two.'],
                ['outfile|o=s', 'Output filename (optional).'],
                ['output_format=s', 'Output format [svg|phyloxml|newick]', {default=> 'svg'}],
                ['overwrite|f', 'Overwrite existing files if any.'],
                ['verbose|debug|v', 'Write status messages to STDERR.'],
                ['name=s', 'Name for tree.'],
                ['description=s', 'Description of tree.'],
        );

# Check the parameters.

if ($opt->verbose) {
    print STDERR "opt = %$opt\n";
    for my $key (keys %$opt) {
        print STDERR "\t$key\t$opt->{$key}\n";
    }
}

#die($usage->text) if @ARGV != 1;

# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
}

my $newickString = '';
my $newickFile = $opt->infile;
my $original_wd = getcwd();
my $workspace_dir;
if ($newickFile) {
    my @fields = split("/", $newickFile);
    if (scalar @fields > 2 and $fields[1] =~ '@') {
        print STDERR "looking for $newickFile in user workspace\n" if $opt->verbose;
        # probably a user workspace path
        my $tmpdir = File::Temp->newdir( "/tmp/TreeAnnotation_XXXXX", CLEANUP => !$debug );
        system("chmod", "755", "$tmpdir");
        print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
        chdir($tmpdir); # do all work in temporary directory
        my ($user, $brc) = split('@', $fields[1]);
        if ($brc and $brc =~ /patricbrc\.org|bvbrc/) {
            my $workspace_newick = $newickFile;
            $newickFile = pop @fields; # grab the last '/'-delimited field in user workspace path
            $workspace_dir = join('/', @fields);
            print STDERR "Parsed user=$user, brc=$brc, file=$newickFile\n" if $debug;
            my $ls_result = `p3-ls '$workspace_newick'`;
            chomp $ls_result;
            if ($ls_result ne $workspace_newick) {
                print STDERR "'$workspace_newick'\n'$ls_result'\n" if $opt->verbose;
                print STDERR "Cannot access $workspace_newick\nPerhaps not logged in as user $user\n";
                exit(1);
            }
            if (-f $newickFile and not $opt->overwrite) {
                print STDERR "Refusing to overwrite local file $newickFile. Use --overwrite (or -f) to enable overwrite.\n";
                exit(1);
            }
            my @command = ('p3-cp');
            push @command, '-f' if $opt->overwrite;
            push @command, "ws:" . $workspace_newick;
            push @command, '.';
            print STDERR "command to copy tree from workspace:\n@command\n" if $opt->verbose;
            my $rc = system(@command);
            die "Failure $rc running @command" unless $rc == 0;
            unless (-f $newickFile) {
                die "Failed to copy $workspace_newick to local file system.";
            }
        }
    }
    print STDERR "read file $newickFile\n" if $debug;
    open F, $newickFile or die "cannot open $newickFile";
    while (<F>) {
        chomp;
        $newickString .= $_;
        last if /;$/; # newick format ends with a semi-colon, end at first tree if multiple (one per line)
    }
    close F;
}
else { # read tree from stdin
    while (<STDIN>) {
        chomp;
        $newickString .= $_;
        last if /;$/; # newick format ends with a semi-colon, end at first tree if multiple (one per line)
    }
}
unless ($newickString) {
    print STDERR "Cannot find input newick string.";
    exit(1);
}
my $tree = new Phylo_Tree($newickString);
if ($opt->databaselink) {
    $tree->set_type($opt->databaselink);
}
if ($opt->name) {
    $tree->set_name($opt->name);
}
if ($opt->description) {
    $tree->set_description($opt->description);
}

# tree_file_base is used to control whether to write output to file or to STDOUT: to STDOUT if it is empty
my $tree_file_base = $newickFile;
if ($opt->outfile) {
    $tree_file_base = $opt->outfile;
}
$tree_file_base =~ s/\.nwk$//;
$tree_file_base =~ s/\.tree$//;

my $database_link = $opt->databaselink;
my $genome_fields = $opt->genomefields;
if ($opt->output_format eq 'json') {
    $database_link = 'genome_id';
    $genome_fields = 'genome_name';
}

#print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n" if $debug;
if ($database_link eq 'genome_id' and $genome_fields) {
    # Get access to PATRIC.
    my $api = P3DataAPI->new();
    my $treeIds = $tree->get_tip_names();
    my $num_tips = scalar @$treeIds;
    my $limit = "limit($num_tips)";
    print STDERR "tree IDs are: ", join(", ", @$treeIds), "\n" if $debug;
    my @escaped_IDs = map { uri_escape $_ } @$treeIds;
    my $query = "in(genome_id,(" . join(",",@escaped_IDs). "))";
    print STDERR "query=$query\n" if $debug;
    my %meta_column;
    my $select = "select($genome_fields,genome_id)";
    if ($opt->verbose) {
        print STDERR "query = $query&$select&$limit\n\n";
    }
    my ($resp, $data) = $api->submit_query('genome', "$query&$select&$limit");
    print STDERR "resp = $resp\n" if $opt->verbose();
    print STDERR "data retrieved = $data\n" if $opt->verbose();
    for my $record (@$data) {
        print STDERR join("||", keys %$record)."\n" if $opt->verbose();
        my $id = $record->{genome_id};
        for my $key (keys %$record) {
            if ($genome_fields =~ /$key/) {
                $record->{$key} =~ tr/'//d; # quotes mess up embedding in javascript
                $meta_column{$key}{$id} = $record->{$key};
            }
        }
    }
    for my $field (keys %meta_column) {
        print STDERR "add annotation $field to tree\n" if $opt->verbose();
        if ($opt->verbose) {
            my $limit = 4;
            for my $id (keys %{$meta_column{$field}}) {
               print STDERR "$id\t$meta_column{$field}{$id}\n";
               last unless $limit--;
            } 
        }
        $tree->add_tip_annotation($field, $meta_column{$field});
        $tree_file_base .= "_$field" if $tree_file_base;
    }
}
if ($opt->annotationtsv) {
    my $metadata_file = $opt->annotationtsv;
    print STDERR "reading metadata from $metadata_file\n" if $debug;
    my %meta_column;
    open F, $metadata_file;
    $_ = <F>;
    chomp;
    my @header = split("\t");
    shift @header; #remove first element
    for my $column_head (@header) {
        $meta_column{$column_head}{'column_head'} = $column_head;
    }
    print STDERR "Header fields: " . join(", ", @header) . "\n" if $debug;
    while (<F>) {
        s/^#//; # remove leading pound sign, if any
        chomp;
        my @fields = split("\t");
        my $id = shift @fields; # remove first element
        print STDERR "got metadatafields for $id\n" if $debug;
        for my $i (0 .. $#header) {
            my $column_head = $header[$i];
            my $val = $fields[$i];
            $meta_column{$column_head}{$id} = $val;
            print STDERR "\t$column_head=$val" if $debug;
        }
        print STDERR "\n" if $debug;
    }
    for my $field (keys %meta_column) {
        print STDERR "add annotation $field to tree\n" if $opt->verbose();
        if ($opt->verbose) {
            my $limit = 4;
            for my $id (keys %{$meta_column{$field}}) {
               print "$id\t$meta_column{$field}{$id}\n";
               last unless $limit--;
            } 
        }
        my $val = $meta_column{$field};
        $val = '' unless $meta_column{$field};
        $tree->add_tip_annotation($field, $val);
        $tree_file_base .= "_$field" if $tree_file_base;
    }
}

if ($opt->name) {
    $tree->set_name($opt->name);
}
if ($opt->description) {
    $tree->set_description($opt->description);
}
if ($opt->midpoint) {
    $tree->midpoint_root();
    $tree_file_base .= "_midpoint" if $tree_file_base;
}
if ($opt->quartet) {
    my @quartet_labels = split(",", $opt->quartet);
    $tree->root_by_quartet(@quartet_labels);
    $tree_file_base .= "_quarted_rooted" if $tree_file_base;
}

my $output_file;
if ($opt->output_format eq 'svg') {
    my $svg_data = $tree->write_svg();
    if ($tree_file_base) { # if data was from a file, write it to a file
        $output_file = $tree_file_base . ".svg";
        open F, ">$output_file";
        print F $svg_data;
        close F;
    }
    else {
        print STDOUT $svg_data;
    }
}
elsif ($opt->output_format eq 'phyloxml') {
    my $data = $tree->write_phyloXML();
    if ($tree_file_base) { # if data was from a file, write it to a file
        $output_file = $tree_file_base . ".phyloxml";
        open F, ">$output_file";
        print STDERR F $data;
        close F;
    }
}
elsif ($opt->output_format eq 'newick') {
    my $newick_data = $tree->write_newick();
    if ($tree_file_base) { # if data was from a file, write it to a file
        $output_file = $tree_file_base . ".nwk";
        if (-f $output_file) {
            print STDERR "File $output_file exists, not overwriting.\n";
            exit(0);
        }
        open F, ">$output_file";
        print STDERR F $newick_data;
        close F;
    }
    else {
        print STDOUT $newick_data;
    }
}
elsif ($opt->output_format eq 'json') {
    # need to verify or retrieve taxon_id, taxon_name, taxon_rank
    my ($taxon_id, $taxon_rank, $taxon_name);
    if ($opt->taxon_name and $opt->taxon_rank) {
        $taxon_name = $opt->taxon_name;
        $taxon_rank = $opt->taxon_rank;
        my $command = ["p3-all-taxonomies", "--eq", "taxon_name,$taxon_name", "--eq", "taxon_rank,$taxon_rank", "-a", "taxon_id"];
        my @stdout;
        print STDERR "running command: ", join(" ", @$command), "\n" if $debug;
        run3( $command, undef, \@stdout);
        $taxon_id = $stdout[1];
        chomp($taxon_id);
        print STDERR "results: ", join("\n", @stdout), "\n" if $debug;
    }
    else {
        if ($opt->taxon_id) {
            print STDERR "Using opt->taxon_id: ", $opt->taxon_id, " \n" if $debug;
            $taxon_id = $opt->taxon_id
        }
        else {
            # compute best-fitting taxon given taxonomies of genomes
            # get lowest taxon (in taxon_lineages) that covers 80% of genomes in tree
            print STDERR "Need to estimate taxon from genome IDs\n" if $debug;
            $taxon_id = estimate_taxon($tree->get_tip_names());
        }
        my $command = ["p3-all-taxonomies", "--eq", "taxon_id,$taxon_id", "-a", "taxon_rank,taxon_name"];
        my @stdout;
        print STDERR "running command: ", join(" ", @$command), "\n" if $debug;
        run3( $command, undef, \@stdout);
        my ($id, $taxon_rank, $taxon_name) = split("\t", $stdout[1]);
        print STDERR "results: ", join("\n", @stdout), "\n" if $debug;
    }
    
    my $data = $tree->write_json($taxon_name, $taxon_rank);
    if ($tree_file_base) {
        $output_file = $tree_file_base . ".json";
        open F, ">$output_file";
        print STDERR "now call write_json($taxon_name, $taxon_rank)\n" if $debug;
        print F $data;
        close F;
    }
    else {
        print STDOUT $data;
    }
}

if ($workspace_dir) {
    # copy phyoxml file to user's workspace
    my @command = ("p3-cp", "-m", "svg=svg", "-m", "phyloxml=phyloxml", "-m", "nwk=newick");
    push(@command, " -f") if $opt->overwrite;
    push(@command, $output_file, "ws:" . $workspace_dir);
    print STDERR "commnd to copy back to workspace:\n@command\n" if $opt->verbose;
    my $rc = system(@command);
    die "Failure $rc running @command" unless $rc == 0;
    chdir($original_wd); # change back to the starting working directory
}

sub estimate_taxon {
    my $genome_ids = shift;
    my $prop_taxa_required = 0.8;
    print STDERR "estimate_taxon from these genomes: ", join(", ", @$genome_ids), "\n\n" if $debug;
    #my $command = ["p3-get-genome-data", "--nohead", "-a", "taxon_lineage_ids"];
    my $command = ["p3-all-genomes", "--in", 'genome_id,'. join(",", @$genome_ids), '-a', 'taxon_lineage_ids'];
    print STDERR "command = ", join(" ", @$command), "\n" if $debug;
    my @stdout;
    run3( $command, $genome_ids, \@stdout);
    print STDERR "num output rows = ", scalar @stdout, "\n" if $debug;
    
    my %taxon_id_count;
    my %taxon_id_level; # index of rank
    for my $row (@stdout) {
        chomp($row);
        my ($taxon_id, $id_lineage) = split('\t', $row);
        my @taxon_ids = split("::", $id_lineage);
        pop @taxon_ids;
        shift @taxon_ids;

        print STDERR "data_row: ", join("||", @taxon_ids), "\n", if $debug;
        for my $i (0.. $#taxon_ids) {
            $taxon_id_count{$taxon_ids[$i]}++;
            $taxon_id_level{$taxon_ids[$i]} = $i;
        }
    }
    my $approximate_taxon_id;
    my $required_taxon_count = scalar(@$genome_ids) * $prop_taxa_required;
    print STDERR "need $required_taxon_count members to recognize a taxon\n" if $debug;
    #for my $taxon_id (sort {($taxon_id_count{$a}+$taxon_id_level{$a}) <=> ($taxon_id_count{$b}+$taxon_id_level{$b})} keys %taxon_id_count) {
    for my $taxon_id (sort {($taxon_id_count{$a}) <=> ($taxon_id_count{$b})} keys %taxon_id_count) {
        #favor lower leve (rank) taxa that meet criterion
        print STDERR "taxon $taxon_id, count= $taxon_id_count{$taxon_id}\n" if $debug;
        $approximate_taxon_id = $taxon_id;
        if ($taxon_id_count{$taxon_id} > $required_taxon_count) {
            print STDERR "exceeded required count with taxon $taxon_id, count=$taxon_id_count{$taxon_id}\n" if $debug;
            last
        }
    }
    print STDERR "estimate_taxon returing $approximate_taxon_id\n" if $debug;
    return ($approximate_taxon_id);
}

