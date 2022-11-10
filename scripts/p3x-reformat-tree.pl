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

my($opt, $usage) = P3Utils::script_opts('newickFile',
                ['annotationtsv|a=s', 'Name of a TSV file containing annotation for tips on the tree'],
                ['databaselink|link|l=s', 'Name of database field that tree identifiers map to (feature_id, genome_id, patric_id).'],
                ['genomefields|g=s', 'Comma-separated list of genome fields to annotate each tree tip.', {default => 'genome_name,family,order'}],
                ['taxon_id|t=s', 'NCBI taxon fields to annotate each tree tip.', {default => 'genome_name,family,order'}],
                ['midpoint|m', 'Tree will be rooted in the middle of the longest path.'],
                ['quartet|q=s', 'Four* tip labels, comma-separated. Tree will be rooted below first node subtending any two.'],
                ['output_name=s', 'Output filename (will have ".svg" appended if needed.'],
                ['output_format=s', 'Output format [svg|phyloxml|newick]', {default=> 'svg'}],
                ['overwrite|f', 'Overwrite existing files if any.'],
                ['verbose|debug|v', 'Write status messages to STDERR.'],
                ['name=s', 'Name for tree.'],
                ['description=s', 'Description of tree.'],
        );

# Check the parameters.

if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}

die($usage->text) if @ARGV != 1;

# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}

my $newickFile = shift;
my $original_wd = getcwd();
my $workspace_dir;
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
            print "'$workspace_newick'\n'$ls_result'\n" if $opt->verbose;
            print "Cannot access $workspace_newick\nPerhaps not logged in as user $user\n";
            exit(1);
        }
        if (-f $newickFile and not $opt->overwrite) {
            print "Refusing to overwrite local file $newickFile. Use --overwrite (or -f) to enable overwrite.\n";
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
unless (-f $newickFile) {
    print "Cannot find file $newickFile.";
    exit(1);
}
my $tree = new Phylo_Tree($newickFile, $opt->databaselink);

my $tree_file_base = $newickFile;
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
    print STDERR "data retrieved = $data\n" if $opt->verbose();
    for my $record (@$data) {
        print join("||", keys %$record)."\n" if $opt->verbose();
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
               print "$id\t$meta_column{$field}{$id}\n";
               last unless $limit--;
            } 
        }
        $tree->add_tip_annotation($field, $meta_column{$field});
        $tree_file_base .= "_$field";
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
    $tree_file_base .= "_midpoint";
}
if ($opt->quartet) {
    my @quartet_labels = split(",", $opt->quartet);
    $tree->root_by_quartet(@quartet_labels);
    $tree_file_base .= "_quarted_rooted";
}

my $output_file;
if ($opt->output_format eq 'svg') {
    $output_file = $tree_file_base . ".svg";
    open F, ">$output_file";
    my $svg_data = $tree->write_svg();
    print F $svg_data;
    close F;
}
elsif ($opt->output_format eq 'phyloxml') {
    $output_file = $tree_file_base . ".phyloxml";
    open F, ">$output_file";
    my $data = $tree->write_phyloXML();
    print F $data;
    close F;
}
elsif ($opt->output_format eq 'newick') {
    $output_file = $tree_file_base . ".nwk";
    if (-f $output_file) {
        print "File $output_file exists, not overwriting.\n";
        exit(0);
    }
    open F, ">$output_file";
    my $newick_data = $tree->write_newick();
    print F $newick_data;
    close F;
}
elsif ($opt->output_format eq 'json') {
    # need to verify or retrieve taxon_id, taxon_name, taxon_rank
    my ($taxon_id, $taxon_rank, $taxon_name);
    if ($opt->taxon_id) {
        $taxon_id = $opt->taxon_id
    }
    else {
        # compute best-fitting taxon given taxonomies of genomes
        # get lowest taxon (in taxon_lineages) that covers 80% of genomes in tree
        $taxon_id = estimate_taxon($tree->get_tip_names());
    }
    my $command = ["p3-all-taxonomies", "--eq", "taxon_id,$taxon_id", "-a", "taxon_rank,taxon_name"];
    my @stdout;
    run3( $command, undef, \@stdout);
    my ($id, $taxon_rank, $taxon_name) = split("\t", $stdout[1]);
    
    $output_file = $tree_file_base . ".json";
    open F, ">$output_file";
    print "now call write_json($taxon_name, $taxon_rank)\n" if $debug;
    my $data = $tree->write_json($taxon_name, $taxon_rank);
    print F $data;
    close F;

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
    my $command = ["p3-get-genome-data", "--nohead", "-a", "taxon_lineage_ids"];
    my @stdout;
    run3( $command, $genome_ids, \@stdout);
    
    my %taxon_id_count;
    my %taxon_id_name;
    my %taxon_id_level; # index of rank
    for my $row (@stdout) {
        chomp($row);
        my ($taxon_id, $id_lineage, $name_lineage) = split('\t', $row);
        my @taxon_ids = split("::", $id_lineage);
        my @taxon_names = split("::", $name_lineage);
        for my $i (0.. $#taxon_ids) {
            $taxon_id_count{$taxon_ids[$i]}++;
            $taxon_id_level{$taxon_ids[$i]} = $i;
        }
    }
    my $approximate_taxon_id;
    for my $taxon_id (sort {$taxon_id_count{$a}+$taxon_id_level{$a} <=> $taxon_id_count{$b}+$taxon_id_level{$b}} keys %taxon_id_count) {
        #favor lower leve (rank) taxa that meet criterion
        $approximate_taxon_id = $taxon_id;
        last if ($taxon_id_count{$taxon_id} > scalar(@$genome_ids) * $prop_taxa_required);
    }
    return ($approximate_taxon_id);
}

