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
                ['midpoint|m', 'Tree will be rooted in the middle of the longest path.'],
                ['quartet|q=s', 'Four* tip labels, comma-separated. Tree will be rooted below first node subtending two.'],
                ['output_name=s', 'Output filename (will have ".svg" appended if needed.'],
                ['overwrite|o', 'Overwrite existing files if any.'],
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

my $newickFile = shift;

# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}

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
#print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n" if $debug;
if ($opt->databaselink eq 'genome_id' and $opt->genomefields) {
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
    my $fields = $opt->genomefields();
    my $select = "select($fields,genome_id)";
    if ($opt->verbose) {
        print STDERR "query = $query&$select&$limit\n\n";
    }
    my ($resp, $data) = $api->submit_query('genome', "$query&$select&$limit");
    print STDERR "data retrieved = $data\n" if $opt->verbose();
    for my $record (@$data) {
        print join("||", keys %$record)."\n" if $opt->verbose();
        my $id = $record->{genome_id};
        for my $key (keys %$record) {
            if ($fields =~ /$key/) {
                $record->{$key} =~ tr/'//d; # quotes mess up embedding in javascript
                $meta_column{$key}{$id} = $record->{$key};
            }
        }
    }
    for my $field (keys %meta_column) {
        print STDERR "add alias $field to tree\n" if $opt->verbose();
        if ($opt->verbose) {
            my $limit = 4;
            for my $id (keys %{$meta_column{$field}}) {
               print "$id\t$meta_column{$field}{$id}\n";
               last unless $limit--;
            } 
        }
        $tree->add_tip_alias($field, $meta_column{$field});
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
}
if ($opt->quartet) {
    my @quartet_labels = split(",", $opt->quartet);
    $tree->root_by_quartet(@quartet_labels);
}

my $svg_file = $newickFile;
$svg_file =~ s/\.nwk$//;
$svg_file =~ s/\.tree$//;
$svg_file .= ".svg";
open F, ">$svg_file";
my $svg_data = $tree->write_svg();
print F $svg_data;
close F;

if ($workspace_dir) {
    # copy phyoxml file to user's workspace
    my @command = ("p3-cp", "-m", "svg=svg");
    push(@command, " -f") if $opt->overwrite;
    push(@command, $svg_file, "ws:" . $workspace_dir);
    print STDERR "commnd to copy back to workspace:\n@command\n" if $opt->verbose;
    my $rc = system(@command);
    die "Failure $rc running @command" unless $rc == 0;
    chdir($original_wd); # change back to the starting working directory
}
