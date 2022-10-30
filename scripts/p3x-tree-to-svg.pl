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
                ['provenance|p=s', 'Provenance of annotation. (default: BVBRC)', { default => 'BVBRC' }],
                ['featurefields|f=s', 'Comma-separated list of feature fields to annotate each tree tip.', {default => 'product,accession'}],
                ['genomefields|g=s', 'Comma-separated list of genome fields to annotate each tree tip.', {default => 'species,strain,geographic_group,isolation_country,host_group,host_common_name,collection_year,genus,mlst'}],
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

my $link = $opt->databaselink;
my $tree = new Phylo_Tree($newickFile, $link);
#print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n" if $debug;

if ($opt->name) {
    $tree->set_name($opt->name);
}
if ($opt->description) {
    $tree->set_description($opt->description);
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
