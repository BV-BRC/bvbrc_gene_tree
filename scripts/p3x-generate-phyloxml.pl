#!/usr/bin/env perl
=head1 Generate PhyloXML version of a tree with metadata for tip nodes

    p3x-generate-phyloxml.pl [options] newickFile

This script reads a newick tree file and writes the tree in PhyloXML format.
Optionally, it can add annotation describing the tips of the tree.
Annotation can be conveyed in a TSV file or fetched via the BVBRC data api.

=head2 Parameters

The one positional parameter is the path of the newick file.

The command-line options are as follows.

=over 4

=item annotationTsv

A tab-delimited file where the first column contains the tip IDs on the tree.
The column headings are the names of the fields.
The values are annotations for each field (column) for each tip (row).
An optional second header row with the first field 'provenance' will be used to prefix the phy:Property fields in the PhyloXML structure (othersize the value of annotationProvenance is used).

=item databaseLink 

Identifies how the tree tip identifiers are to be interpreted when linking to the BVBRC database.
Options are 'patric_id', 'genome_id', 'feature_id'.

=item featureFields

Comma separated list of feature-level fields requested for annotating tips.

=item genomeFields

Comma separated list of genome-level fields requested for annotating tips.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Phylo_Tree; # should be in lib directory
use URI::Escape;

#$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('newickFile',
                ['annotationtsv|a=s', 'Name of a TSV file containing annotation for tips on the tree'],
                ['databaselink|link|l=s', 'Name of database field that tree identifiers map to (feature_id, genome_id, patric_id).'],
                ['provenance|p=s', 'Provenance of annotation. (default: BVBRC)', { default => 'BVBRC' }],
                ['featurefields|f=s', 'Comma-separated list of feature-level fields to annotate each tree tip.'],
                ['genomefields|g=s', 'Comma-separated list of genome-level fields to annotate each tree tip.'],
                ['verbose|debug|v', 'Write status messages to STDERR'],
        );
# Check the parameters.
my $newickFile = $ARGV[0];
if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}
if (! $newickFile) {
    print "No input file specified.";
    exit(1);
}
# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}

my $tree = new Phylo_Tree($newickFile);
#print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n" if $debug;

my %meta_column; #first key is column (field name), second key is row (tip ID)
if ($opt->annotationtsv) {
    my $metadata_file = $opt->annotationtsv;
    print STDERR "readiing metadata from $metadata_file\n" if $debug;
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
}

if ($opt->databaselink) {
    my $link = $opt->databaselink;
    # Get access to PATRIC.
    my $api = P3DataAPI->new();
    my $treeIds = $tree->get_tip_ids();
    print STDERR "tree IDs are: ", join(", ", @$treeIds), "\n" if $debug;
    my @escaped_IDs = map { uri_escape $_ } @$treeIds;
    my $query = "in($link,(" . join(",",@escaped_IDs). "))";
    print STDERR "query=$query\n" if $debug;
    if ($link =~ /feature_id|patric_id/ and $opt->featurefields) { # get feature annotation
        my $featureFields = $opt->featurefields;
        $featureFields .= ",genome_id" if $opt->genomefields and $featureFields !~ /genome_id/;
        my $select = "select($featureFields,$link)"; 
        print STDERR "select=$select\n" if $debug;
        my %genome_to_links;
        my ($resp, $data) = $api->submit_query('genome_feature', "$query&$select");
        for my $record (@$data) {
            #print join("||", keys %$record), "\n" if $debug;
            my $id = $record->{$link};
            for my $key (keys %$record) {
                next if $key eq $link;
                if ($featureFields =~ /$key/) {
                    $meta_column{$key}{$id} = $record->{$key};
                }
                if ($key eq 'genome_id') {
                    push(@{$genome_to_links{$record->{genome_id}}}, $id); # accumulate all links per genome_id
                }
            }
        }
        if (exists $meta_column{'genome_id'}) {
            my $unique_genomes = join(",", sort keys %genome_to_links);
            $query = "in(genome_id,($unique_genomes))";
            my $genomeFields = $opt->genomefields();
            $select = "select($genomeFields,genome_id)";
            print STDERR "query db for genome fields:\n$query\n$select\n" if $debug;
            my ($resp, $data) = $api->submit_query('genome', "$query&$select");
            for my $record (@$data) {
                #print join("||", keys %$record), "\n" if $debug;
                my $genome_id = $record->{genome_id};
                for my $key (keys %$record) {
                    if ($genomeFields =~ /$key/) {
                        for my $tree_id (@{$genome_to_links{$genome_id}}) {
                            $meta_column{$key}{$tree_id} = $record->{$key};
                        }
                    }
                }
            }
        }
    }
    elsif ($link eq 'genome_id' and $opt->genomefields) { # get genome annotation
        my $fields = $opt->genomefields();
        my $select = "select($fields,genome_id)";
        my ($resp, $data) = $api->submit_query('genome', "$query&$select");
        for my $record (@$data) {
            #print join("||", keys %$record), "\n" if $debug;
            my $id = $record->{genome_id};
            for my $key (keys %$record) {
                if ($fields =~ /$key/) {
                    $meta_column{$key}{$id} = $record->{$key};
                }
            }
        }
    }
}

for my $column (sort keys %meta_column) {
    print "tree->add_phyloxml_tip_properties for $column\n" if $debug;
    my $provenance = "BVBRC";
    $provenance = $opt->provenance if $opt->provenance;
    $provenance = $meta_column{$column}{provenance} if exists $meta_column{$column}{provenance};
    $tree->add_tip_phyloxml_properties($meta_column{$column}, $column, $provenance);
}

my $phyloxml_file = $newickFile;
$phyloxml_file =~ s/\.nwk$//;
$phyloxml_file =~ s/\.tree$//;
$phyloxml_file .= ".xml";
open F, ">$phyloxml_file";
my $phyloxml_data = $tree->write_phyloXML();
print F $phyloxml_data;
close F;
