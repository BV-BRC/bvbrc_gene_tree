#!/usr/bin/env perl
use Phylo_Tree; # should be in lib directory
use strict;

my $tree_file = shift;
my $metadata_file = shift;

our $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
if ($debug) {
    print STDERR "debug = $debug\n" if $debug;
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}

Phylo_Tree::set_debug($debug);
my $tree = new Phylo_Tree($tree_file);

print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n";

if ($metadata_file) {
    print STDERR "read metadata from $metadata_file\n";
    my %meta_column;
    open F, $metadata_file;
    $_ = <F>;
    chomp;
    my @header = split("\t");
    shift @header; #remove first element
    for my $column_head (@header) {
        $meta_column{$column_head}{'column_head'} = $column_head;
    }
    print STDERR "Header = " . join("^^", @header) . "\n";
    while (<F>) {
        s/^#//; # remove leading pound sign, if any
        chomp;
        my @fields = split("\t");
        my $id = shift @fields; # remove first element
        print STDERR "got metadatafields for $id\n";
        for my $i (0 .. $#header) {
            my $column_head = $header[$i];
            my $val = $fields[$i];
            $meta_column{$column_head}{$id} = $val;
        }
    }
    shift @header; # remove first element
    for my $column_head (@header) {
        print "In newick_to_phyloxml: call tree->add_tip_properties for $column_head: $meta_column{$column_head}\n";
        $tree->add_tip_properties($meta_column{$column_head});
    }
}

my $phyloxml_data = $tree->write_phyloXML();
my $phyloxml_file = $tree_file;
$phyloxml_file =~ s/\.nwk$//;
$phyloxml_file =~ s/\.tree$//;
$phyloxml_file .= ".xml";
open F, ">$phyloxml_file";
print F $phyloxml_data;
close F;
