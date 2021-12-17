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
    open F, $metadata_file;
    $_ = <F>;
    chomp;
    my @header = split("\t");
    print STDERR "Header = " . join("^^", @header) . "\n";
    while (<F>) {
        chomp;
        my @fields = split("\t");
        my $id = $fields[0];
        for my $i (1 .. $#header) {
            my $key = $header[$i];
            my $val = $fields[$i];
            print "In newick_to_phyloxml: i=$i call tree->add_tip_metadata($id, $key, $val)\n";
            if ($val) {
                $tree->add_tip_metadata($id, $key, $val);
            }
            else {
                print "skipping because value is empty\n";
            }
        }
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
