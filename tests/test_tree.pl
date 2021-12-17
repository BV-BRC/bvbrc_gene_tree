use strict;
use Phylo_Tree;

Phylo_Node::set_debug(1);

my $newick_file = shift;

my $tree = new Phylo_Tree($newick_file);

$tree->list_tips();

my $nwk = $tree->write_newick();
print "As written from data stucture:\n$nwk\n";

print "Matches input = ", $nwk eq $tree->get_input_newick(), "\n";
