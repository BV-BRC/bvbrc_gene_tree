use strict;
use Phylo_Tree;


my ($newickFile, $quartet) = @ARGV;

my $tree = new Phylo_Tree($newickFile);
my @quartet_labels = split(",", $quartet);
$tree->root_by_quartet(@quartet_labels);
my $output = $tree->write_newick();

print "$output\n";
