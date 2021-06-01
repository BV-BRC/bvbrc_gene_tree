package Phylo_Tree;
use Phylo_Node;
use strict;
use warnings;
our $debug = 0; #static variable within this class
sub set_debug { $debug = shift() ? 1 : 0}

sub new {
    my ($class, $newick) = @_;
    if ($debug) {
        print(STDERR "in Phylo_Tree constructor\n");
        print(STDERR " class = $class\n");
        print(STDERR " input = $newick\n");
        print(STDERR " args = ", ", ".join(@_), ".\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_root} = {};
    $self->{_annot} = {};
    $self->{_ids} = [];
    $self->{_tips} = {};
    if ($newick) {
        $self->read_newick($newick)
    }
    return $self;
}

sub get_ntips { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { my $self = shift; return $self->{_length}}
sub register_tip { 
    my ($self, $name, $node) = @_; 
    print STDERR "register tip:\t$name\t$node\n" if $debug;
    $self->{_tips}->{$name} = $node;
    push @{$self->{_ids}}, $name;
}

sub list_tips {
    my $self = shift;
    my $retval = "Number of tips = " . scalar @{$self->{_ids}};
    $retval .= "\n";
    for my $name (@{$self->{_ids}}) {
        my $node = $self->{_tips}->{$name};
        $retval .= "$name\t$node->{_level}\n";
    }
    $retval;
}

sub read_newick {
    my $self = shift;
    my $newick = shift; #either a newick string or a filname
    if ( -f $newick ) {
        open F, $newick or die "Cannot open $newick";
        $newick = "";
        while(<F>) {
            chomp;
            $newick .= $_;
            last if /\);$/;
        }
    }
    $self->{'_newick'} = $newick;
    $self->{'_root'} = new Phylo_Node($newick, $self, 0);
}

sub write_newick {
    my $self = shift;
    my $retval = $self->{_root}->write_newick();
    return $retval . ";"
}

sub get_input_newick {
    my $self = shift;
    $self->{_newick}
}

sub add_properties {
    my ($self, $metadata, $namespace) = @_;
    $metadata->{namespace} = "BVBRC";
    $self->{_root}->add_properties($metadata);
}

sub write_phyloXML {
    my ($self) = @_;
    my $retval = "";
    $retval .= '<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd" xmlns="http://www.phyloxml.org">
 <phylogeny rooted="true" rerootable="true">
';
    $retval .= $self->{_root}->write_phyloXML(' ');
    $retval .= " </phylogeny>\n</phyloxml>\n";
}

1
