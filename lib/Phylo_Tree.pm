package Phylo_Tree;
use Phylo_Node;
use strict;
use warnings;
our $debug = 0; #static variable within this class
sub set_debug { $debug = shift;
    Phylo_Node::set_debug($debug)}

sub new {
    my ($class, $newick, $support_type) = @_;
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
    $self->{_tips} = {};
    $self->{_support_type} = {$support_type} if $support_type;
    $self->{'_property_datatype'} = {};
    $self->{'_property_applies_to'} = {};
    if ($newick) {
        $self->read_newick($newick)
    }
    return $self;
}

sub get_ntips { my $self = shift; return scalar(keys %{$self->{_tips}})}
sub get_length { my $self = shift; return $self->{_length}}
sub get_support_type { my $self = shift; return defined $self->{_support_type} ? $self->{_support_type} : "support" }
sub register_tip { 
    my ($self, $name, $node) = @_; 
    print STDERR "register tip:\t$name\t$node\n" if $debug > 2;
    $self->{_tips}{$name} = $node;
}

sub list_tips {
    my $self = shift;
    my $retval = "Number of tips = " . scalar keys %{$self->{_tips}};
    $retval .= "\n";
    for my $name (keys %{$self->{_tips}}) {
        my $node = $self->{_tips}{$name};
        $retval .= "$name\t$node->{_level}\n";
    }
    $retval;
}

sub add_tip_metadata {
    my ($self, $tip_name, $field_name, $value, $provenance) = @_;
    print STDERR "in Phylo_Tree.add_tip_metadata($self, $tip_name, $field_name, $value)\n";
    unless (exists $self->{_tips}{$tip_name}) {
        print STDERR "Hey! $tip_name not found in $self->{_tips}\n";
        #warn "Tried to add metadata to non-existent tip $tip_name";
        #return;
    }
    my $node = $self->{_tips}{$tip_name};
    print STDERR "In Phylo_Tree.add_tip_metadata: about to call node->add_single_property($field_name, $value) on $node\n" if $debug > 2;
    $node->add_single_property($field_name, $value);
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

sub add_tip_properties {
    my ($self, $prop_hashref) = @_;
    print STDERR "in:add_tip_properties($self, $prop_hashref)\n";
    print STDERR "prop_hashref = %{$prop_hashref}\n";
    my ($applies_to, $datatype, $ref) = ('node', 'xsd:string', undef);
    $applies_to = $prop_hashref->{'QName:applies_to'} if $prop_hashref->{'QName:applies_to'};
    $datatype = $prop_hashref->{'QName:datatype'} if $prop_hashref->{'QName:datatype'};
    $ref = $prop_hashref->{'column_head'} if $prop_hashref->{'column_head'};
    $ref = $prop_hashref->{'QName:ref'} if $prop_hashref->{'QName:ref'};
    $self->{'_property_datatype'}{$ref} = $datatype;
    $self->{'_property_applies_to'}{$ref} = $applies_to;
    for my $tip_name (keys %{$self->{_tips}}) {
        print STDERR "try node $tip_name\n";
        if (exists $prop_hashref->{$tip_name}) {
            my $value = $prop_hashref->{$tip_name};
            my $node = $self->{_tips}{$tip_name};
            $node->add_property($ref, $value);
        }
    }
}

sub add_bulk_properties {
    my ($self, $metadata, $namespace) = @_;
    $metadata->{namespace} = "BVBRC";
    $self->{_root}->add_bulk_properties_recursive($metadata);
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
