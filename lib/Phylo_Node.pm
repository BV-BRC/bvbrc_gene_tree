package Phylo_Node;
use strict;
use warnings;
our $debug = 0;

sub set_debug { $debug = shift}

sub new {
    my ($class, $newick, $owner, $level) = @_;
    if ($debug > 2) {
        print(STDERR "in Phylo_Node constructor");
        print(STDERR " input = ", substr($newick, 0, 7), " len=", length($newick));
        print(STDERR " level = ", $level, ".\n");
    }
    die "need to pass newick string" unless $newick;
    my $self = {};
    bless $self, $class;
    $self->{_tree} = $owner;
    $self->{_level} = $level;
    $self->parse_newick($newick);
    return $self;
}

our %quote_char = ('"' => 1, "'", => 1);

sub parse_newick {
    my ($self, $newick) = @_;
	my $pos = 0;
    my $active_quote = 0;
    my @subclades = ();
    my $subclade = "";
    my $node_name = "";
    my $branch_length = "";
    my $state = 'none'; #first subclade, then name, then branch length
    my $open = 0;
    #print STDERR "parse newick:\n$newick\n" if $debug;
	for my $char (split('', "$newick;")) { #iterate over each charachter
        $pos++;
        if ($state eq 'none') {
            $state = $char eq '(' ? 'subclade' : 'name'
        }
        if (exists $quote_char{$char}) {
            if ( $active_quote and $char eq $active_quote) { #turn it off
                $active_quote = 0
            }
            else { #turn it on
                $active_quote = $char
            }
        }
        #print STDERR "\t$self->{_level}\t$pos\t$char\t$open\t$active_quote\t$state\n" if $debug;
        if ($state eq 'subclade') {
            if ($active_quote) { #ignode open/close parens and other terminators if in a quoted region
                    $subclade .= $char;
                }
            else { 
                if ($char eq '(') {
                    $open++;
                    if ($open == 1) {
                        $self->{_children} = () unless (exists $self->{_children});
                        print STDERR "  initiate subclade\n" if $debug > 2;
                    }
                    else { # do not put opening paren on subclade
                        $subclade .= $char
                    }
                }
                elsif ($open == 1 and ($char eq ',' or $char eq ')')) {
                    #found separator between subclades or terminator 
                    #print STDERR "Got subclade:\n", $subclade, "\n" if $debug; #$subclade\n";
                    push @subclades, $subclade;
                    $subclade = "";
                }
                else { #avoid putting separator in subclade string
                    $subclade .= $char;
                }
                $open-- if $char eq ')';
                $state = 'name' if $open == 0;
            }
		}
        elsif ($state eq 'name') {
            #if ($char eq ')' and not $node_name) {
            #    continue # terminator from preceding subclade operation
            #}
            if ($active_quote) {
                $node_name .= $char
            }
            elsif ($char eq ';' or $char eq ',' or $char eq ')' or $char eq ':') {
                $state = 'terminated';
                $state = 'branch_length' if $char eq ':';
                print STDERR " name: $node_name,  term = $char\n" if $debug > 2;
                if (exists $self->{_children}) {
                    if ($node_name =~ /^[\d\.eE-]+$/) { # if it is a number
                        $self->{_support} = $node_name
                    }
                    else {
                        $self->{_name} = $node_name;
                    }
                }
                else {
                    $self->{_name} = $node_name;
                    $self->{_tree}->register_tip($node_name, $self) unless (exists $self->{_children});
                }
            }
            else {
				$node_name .= $char;
            }
		}
		elsif ($state eq 'branch_length') {
            if ($char eq ',' or $char eq ';') {
                $state = 'terminated';
                print STDERR " bl: $branch_length, term = $char\n" if $debug > 1;
            }
            else {
				$branch_length .= $char;
			}
		}
    }
    $self->{_branch_length} = $branch_length if $branch_length;
    for my $subclade (@subclades) {
        my $child = new Phylo_Node($subclade, $self->{_tree}, $self->{_level}+1);
        push @{$self->{_children}}, $child;
    }
}

sub write_newick {
	my $self = shift;
    my $retval = "";
	if (exists $self->{'_children'}) {
        $retval = "(";
		my $first = 1;
		for my $child (@{$self->{'_children'}}) {
			$retval .= "," unless $first;
			$first = 0;
			$retval .= $child->write_newick();
		}
        $retval .= ")";
	}
	if ($self->{'_name'}) {
		$retval .= $self->{'_name'};
	}
	if ($self->{'_branch_length'}) {
		$retval .= ":$self->{'_branch_length'}"
	}
    return $retval
}

sub add_phyloxml_property {
    my ($self, $ref, $val, $data_type, $applies_to) = @_;
    $data_type = "xsd:string" unless $data_type;
    $applies_to = "node" unless $applies_to;
    $self->{_properties}{$ref} = $val;
    $self->{_property_applies_to}{$ref} = $applies_to;
    $self->{_property_datatype}{$ref} = $data_type;
    print STDERR "Phylo_Node:add_phyloxml_property just added key=$ref, val=$self->{_properties}{$ref}\n" if $debug > 2;
}

sub write_phyloXML {
    my ($self, $indent) = @_;
    print STDERR "node:write_phyloXML, self keys = ", join(", ", keys %{$self}) if $debug > 2;
    my $retval = $indent . "<clade>\n";
    if (exists $self->{'_name'} and $self->{'_name'}) {
        $retval .= $indent . " <name>$self->{'_name'}</name>\n";
    }
    if (exists $self->{_branch_length}) {
        $retval .= $indent . " <branch_length>" . $self->{_branch_length} . "</branch_length>\n";
    }
    if (exists $self->{_support}) {
        $retval .= $indent . " <confidence type=\"" . $self->{_tree}->get_support_type() . "\">" . $self->{_support} . "</confidence>\n";
    }
    if (exists $self->{_properties}) {
        print STDERR "Properties found: refs = ", join(",", keys %{$self->{_properties}}), "\n" if $debug > 2;
        for my $ref (sort keys %{$self->{_properties}}) {
            my $applies_to = $self->{_property_applies_to}{$ref};
            my $datatype = $self->{_property_datatype}{$ref};
            $retval .= $indent . " <property ref=\"$ref\" datatype=\"$datatype\" applies_to=\"$applies_to\">";
            $retval .= $self->{_properties}{$ref} . "</property>\n";
        }
    }
	if (exists $self->{'_children'}) {
		for my $child (@{$self->{'_children'}}) {
			$retval .= $child->write_phyloXML($indent . " ");
		}
	}
    $retval .= $indent . "</clade>\n";
    return $retval;
}

1
