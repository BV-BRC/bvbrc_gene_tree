package Phylo_Node;
use strict;
use warnings;
our $debug = 0;

sub set_debug { $debug = shift() ? 1 : 0}

sub new {
    my ($class, $newick, $owner, $level) = @_;
    if ($debug) {
        print(STDERR "in Phylo_Node constructor\n");
        print(STDERR " class = $class\n");
        print(STDERR " input = ", substr($newick, 0, 7), " len=", length($newick), "\n");
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
    print STDERR "parse newick:\n$newick\n" if $debug;
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
        print STDERR "\t$self->{_level}\t$pos\t$char\t$open\t$active_quote\t$state\n" if $debug;
        if ($state eq 'subclade') {
            if ($active_quote) { #ignode open/close parens and other terminators if in a quoted region
                    $subclade .= $char;
                }
            else { 
                if ($char eq '(') {
                    $open++;
                    if ($open == 1) {
                        $self->{_children} = () unless (exists $self->{_children});
                        print STDERR "  initiate subclade\n" if $debug;
                    }
                    else { # do not put opening paren on subclade
                        $subclade .= $char
                    }
                }
                elsif ($open == 1 and ($char eq ',' or $char eq ')')) {
                    #found separator between subclades or terminator 
                    print STDERR "Got subclade:\n", $subclade, "\n" if $debug; #$subclade\n";
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
                print STDERR " name: $node_name,  term = $char\n" if $debug;
                $self->{_name} = $node_name;
                $self->{_tree}->register_tip($node_name, $self) unless (exists $self->{_children});
            }
            else {
				$node_name .= $char;
            }
		}
		elsif ($state eq 'branch_length') {
            if ($char eq ',' or $char eq ';') {
                $state = 'terminated';
                print STDERR " bl: $branch_length, term = $char\n" if $debug;
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

1
