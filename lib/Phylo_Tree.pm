package Phylo_Tree;
use Phylo_Node;
use strict;
use warnings;
our $debug = 0; #static variable within this class
sub set_debug { $debug = shift;
    Phylo_Node::set_debug($debug)}

sub new {
    my ($class, $newick) = @_;
    if ($debug) {
        print(STDERR "in Phylo_Tree constructor\n");
        print(STDERR " input = $newick\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_annot} = {};
    $self->{_tips} = ();
    $self->{_interior_nodes} = ();
    if ($newick) {
        $self->read_newick($newick)
    }
    return $self;
}

sub set_type {
    my ($self, $type) = @_;
    $self->{_type} = $type;
}

sub set_support_type {
    my ($self, $type) = @_;
    $self->{_support_type} = $type;
}

sub set_description {
    my ($self, $description) = @_;
    $self->{_description} = $description;
}

sub set_name {
    my ($self, $name) = @_;
    $self->{_name} = $name;
}

sub get_ntips { my $self = shift; return scalar(@{$self->{_tips}})}
sub get_length { my $self = shift; return $self->{_length}}
sub get_support_type { my $self = shift; return defined $self->{_support_type} ? $self->{_support_type} : "support" }
sub register_tip { 
    my ($self, $node) = @_; 
    #print STDERR "register tip:\t$node\n" if $debug > 2;
    for my $existing_node (@{$self->{_tips}}) {
        if ($node == $existing_node) {
            die "registering redundant node: $node";
        }
    }
    push @{$self->{_tips}}, $node;
}

sub register_interior_node { 
    my ($self, $node) = @_; 
    #print STDERR "register node:\t$node\n" if $debug > 2;
    for my $existing_node (@{$self->{_interior_node}}) {
        if ($node == $existing_node) {
            die "registering redundant node: $node";
        }
    }
    push @{$self->{_interior_nodes}}, $node;
}

sub get_tip_names {
    my $self = shift;
    my @retval;
    for my $tip (@{$self->{_tips}}) {
        push @retval, $tip->{_name};
    }
    return \@retval;
}

sub swap_tip_names {
    my $self = shift;
    my $swap_dict = shift;
    my $old_name_field = shift;
    my $annotation_dictionary;
    for my $tip (@{$self->{_tips}}) {
        if (exists $swap_dict->{$tip->{_name}}) {
            my $old_name = $tip->{_name};
            my $new_name = $swap_dict->{$old_name};
            $tip->{_name} = $new_name;
            $annotation_dictionary->{$new_name} = $old_name;
            #$tip->add_phyloxml_property("BVBRC:$old_name_field", $old_name);
            print STDERR "swap_tip_names: $old_name -> $new_name\n" if $debug;
        }
    }
    $self->{_annotation}{$old_name_field} = $annotation_dictionary;
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
            last if /\)\s*;$/;
        }
    }
    $self->{'_newick'} = $newick;
    $self->{'_root'} = new Phylo_Node($self, 0, $newick);
}

sub write_newick {
    my $self = shift;
    my $retval = $self->{_root}->write_newick();
    return $retval . ";\n"
}

sub get_input_newick {
    my $self = shift;
    $self->{_newick}
}

sub midpoint_root {
    my $self = shift;
    my %dist_to_tip;
    my %tip_to_tip;
    my $longest_dist = 0;
    my $distant_tip;
    my $deep_node;
    # first put all tips self dist of zero on dist_to_tip to expedite later tip-to-tip calculations as a placeholder
    for my $tip (@{$self->{_tips}}) {
        $dist_to_tip{$tip}{$tip} = 0;
        print "Tip $tip->{_name}\t$tip\t$tip->{_branch_length}\n" if $debug;
    }
    for my $node (reverse @{$self->{_interior_nodes}}) {
        #print "\nFind dist to subtended tips: $node\n" if $debug;
        $node->{dist_to_tip} = ();
        my @children = @{$node->{_children}};
        #print "Num children = ", scalar @children, "\n" if $debug;
        for my $child (@children) {
            #print "child $child, len dist_to_tip: ", scalar keys %{$dist_to_tip{$child}}, "\n" if $debug;
            for my $tip_node (keys %{$dist_to_tip{$child}}) {
                $dist_to_tip{$node}{$tip_node} = $dist_to_tip{$child}{$tip_node} + $child->{_branch_length};
                #print "Got dist_to_tip: $child, $tip_node = $dist_to_tip{$node}{$tip_node}\n" if $debug;
            }
        }
        my $num_children = scalar(@{$node->{_children}});
        for my $index1 (0..$num_children-2) {
            my $child1 = $node->{_children}->[$index1];
            for my $index2 ($index1+1..$num_children-1) {
                my $child2 = $node->{_children}->[$index2];
                for my $tip1 (keys %{$dist_to_tip{$child1}}) {
                    for my $tip2 (keys %{$dist_to_tip{$child2}}) {
                        my $dist = $dist_to_tip{$node}{$tip1} + $dist_to_tip{$node}{$tip2}; 
                        my @pair = ($tip1, $tip2);
                        #print "tip_to_tip: @pair = $dist\n";
                        if ($dist > $longest_dist) {
                            $longest_dist = $dist;
                            $distant_tip = $tip1;
                            if ($dist_to_tip{$node}{$tip2} > $dist_to_tip{$node}{$tip1}){
                                $distant_tip = $tip2;
                            }
                            $deep_node = $node;
                            my @node_pair = ($tip1, $tip2);
                            my @dist_pair = ($dist_to_tip{$node}{$tip1}, $dist_to_tip{$node}{$tip2});
                            #print "longest dist: @node_pair\t@dist_pair\n" if $debug;
                        }
                    }
                }
            }
        }
    }
    print "longest dist: $longest_dist, tip=$distant_tip\n" if $debug;
    my $half_dist = $longest_dist/2.0;
    my $node = $self->{_root};
    while ($dist_to_tip{$node}{$distant_tip} > $half_dist) {
        for my $child (@{$node->{_children}}) {
            if (exists $dist_to_tip{$child}{$distant_tip}) {
                # follow trail toward distant tip
                $node = $child;
                last;
            }
        }
    }
    my $node_above_mid = $node;
    print "half way: $node_above_mid\t$dist_to_tip{$node_above_mid}{$distant_tip}\t$node_above_mid->{_branch_length}\n" if $debug;
    my $stem = $half_dist - $dist_to_tip{$node_above_mid}{$distant_tip};
    print "half dist: $half_dist, stem=$stem\n" if $debug;
    $self->root_below_node($node_above_mid, $stem);
}

sub root_by_quartet {
    my ($self, @quartet_labels) = @_;
    print STDERR "root_by_quartet: @quartet_labels\n" if $debug;
    for my $i (0..$#quartet_labels) {
        # check for uniqueness
        for my $j ($i+1..$#quartet_labels) {
            if ($quartet_labels[$i] == $quartet_labels[$j]) {
                die "repeated label $quartet_labels[$i] in quartet: @quartet_labels";
            }
        }
    }
    my @quartet_nodes;
    for my $tip (@{$self->{_tips}}) {
        for my $label (@quartet_labels) {
            if ($label == $tip->{_name}) {
                push @quartet_nodes, $tip;
            }
        }
    }
    if ($#quartet_labels != $#quartet_nodes) {
        die "Not all quartet labels found, exiting.";
    }
    # find first interior node subtending two quartet nodes
    my $quartet_subtending_node = undef;
    print STDERR "quartet nodes: @quartet_nodes\n" if $debug;
    print STDERR "traverse interiod nodes looking for quartet ancestors.\n" if $debug;
    for my $node (reverse @{$self->{_interior_nodes}}) {
        my $matches = 0;
        print STDERR " vistit $node: " if $debug;
        for my $child (@{$node->{_children}}) {
            for my $i (0..$#quartet_nodes) {
                if ($child == $quartet_nodes[$i]) {
                    $quartet_nodes[$i] = $node;
                    $matches++;
                    print " match to child $child" if $debug;
                }
            }
        }
        print "\n" if $debug;
        if ($matches > 1) {
            $quartet_subtending_node = $node;
            last;
        }
    }
    print STDERR "search for node subtending two of quartet yielded $quartet_subtending_node\n" if $debug;
    if ($quartet_subtending_node == $self->{_root}) {
        warn "root_by_quartet found the original root, this probably should not happen";
    }
    else {
        $self->root_below_node($quartet_subtending_node);
    }
}

sub root_below_node {
    my ($self, $node_above_root, $stem) = @_; # create new root on branch below node, at distance stem if given
    print STDERR "root_below_node $node_above_root, stem=$stem\n" if $debug;
    $stem = $node_above_root->{_branch_length}/2 unless $stem;
    my $stem_remainder = $node_above_root->{_branch_length} - $stem;
    if (0 and $debug) {
        print "Verify order of interior nodes:\n";
        my %visited;
        for my $tip (@{$self->{_tips}}) {
                $visited{$tip} = 1;
        }
        for my $int_node (reverse @{$self->{_interior_nodes}}) {
            $int_node->describe();
            for my $child (@{$int_node->{_children}}) {
                die "child $child not seen before" unless $visited{$child};
            }
            $visited{$int_node} = 2;
        }
        print "Done verifying.\n";
    }
    my @path_to_root;
    my $path_node = $node_above_root;
    for my $int_node (reverse @{$self->{_interior_nodes}}) {
        # take advantage of interior node list being partially ordered, deeper nodes at beginning of list
        for my $child (@{$int_node->{_children}}) {
            if ($child == $path_node) {
                push @path_to_root, $int_node;
                $path_node = $int_node;
                last;
            }
        }
    }
    print "path_to_root: ", @path_to_root, "\n" if $debug;
    if ($debug) {
        for my $node (@path_to_root) {
            $node->describe();
        }
    }
                    

    my $new_root = new Phylo_Node($self);
    print "new node = $new_root\n" if $debug;
    $node_above_root->set_branch_length($stem);
    $new_root->add_child($node_above_root);
    my $node_below_root = shift @path_to_root; #special case of first re-oriented node, different than in subsequent loop
    $node_below_root->remove_child($node_above_root);
    my $prev_bl = $node_below_root->{_branch_length};
    $node_below_root->set_branch_length($stem_remainder);
    $new_root->add_child($node_below_root);
    my $prev_node = $node_below_root;
    print "\nnode_below_root=$node_below_root, children = @{$node_below_root->{_children}}, length = ", scalar(@{$node_below_root->{_children}}), "\n" if $debug;
    print "Now traverse path_to_root in reverse: ", reverse(@path_to_root), "\n" if $debug;
    for my $node (@path_to_root) {
        print "\nreconnect node=$node, children = @{$node->{_children}}, length = ", scalar(@{$node->{_children}}), "\n" if $debug;
        print "This is the old root.\n" if $node == $self->{_root} and $debug;
        #last if $node == $self->{_root};
        $node->remove_child($prev_node);
        $prev_node->add_child($node);
        print "node=$node, children = ", @{$node->{_children}}, "\n" if $debug;
        print "node=$node, children = @{$node->{_children}}, length = ", scalar(@{$node->{_children}}), "\n" if $debug;
        my $temp = $node->{_branch_length};
        $node->set_branch_length($prev_bl);
        $prev_bl = $temp;
        $prev_node = $node;
    }
    $self->{_root} = $new_root;
    $self->register_nodes();
}

sub register_nodes {
    my $self = shift;
    $self->{_tips} = ();
    $self->{_interior_nodes} = ();
    $self->{_root}->register();
}

sub add_tip_phyloxml_properties {
    my ($self, $prop_hashref, $ref, $default_provenance) = @_;
    print STDERR "in:add_tip_phyloxml_properties($self, $prop_hashref, $ref, $default_provenance)\n" if $debug;
    print STDERR "prop_hashref = %{$prop_hashref}\n" if $debug > 1;
    print STDERR "prop_hashref keys = ", join(" ", keys %$prop_hashref), "\n\n" if $debug > 2;
    $default_provenance = "BVBRC" unless $default_provenance;
    my ($applies_to, $datatype) = ('node', 'xsd:string');

    $applies_to = $prop_hashref->{'QName:applies_to'} if $prop_hashref->{'QName:applies_to'};
    $datatype = $prop_hashref->{'QName:datatype'} if $prop_hashref->{'QName:datatype'};
    $ref = $prop_hashref->{'column_head'} if $prop_hashref->{'column_head'};
    $ref = $prop_hashref->{'QName:ref'} if $prop_hashref->{'QName:ref'};
    $ref = "$default_provenance:$ref" unless $ref =~ /(.+):(.+)/;
    for my $tip (@{$self->{_tips}}) {
        my $tip_name = $tip->{_name};
        if (exists $prop_hashref->{$tip_name}) {
            my $value = $prop_hashref->{$tip_name};
            if ($value) {
                $tip->add_phyloxml_property($ref, $value, $datatype, $applies_to);
            }
        }
    }
}

sub get_phyloxml_properties {
    my ($self, $tip_label) = @_;
    # called from a Phylo_Node, allows keeping metadata in Phylo_Tree annotation - shared storage with svg output
    # allows has to have enties titled 'applies_to' and 'provenance' and 'data_type' to enable special cases other than 'node', 'BVBRC' and 'xsd:string'
    print STDERR "get_phyloxml_properties for $tip_label\n" if $debug;
    my @retval = ();
    for my $ref (sort keys %{$self->{_annotation}}) {
        my $provenance = "BVBRC"; #default
        my $data_type = "xsd:string";
        my $applies_to = "node";
        if (exists $self->{_annotation}->{$ref}->{'provenance'}) {
            $provenance = $self->{_annotation}->{$ref}->{'provenance'}; #allow specifying provenance per ref (field)
        }
        $applies_to = $self->{_annotation}->{$ref}->{applies_to} if $self->{_annotation}->{$ref}->{applies_to};
        $data_type = $self->{_annotation}->{$ref}->{data_type} if $self->{_annotation}->{$ref}->{data_type};
        if (exists $self->{_annotation}->{$ref}->{$tip_label}) {
            my $property = "<property ref=\"$provenance:$ref\" datatype=\"$data_type\" applies_to=\"$applies_to\">$self->{_annotation}->{$ref}->{$tip_label}</property>";
            push @retval, $property;
        }
    }
    return \@retval;
}

sub write_phyloXML {
    my $self = shift;
    my $retval = "";
    $retval .= '<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd" xmlns="http://www.phyloxml.org">
 <phylogeny rooted="true" rerootable="true"';
    $retval .= " type=\"$self->{_type}\"\n" if $self->{_type};
    $retval .= " description=\"$self->{_description}\"\n" if $self->{_description};
    $retval .= " support_type=\"$self->{_support_type}\"\n" if $self->{_support_type};
    $retval .= ">\n";
    $retval .= " <name>$self->{_name}</name>\n" if $self->{_name};
    $retval .= " <description>$self->{_description}</description>\n" if $self->{_description};
    $retval .= $self->{_root}->write_phyloXML(' ');  # recursively write root and all descendants
    $retval .= " </phylogeny>\n</phyloxml>\n";
}

sub add_tip_annotation {
    # an annotation is a dictionary where keys are tip IDs to an alternative label
    my ($self, $annotation_title, $annotation_dictionary) = @_;
    $self->{_annotation}{$annotation_title} = $annotation_dictionary;
    if ($debug) {
        print "add tip annotation: $annotation_title\t$self->{_annotation}{$annotation_title}\n";
        if ($debug) {
            my $limit = 3;
            for my $id (sort keys %{$self->{_annotation}{$annotation_title}}) {
                print "$id\t$self->{_annotation}->{$annotation_title}->{$id}\n";
                last unless --$limit;
            }
        }
    }
}

sub write_annotation_js {
    my $self = shift;
    my $retval = "<script type='text/javascript'>\n";
    for my $field (sort keys %{$self->{_annotation}}) {
        $retval .= "annotation['$field'] = {";
        print "write_annotation_js for $field, $self->{_annotation}{$field}" if $debug;
        print "write_annotation_js for $field, $self->{_annotation}->{$field}" if $debug;
        print ", len = ", scalar keys %{$self->{_annotation}->{$field}}, "\n" if $debug;
        my $limit = 4;
        for my $id (sort keys %{$self->{_annotation}->{$field}}) {
            if ($self->{_annotation}->{$field}->{$id}) {
                $retval .= "'$id':'$self->{_annotation}->{$field}->{$id}', ";
                print "annotation: $id -> $self->{_annotation}->{$field}->{$id}\n" if $debug and $limit-- > 0;
            }
        }
        $retval .= "}\n";
    }
    $retval .= "</script>\n";
    return $retval;
}

sub write_svg {
    my $self = shift;
    my $width = 800;
    my $delta_y = 15;
    my $current_y = $delta_y;
    # set the y positions of the tips
    for my $tip (@{$self->{_tips}}) {
        # explicitly specify y_pos of tips in order
        $tip->{_ypos} = $current_y;
        $current_y += $delta_y;
    }
    my $height = $current_y;
    for my $node (reverse @{$self->{_interior_nodes}}) {
        $node->embed_interior_node_y();
    }
    $self->{x_scale} = 1.0;
    $self->{_root}->embed_x(); # embed first to find max x (tip farthest from root)
    my $max_x = 0;
    for my $tip (@{$self->{_tips}}) {
        $max_x = Phylo_Node::max($tip->{_xpos}, $max_x);
        printf("%s\t%.4f\t%.4f\n", $tip->{_name}, $tip->{_xpos}, $tip->{_ypos}) if $debug;
    }
    #$self->{y_scale} = $height/$current_y;
    $self->{x_scale} = 0.5*$width/$max_x;
    $self->{_root}->embed_x(); # re-embed now that scaling is set
    if ($debug) {
        print "Num interior nodes is " . scalar @{$self->{_interior_nodes}}, "\n";
        for my $node (@{$self->{_interior_nodes}}) {
            printf("%s\t%.4f\t%.4f\t%.4f\n", $node, $node->{_branch_length}, $node->{_xpos}, $node->{_ypos});
        }
        print "root ypos = $self->{_root}->{_ypos}\n";
        print "max ypos = $current_y\n";
        print "max xpos = $max_x\n";
        printf("x_scale = %.4f\n", $self->{x_scale});
    }
    my $retval = <<END;
<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' style='border: 1px solid rgb(144, 144, 144); font-size: 12px;'>
<style>
    .highlight { fill: rgb(255,0,0); }       
    .tip_name { font-family: Arial, sans-serif; font-size: 10px}
    .click_circle { cursor: pointer; opacity: 0; r: 5 }
    .branch {fill: none; stroke-width: 3; stroke: black }
</style>
<script type='text/javascript'>
    toggle_hilight=function(target) {
       target.classList.toggle('highlight')
       var affected_tips = target.querySelectorAll('.tip_name', '.subtree');
       for (const tip of affected_tips) {
           if (target.matches('.highlight')) {tip.classList.add('highlight');}
           else {tip.classList.remove('highlight');}
        }
       if (on_selection_change_callback != 1) { 
           var selected_tip_id_list = [];
           tip_elements = document.querySelectorAll('.tip_name');
           for (const tip of tip_elements) {
                if (tip.matches('.highlight')) { selected_tip_id_list.push(tip.id); }
           }
           on_selection_change_callback(selected_tip_id_list); 
       }
  }

   var on_selection_change_callback = 1;
   set_on_selection_change_callback=function(callback) {
       on_selection_change_callback = callback;
   }
   show_id_list_to_console=function(id_list) {
         console.log("Selected IDs: "+id_list.join());
   }
   set_on_selection_change_callback(show_id_list_to_console);

var annotation = {};
relabel_tips = function(dict) {
    console.log("function relabel_tips");
    for (var id in dict) {
    document.getElementById(id).innerHTML = dict[id];
    }
}
restore_ids = function() {
    tip_elements = document.querySelectorAll('.tip_name');
    for (const tip of tip_elements) {
        tip.innerHTML = tip.id;
    }
}
  
</script>
END
    if (exists $self->{_annotation})
    {
        $retval .= $self->write_annotation_js();
        my ($label_x, $label_y, $label_delta_y) = ($width - 75, 15, 15);
        $retval .= "<text x='$label_x' y='$label_y'>Tip Label</text>\n";
        $label_y += $label_delta_y;
        for my $field (keys %{$self->{_annotation}}) {
            $retval .= "<text x='$label_x' y='$label_y' cursor='pointer' onclick=\"relabel_tips(annotation[\'$field\'])\">$field</text>\n";
            $label_y += $label_delta_y;
        }
        $retval .= "<text x='$label_x' y='$label_y' cursor='pointer' onclick='restore_ids()'>genome_id</text>\n";
        $retval .= sprintf("<path fill='none' stroke-width='1.5' stroke='gray' d='M %d,%d V %d H %d V %d Z'/>\n", $label_x-4, 4, $label_y+4, $width-4,4);
    }

    $retval .= "<g transform='translate(25,0)'>";
    $retval .= $self->{_root}->write_svg(0, $self->{_root}->{_ypos}+250);
    $retval .= "</g>\n</svg>\n";
}

sub write_json {
    my ($self, $taxon_name, $taxon_rank) = @_;
    # json fromat has 3 sections: info, labels, and tree
    my $num_tips = $self->get_ntips();
    my $retval = "{\n";
    $retval .= "\t'info': {\n\t\t'count': $num_tips,\n\t\t'taxon_name': '$taxon_name',\n\t\t'taxon_rank': '$taxon_rank'\n\t},\n";

    $retval .= "\t'labels': {\n";
    my @label_rows;
    for my $id (sort keys %{$self->{_annotation}->{genome_name}}) {
        push @label_rows, "\t\t'$id': '$self->{_annotation}->{genome_name}->{$id}'";
    }
    $retval .= join(",\n", @label_rows);
    $retval .= "\n\t},\n";

    $retval .= "\t'tree': '" . $self->write_newick() . "'\n";
    $retval .= "}\n";
    return $retval;    
}
1
