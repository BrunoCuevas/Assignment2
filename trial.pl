use warnings;
use strict;
use LWP::Simple;
use JSON;
use gene_library;
my $gene_lib = gene_library->new(
	'number_of_networks' => 0
);
$gene_lib -> input_genes_from_kegg('kegg_ids.txt');
for (sort keys $gene_lib -> genes) {
	print $_;
	(defined $gene_lib -> genes->{$_}->PATHWAY->[0]) and
	print "\t:\t".$gene_lib -> genes->{$_}->PATHWAY->[0]->KEGG_id;
	print "\n";
}
$gene_lib -> createnetworks;
