use warnings;
use strict;
use LWP::Simple;
use JSON;
use gene_library;
my $gene_lib = gene_library->new();
$gene_lib -> input_genes_from_kegg('kegg_ids.txt');
