package KEGG_Pathway;
use Moose;
has 'KEGG_id' => (
	'is' => 'rw',
	'isa' => 'Str'
);
has 'KEGG_name' => (
	'is' => 'rw',
	'isa' => 'Str'
);
1;
