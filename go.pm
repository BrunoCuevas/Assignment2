package go ;
use Moose;
has 'GO_id' => (
	'is' => 'rw',
	'isa' => 'Str'
);
has 'GO_category' => (
	'is' => 'rw',
	'isa' => 'Str'
);
has 'GO_term' => (
	'is' => 'rw',
	'isa' => 'Str'
);
1;
