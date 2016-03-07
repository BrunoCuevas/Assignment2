package gene;
use Moose;
use LWP::Simple;
has 'KEGG_ID' => (
	is => 'rw',
	isa => 'Str'
);
has 'UNIPROT_ID' => (
	is => 'rw',
	isa => 'Str'
);
has 'SEQ' => (
	is => 'rw',
	isa => 'Str',
);
sub return_uniprot {
	my ($self) = @_;
	return $self->UNIPROT_ID;
}
sub return_keggid {
	my ($self) = @_;
	return $self->KEGG_ID;
}
sub return_fasta {
	my ($self) = @_;
	my $fasta;
		$fasta = '>'.$self->UNIPROT_ID."\n".$self->SEQ;
	return $fasta;
}
1;
