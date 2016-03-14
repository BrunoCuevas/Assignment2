package gene;
use Moose;
use LWP::Simple;
use go;
use KEGG_Pathway;
has 'KEGG_ID' => (
	is => 'rw',
	isa => 'Str',
	required => 0
);
has 'UNIPROT_ID' => (
	is => 'rw',
	isa => 'Str',
);
has 'SEQ' => (
	is => 'rw',
	isa => 'Str',
	required => 0
);
has 'PATHWAY' => (
	is => 'rw',
	isa => 'ArrayRef[KEGG_Pathway]',
	required => 0

);
has 'GO' => (
	is => 'rw',
	isa => 'HashRef[go]'
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

# sub add_interaction {
# 	my ($self, $gene) = @_ ;
# 	my @interacts = @{$self -> InteractsWith};
# 	push(@interacts, $gene)
# 		or die "ERROR : Not able to input new gene";
# 	$self->InteractsWith(\@interacts);
# }
1;
