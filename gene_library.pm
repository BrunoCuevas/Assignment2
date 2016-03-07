package gene_library;
use Moose;
use strict;
use warnings;
use LWP::Simple;
use gene;
use JSON;
has 'genes' => (
	is => 'rw',
	isa => 'HashRef[gene]',
	required => 0
);

sub input_genes_from_kegg {
	if (@_) {
		my ($self, $path) = @_;
		if ($path) {
			if (open(INFILE, $path)) {
				#
				#
				#
				my $page					;
				my @kegg_file = <INFILE> 	;
				my $json = JSON->new		;
				my $ref_content				;
				my $refs					;
				my $uniprotreff				;
				my %hashgene				;
				my $seq						;
				foreach my $keggref (@kegg_file) {
					print "working in $keggref\n";
					$page = get ('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/dblinks.json')
						or next;
					$ref_content	=	$json->decode($page);
					$refs			=	$ref_content -> [0]	;
					$uniprotreff	=	$refs->{UniProt}->[0]	;
					$page = get ('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/aaseq.json')
						or next;
					$ref_content	=	$json->decode($page);
					$seq			=	$ref_content->[0];

					$hashgene{$keggref} = gene-> new(
						'KEGG_ID' => $keggref,
						'UNIPROT_ID' => $uniprotreff,
						'SEQ' => $seq
					);
				}
				$self -> genes(\%hashgene) ;
			}
		} else {
			die "ERROR : No path\n";
		}

	} else {
		die "ERROR : No input\n";
	}
}

1;
