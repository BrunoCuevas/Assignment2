#!/usr/bin/perl
#7th March 2016
#
=head1 'network class'

=head1 'DESCRIPTION'

This class will allow us to store interaction networks for the
second assignment

=head1 'SYNOPSIS'

=cut

package network	;
use warnings	;
use strict		;
use Moose		;
use gene		;

has 'network_id' => (
	is => 'rw',
	isa => 'Str'
);
has 'go_anotation' => (
	is => 'rw',
	isa => 'ArrayRef[go]',
	required => 0
);
has 'interaction' => (
	is => 'rw',
	isa => 'HashRef[interaction]',
	required => 0
);
package interaction ;
use warnings		;
use strict			;
use Moose			;
use gene			;

has 'id' => (
	is => 'rw',
	isa => 'Str'
);
has 'gene1' => (
	is => 'rw',
	isa => 'gene'
);
has 'gene2' => (
	is => 'rw',
	isa => 'gene'
);

1;
