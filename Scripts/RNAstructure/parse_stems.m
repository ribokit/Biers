function stems = parse_stems( structure )

bps = convert_structure_to_bps( structure );
stems = {};
stems = parse_stems_from_bps( bps )
