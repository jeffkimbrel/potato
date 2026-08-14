# POTATO Package Functions

## Checked Functions

- initialize_potato_sack 
- create_sack
	- find_potato_sack
	- load_potato_config
		- validate_database_configs
	- load_potatoes
		- load_potato_v2
			- PotatoV2 (S7 class)
- add_genomes
	- jakomics_to_genome_file
- run_kofam
	- create_kofam_hal
	- compute_potato_hash
	- get_potato_hashes
	- find_conda
	- run_kofam_cmd (embedded within run_kofam)
	- kofam_hits_to_tibble
- run_blast
	- create_blast_db
	- run_blast_cmd (embedded within run_blast)
	- blast_hits_to_tibble
- run_hmm
	- create_hmm_profile
	- hmm_hits_to_tibble
- plot_annotation_coverage
- print_provenance
- score_pathways
	- score_single_pathway_network
	- is_node_detected_network
- print_potato
- validate_potato
	- validate_multi_pathway
- get_verification_status
- plot_v2
- view_pathway_detail
	- build_graph_v2
	- plot_v2_interactive
	- get_compound_name (embedded in view_pathway_detail)
- potato_theme
- get_gene_results
- get_pathway_scores

## Planned but not in use

- align_coordinates
- update_potato_coordinates

## Removed

- annotate_sack_simple
- export_potato_dot
- score_single_pathway (replaced by score_single_pathway_network)
- print_multi_pathway_network
- print_pathway_compact
- format_gene_compact
- format_gene
- get_step_output_compounds
- print_validation
- validate_single_pathway
- check_verified_status (rolled into get_verification_status)
- build_visnetwork
- prepare_potato_for_plotting
- plot_potato_interactive2
- calculate_node_layout
- build_bipartite_graph
- build_visnetwork_with_gene_connectors
- create_step_layout
- get_detection_terms
- get_marker_genes
- get_node_status
- build_potato_graph
- print_pathway_detail
- load_test_potato


## Not Checked

- find_near_miss_pathways
- get_enzyme_nodes
- is_node_detected
- normalize_compound_name
- plot_all_pathways_heatmap
- plot_near_miss_pathways
- plot_pathway_prevalence
- plot_pathway_uniqueness
- summarize_missing_genes
