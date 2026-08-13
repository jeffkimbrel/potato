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

## Removed

- annotate_sack_simple
- export_potato_dot


## Not Checked

- align_coordinates
- print_provenance
- build_bipartite_graph
- build_graph_v2
- build_potato_graph
- build_visnetwork
- build_visnetwork_with_gene_connectors
- calculate_node_layout
- check_verified_status
- create_step_layout
- find_near_miss_pathways
- format_gene
- format_gene_compact
- get_detection_terms
- get_enzyme_nodes
- get_gene_results
- get_marker_genes
- get_node_status
- get_pathway_scores
- get_step_output_compounds
- get_verification_status
- is_node_detected
- is_node_detected_network
- jakomics_to_genome_file
- load_test_potato
- normalize_compound_name
- plot_all_pathways_heatmap
- plot_near_miss_pathways
- plot_pathway_prevalence
- plot_pathway_uniqueness
- plot_potato_interactive2
- plot_v2
- plot_v2_interactive
- potato_theme
- prepare_potato_for_plotting
- print_multi_pathway_network
- print_pathway_compact
- print_pathway_detail
- print_potato
- print_validation
- score_pathways
- score_single_pathway
- score_single_pathway_network
- summarize_missing_genes
- update_potato_coordinates
- validate_multi_pathway
- validate_potato
- validate_single_pathway
- view_pathway_detail
