#' Score pathway presence/absence across all genomes
#'
#' Applies quality thresholds to annotation hits and calculates pathway-level
#' completion scores. For multi-pathway networks, scores each pathway independently.
#'
#' Handles OR branches (alternative genes), required vs optional genes, and
#' multi-pathway networks where genes are shared across pathways.
#'
#' Threshold priority:
#' - Kofam: Uses per-gene threshold from KEGG (can override with kofam_threshold)
#' - HMM: Uses per-profile TC (trusted cutoff) if available, otherwise hmm_evalue
#' - BLAST: Uses global blast_evalue and blast_bitscore
#'
#' @param sack PotatoSack object with annotation results
#' @param kofam_threshold Score threshold for kofam hits (NULL = use per-gene threshold)
#' @param blast_evalue E-value threshold for BLAST hits (default: 1e-10)
#' @param blast_bitscore Bitscore threshold for BLAST hits (default: 50)
#' @param hmm_evalue E-value threshold for HMM hits without TC (default: 1e-10)
#'
#' @returns Modified PotatoSack with scores in @scores. For multi-pathway networks,
#'   scores tibble includes 'pathway' and 'pathway_name' columns with one row per
#'   pathway per genome. Scoring includes both all-steps metrics (total_steps_detected,
#'   total_steps, fraction, present) and required-only metrics (essential_total_steps_detected,
#'   essential_steps, essential_fraction, essential_pathway_present).
#' @export

score_pathways <- function(sack,
                           kofam_threshold = NULL,
                           blast_evalue = 1e-10,
                           blast_bitscore = 50,
                           hmm_evalue = 1e-10) {

  # Check that we have annotation results
  if (is.null(sack@results)) {
    cli::cli_abort("No annotation results found. Run annotation tools first.")
  }

  cli::cli_alert_info("Scoring pathways across {nrow(sack@results)} genome{?s}...")

  # Get genome names
  genome_names <- sack@results$genome

  # Initialize scores list
  all_scores <- list()

  # Process each genome
  for (i in seq_along(genome_names)) {
    genome_name <- genome_names[i]

    # Collect all hits for this genome across all tools
    genome_hits <- list()

    # Kofam hits
    if ("kofam" %in% names(sack@results)) {
      kofam_data <- sack@results$kofam[[i]]
      if (!is.null(kofam_data) && nrow(kofam_data) > 0) {
        # Filter by threshold
        if (is.null(kofam_threshold)) {
          # Use kofam's own threshold
          kofam_filtered <- kofam_data[kofam_data$score >= kofam_data$threshold, ]
        } else {
          kofam_filtered <- kofam_data[kofam_data$score >= kofam_threshold, ]
        }
        genome_hits$kofam <- kofam_filtered
      }
    }

    # BLAST hits
    if ("blast" %in% names(sack@results)) {
      blast_data <- sack@results$blast[[i]]
      if (!is.null(blast_data) && nrow(blast_data) > 0) {
        # Filter by thresholds
        blast_filtered <- blast_data[
          blast_data$evalue <= blast_evalue &
          blast_data$bitscore >= blast_bitscore, ]
        genome_hits$blast <- blast_filtered
      }
    }

    # HMM hits
    if ("hmm" %in% names(sack@results)) {
      hmm_data <- sack@results$hmm[[i]]
      if (!is.null(hmm_data) && nrow(hmm_data) > 0) {
        # Filter by threshold: use TC if available, otherwise use e-value
        hmm_filtered <- hmm_data[
          (!is.na(hmm_data$tc_threshold) & hmm_data$score >= hmm_data$tc_threshold) |
          (is.na(hmm_data$tc_threshold) & hmm_data$evalue <= hmm_evalue), ]
        genome_hits$hmm <- hmm_filtered
      }
    }

    # Score each potato for this genome
    for (potato in sack@potatoes) {
      # Score each pathway independently
      for (pathway_id in names(potato@pathways)) {
        pathway <- potato@pathways[[pathway_id]]

        pathway_score <- score_single_pathway_network(
          potato_id = potato@id,
          potato_name = potato@name,
          pathway_id = pathway_id,
          pathway = pathway,
          global_nodes = potato@genes,
          genome_name = genome_name,
          genome_hits = genome_hits
        )

        all_scores[[length(all_scores) + 1]] <- pathway_score
      }
    }
  }

  # Combine into tibble
  sack@scores <- dplyr::bind_rows(all_scores)

  # Add provenance tracking
  n_pathways <- sum(sapply(sack@potatoes, function(p) {
    if (!is.null(p@pathways) && is.list(p@pathways)) {
      length(p@pathways)
    } else {
      1
    }
  }))

  sack@provenance$scoring <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    thresholds = list(
      kofam_threshold = kofam_threshold,
      blast_evalue = blast_evalue,
      blast_bitscore = blast_bitscore,
      hmm_evalue = hmm_evalue
    ),
    n_pathways = n_pathways,
    n_genomes = length(genome_names)
  )

  cli::cli_alert_success("Scored {length(sack@potatoes)} pathway{?s} across {length(genome_names)} genome{?s}")

  sack
}



#' Score a single pathway in a multi-pathway network (internal)
#' @noRd
score_single_pathway_network <- function(potato_id, potato_name, pathway_id,
                                         pathway, global_nodes, genome_name,
                                         genome_hits) {

  # V2 schema: extract genes from edges
  pathway_edges <- pathway$edges

  if (is.null(pathway_edges) || length(pathway_edges) == 0) {
    # Empty pathway
    return(list(
      genome = genome_name,
      potato = potato_id,
      potato_name = potato_name,
      pathway = pathway_id,
      pathway_name = pathway$name,
      total_steps_detected = 0,
      total_steps = 0,
      fraction = 0,
      min_fraction = 0.75,
      present = FALSE,
      essential_total_steps_detected = NA_integer_,
      essential_steps = NA_integer_,
      essential_fraction = NA_real_,
      essential_pathway_present = NA
    ))
  }

  # Extract unique gene IDs from edges and their attributes
  gene_info <- list()
  for (edge in pathway_edges) {
    from_id <- edge$from
    to_id <- edge$to

    # Check if from/to are genes (not compounds starting with C)
    if (!is.null(from_id) && !grepl("^C\\d+", from_id)) {
      if (is.null(gene_info[[from_id]])) {
        gene_info[[from_id]] <- list(
          required = FALSE,
          marker = FALSE
        )
      }
      if (!is.null(edge$required) && edge$required) {
        gene_info[[from_id]]$required <- TRUE
      }
      if (!is.null(edge$marker) && edge$marker) {
        gene_info[[from_id]]$marker <- TRUE
      }
    }

    if (!is.null(to_id) && !grepl("^C\\d+", to_id)) {
      if (is.null(gene_info[[to_id]])) {
        gene_info[[to_id]] <- list(
          required = FALSE,
          marker = FALSE
        )
      }
      if (!is.null(edge$required) && edge$required) {
        gene_info[[to_id]]$required <- TRUE
      }
      if (!is.null(edge$marker) && edge$marker) {
        gene_info[[to_id]]$marker <- TRUE
      }
    }
  }

  if (length(gene_info) == 0) {
    # No genes in pathway
    return(list(
      genome = genome_name,
      potato = potato_id,
      potato_name = potato_name,
      pathway = pathway_id,
      pathway_name = pathway$name,
      total_steps_detected = 0,
      total_steps = 0,
      fraction = 0,
      min_fraction = 0.75,
      present = FALSE,
      essential_total_steps_detected = NA_integer_,
      essential_steps = NA_integer_,
      essential_fraction = NA_real_,
      essential_pathway_present = NA
    ))
  }

  # Build merged genes: global detection methods + pathway-specific attributes
  merged_genes <- list()
  for (gene_id in names(gene_info)) {
    # Find global gene
    global_gene <- Find(function(g) g$id == gene_id, global_nodes)

    if (is.null(global_gene)) {
      cli::cli_warn("Pathway '{pathway_id}' references gene '{gene_id}' not found in global genes")
      next
    }

    # Merge: global databases + pathway-specific required/marker
    merged_gene <- global_gene
    merged_gene$required <- gene_info[[gene_id]]$required
    merged_gene$marker <- gene_info[[gene_id]]$marker

    merged_genes[[length(merged_genes) + 1]] <- merged_gene
  }

  # V2: no steps, just count genes detected
  # Check detection for each gene (handles complexes internally)
  genes_detected <- sapply(merged_genes, function(gene) {
    is_node_detected_network(gene, potato_id, genome_hits, global_nodes)
  })

  # Track which genes are detected for gap-based scoring
  detected_gene_ids <- sapply(merged_genes[genes_detected], function(g) g$id)

  # Calculate completion (all genes)
  total_steps_detected <- sum(genes_detected)
  total_steps <- length(merged_genes)
  fraction <- if (total_steps > 0) total_steps_detected / total_steps else 0

  # Determine presence based on min_fraction threshold
  min_fraction <- pathway$scoring$min_fraction
  if (is.null(min_fraction)) min_fraction <- 0.75

  present_fraction <- fraction >= min_fraction

  # Gap-based scoring (if max_gaps and input/output defined)
  max_gaps <- pathway$scoring$max_gaps
  has_gap_scoring <- !is.null(max_gaps) &&
                     !is.null(pathway$input) && length(pathway$input) > 0 &&
                     !is.null(pathway$output) && length(pathway$output) > 0

  if (has_gap_scoring) {
    gap_result <- score_gaps(
      pathway_edges = pathway_edges,
      input_compounds = pathway$input,
      output_compounds = pathway$output,
      detected_gene_ids = detected_gene_ids,
      max_gaps = max_gaps,
      global_nodes = global_nodes
    )
    present_gaps <- gap_result$present
    gap_count <- gap_result$gap_count
    best_path <- gap_result$path
  } else {
    present_gaps <- NA
    gap_count <- NA_integer_
    best_path <- NA_character_
  }

  # Overall presence: pass if EITHER fraction OR gaps method succeeds
  if (has_gap_scoring) {
    present <- present_fraction | present_gaps
  } else {
    present <- present_fraction
  }

  # Calculate completion for required genes only
  required_genes_mask <- sapply(merged_genes, function(g) g$required %||% FALSE)
  required_genes <- merged_genes[required_genes_mask]

  if (length(required_genes) > 0) {
    essential_total_steps_detected <- sum(genes_detected[required_genes_mask])
    essential_steps <- length(required_genes)
    essential_fraction <- essential_total_steps_detected / essential_steps
    essential_pathway_present <- essential_fraction >= min_fraction
  } else {
    # No required genes defined
    essential_total_steps_detected <- NA_integer_
    essential_steps <- NA_integer_
    essential_fraction <- NA_real_
    essential_pathway_present <- NA
  }

  # Return score
  list(
    genome = genome_name,
    potato = potato_id,
    potato_name = potato_name,
    pathway = pathway_id,
    pathway_name = pathway$name %||% pathway_id,
    total_steps_detected = total_steps_detected,
    total_steps = total_steps,
    fraction = fraction,
    min_fraction = min_fraction,
    present = present,
    present_fraction = present_fraction,
    present_gaps = present_gaps,
    gap_count = gap_count,
    max_gaps = max_gaps,
    essential_total_steps_detected = essential_total_steps_detected,
    essential_steps = essential_steps,
    essential_fraction = essential_fraction,
    essential_pathway_present = essential_pathway_present
  )
}


#' Check if a gene is detected in genome hits for network pathways (internal)
#' Handles both regular genes and complexes
#' @noRd

is_node_detected_network <- function(gene, potato_id, genome_hits, global_nodes = NULL) {

  # Handle complexes: check if ALL components are detected
  if (!is.null(gene$type) && gene$type == "complex") {
    if (is.null(gene$components) || is.null(global_nodes)) {
      return(FALSE)
    }

    # Check each component
    for (component_id in gene$components) {
      component_gene <- Find(function(g) g$id == component_id, global_nodes)
      if (is.null(component_gene)) {
        return(FALSE)
      }
      # Recursively check component detection
      if (!is_node_detected_network(component_gene, potato_id, genome_hits, global_nodes)) {
        return(FALSE)
      }
    }
    # All components detected
    return(TRUE)
  }

  # Regular gene: check each database type
  databases <- gene$databases

  if (is.null(databases)) return(FALSE)

  # For network pathways, hits are stored with potato_id but we need to match
  # by node_id (gene ID) since genes are shared across pathways

  # Check kofam
  if (!is.null(databases$kofam) && !is.null(genome_hits$kofam)) {
    kofam_hits <- genome_hits$kofam
    # Match by potato_id and node_id
    gene_hits <- kofam_hits[
      kofam_hits$potato == potato_id &
      kofam_hits$node_id == gene$id, ]
    if (nrow(gene_hits) > 0) return(TRUE)
  }

  # Check blast
  if (!is.null(databases$blast) && !is.null(genome_hits$blast)) {
    blast_hits <- genome_hits$blast
    gene_hits <- blast_hits[
      blast_hits$potato == potato_id &
      blast_hits$node_id == gene$id, ]
    if (nrow(gene_hits) > 0) return(TRUE)
  }

  # Check hmm
  if (!is.null(databases$hmm) && !is.null(genome_hits$hmm)) {
    hmm_hits <- genome_hits$hmm
    gene_hits <- hmm_hits[
      hmm_hits$potato == potato_id &
      hmm_hits$node_id == gene$id, ]
    if (nrow(gene_hits) > 0) return(TRUE)
  }

  return(FALSE)
}


#' Score pathway using gap-based method (internal)
#' Binary yes/no: can ANY path exist from input→output with ≤ max_gaps missing genes?
#' @noRd
score_gaps <- function(pathway_edges, input_compounds, output_compounds,
                       detected_gene_ids, max_gaps, global_nodes) {

  # Build adjacency list for DAG traversal
  adj_list <- list()

  for (edge in pathway_edges) {
    from <- edge$from
    to <- edge$to

    if (is.null(adj_list[[from]])) {
      adj_list[[from]] <- character(0)
    }
    adj_list[[from]] <- c(adj_list[[from]], to)
  }

  # Helper: check if node is a compound
  is_compound <- function(node_id) {
    grepl("^C\\d+|^[a-z_]+_(external|internal)$", node_id)
  }

  # Helper: check if gene/complex is detected
  is_gene_detected <- function(gene_id) {
    # First check if it's directly detected
    if (gene_id %in% detected_gene_ids) {
      return(TRUE)
    }

    # Check if it's a complex - need to resolve components
    gene <- Find(function(g) g$id == gene_id, global_nodes)
    if (is.null(gene)) {
      return(FALSE)
    }

    # If it's a complex, check components
    if (!is.null(gene$type) && gene$type == "complex") {
      if (is.null(gene$components)) {
        return(FALSE)
      }
      # All components must be detected
      return(all(gene$components %in% detected_gene_ids))
    }

    # Regular gene - already checked above
    return(FALSE)
  }

  # Find all paths from input to output using DFS
  best_gap_count <- Inf
  best_path <- NULL

  for (start in input_compounds) {
    for (end in output_compounds) {

      # DFS to find all paths
      find_paths_dfs <- function(current, target, visited, path, genes_in_path) {
        if (current == target) {
          # Reached target - count gaps
          missing_genes <- genes_in_path[!sapply(genes_in_path, is_gene_detected)]
          gap_count <- length(missing_genes)

          if (gap_count < best_gap_count) {
            best_gap_count <<- gap_count
            best_path <<- paste(path, collapse = " -> ")
          }
          return()
        }

        if (current %in% visited) {
          return()  # Cycle detection
        }

        visited <- c(visited, current)
        path <- c(path, current)

        # Get neighbors
        neighbors <- adj_list[[current]]
        if (is.null(neighbors)) {
          return()
        }

        for (neighbor in neighbors) {
          # If neighbor is a gene, add to genes_in_path
          new_genes <- genes_in_path
          if (!is_compound(neighbor)) {
            new_genes <- c(new_genes, neighbor)
          }

          find_paths_dfs(neighbor, target, visited, path, new_genes)
        }
      }

      # Start DFS from this input compound
      find_paths_dfs(start, end, character(0), character(0), character(0))
    }
  }

  # Return result
  if (is.finite(best_gap_count)) {
    list(
      present = best_gap_count <= max_gaps,
      gap_count = best_gap_count,
      path = best_path
    )
  } else {
    # No path found
    list(
      present = FALSE,
      gap_count = NA_integer_,
      path = NA_character_
    )
  }
}