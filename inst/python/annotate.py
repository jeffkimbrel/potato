"""
Annotation engine for POTATO - runs bioinformatics tools via jakomics
"""
import os
import tempfile

# Import jakomics utilities
from jakomics import kegg


def extract_ko_list(potato_data):
    """Extract all KO IDs from a potato JSON for kofam searching"""
    ko_list = []
    for node in potato_data.get('nodes', []):
        databases = node.get('databases', {})
        if 'kofam' in databases:
            ko_list.extend(databases['kofam'])
    return list(set(ko_list))  # Deduplicate


def run_kofam_annotation(faa_path, potato_data, kofam_config, conda_env=None):
    """
    Run kofam annotation for a genome against a potato

    Args:
        faa_path: Path to protein FASTA file
        potato_data: Loaded potato JSON dictionary
        kofam_config: Dictionary with kofam configuration (profiles_dir or hal_path, ko_list)
        conda_env: Optional conda environment name

    Returns:
        Dictionary mapping gene_id -> list of KO hits
    """

    # Extract KOs from potato
    ko_list = extract_ko_list(potato_data)

    if len(ko_list) == 0:
        return {}

    # Get hal_path (accept either profiles_dir or hal_path)
    hal_path = kofam_config.get('hal_path') or kofam_config.get('profiles_dir')

    # Create temp directory for kofam
    with tempfile.TemporaryDirectory(prefix="potato_kofam_") as temp_dir:

        # Run kofam via jakomics
        hits = kegg.run_kofam(
            faa_path=faa_path,
            hal_path=hal_path,
            temp_dir=temp_dir,
            ko_list=kofam_config['ko_list'],
            cpus=1,
            t_scale=1.0,
            score_as_ratio=False,
            echo=False,
            run=True,
            conda_env=conda_env
        )

        # Parse hits
        parsed = kegg.parse_kofam_hits(hits)

        # Convert to gene -> KO mapping
        gene_ko_map = {}
        for ko, hit_list in parsed.items():
            for hit in hit_list:
                gene_id = hit.gene
                if gene_id not in gene_ko_map:
                    gene_ko_map[gene_id] = []
                gene_ko_map[gene_id].append({
                    'ko': ko,
                    'score': hit.score,
                    'evalue': hit.evalue,
                    'threshold': hit.threshold
                })

        return gene_ko_map


def annotate_genome_kofam(faa_path, potatoes, kofam_config, conda_env=None):
    """
    Annotate a single genome against multiple potatoes using kofam

    Args:
        faa_path: Path to protein FASTA file
        potatoes: Dictionary of potato_id -> potato_data
        kofam_config: Dictionary with kofam configuration
        conda_env: Optional conda environment name

    Returns:
        Dictionary with annotation results
    """

    results = {
        'genome': os.path.basename(faa_path),
        'potatoes': {}
    }

    # Collect all unique KOs across all potatoes
    all_kos = set()
    for potato_data in potatoes.values():
        all_kos.update(extract_ko_list(potato_data))

    if len(all_kos) == 0:
        return results

    # Get hal_path (accept either profiles_dir or hal_path)
    hal_path = kofam_config.get('hal_path') or kofam_config.get('profiles_dir')

    # Run kofam once for all KOs
    with tempfile.TemporaryDirectory(prefix="potato_kofam_") as temp_dir:
        hits = kegg.run_kofam(
            faa_path=faa_path,
            hal_path=hal_path,
            temp_dir=temp_dir,
            ko_list=kofam_config['ko_list'],
            cpus=1,
            t_scale=1.0,
            score_as_ratio=False,
            echo=False,
            run=True,
            conda_env=conda_env
        )

        parsed = kegg.parse_kofam_hits(hits)

        # Build gene -> KO mapping
        gene_ko_map = {}
        for ko, hit_list in parsed.items():
            for hit in hit_list:
                gene_id = hit.gene
                if gene_id not in gene_ko_map:
                    gene_ko_map[gene_id] = {}
                gene_ko_map[gene_id][ko] = {
                    'score': hit.score,
                    'evalue': hit.evalue,
                    'threshold': hit.threshold
                }

        # Match against each potato
        for potato_id, potato_data in potatoes.items():
            potato_kos = set(extract_ko_list(potato_data))

            # Find which potato genes were detected
            detected_nodes = []
            for node in potato_data.get('nodes', []):
                node_kos = set(node.get('databases', {}).get('kofam', []))

                # Check if any gene in the genome has this node's KOs
                for gene_id, ko_hits in gene_ko_map.items():
                    if node_kos.intersection(ko_hits.keys()):
                        detected_nodes.append({
                            'node_id': node['id'],
                            'step': node['step'],
                            'gene_id': gene_id,
                            'kos': list(node_kos.intersection(ko_hits.keys())),
                            'scores': {ko: ko_hits[ko] for ko in node_kos.intersection(ko_hits.keys())}
                        })

            results['potatoes'][potato_id] = {
                'detected_nodes': detected_nodes,
                'total_nodes': len(potato_data.get('nodes', []))
            }

    return results
