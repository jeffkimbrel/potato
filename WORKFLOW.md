# POTATO v0.5.0 - Workflow Diagram

## User Workflow

```mermaid
flowchart TD
    Start([User starts new project]) --> Init[initialize_potato_sack path]
    Init --> Structure[Creates folder structure:<br/>- potatoes/<br/>- genomes/<br/>- results/<br/>- .claude/<br/>- potato_config.yaml]
    
    Structure --> Create[create_sack path]
    Create --> Sack{PotatoSack<br/>S7 Object}
    
    Sack --> Validate[validate_potato<br/>potato]
    Validate --> Valid{All potatoes<br/>valid?}
    
    Valid -->|Yes| AddGenomes[add_genomes sack, path]
    Valid -->|No| Fix[Fix potato JSON]
    Fix --> Validate
    
    AddGenomes --> Registered[Genomes registered<br/>in sack@genomes]
    
    Registered --> Save[saveRDS sack, path.rds]
    Save --> RDS[(sack.rds<br/>saved)]
    
    RDS --> LoadRDS[readRDS path.rds]
    LoadRDS --> Sack
    
    style Start fill:#e1f5e1
    style Sack fill:#ffe6e6
    style RDS fill:#e6f3ff
```

## Function Call Graph

```mermaid
flowchart LR
    subgraph User Functions
        A[initialize_potato_sack]
        B[create_sack]
        C[validate_potato]
        D[add_genomes]
    end
    
    subgraph Config
        F[load_potato_config]
        G[validate_database_configs]
        H[find_potato_sack]
    end
    
    subgraph Potato Loading
        I[load_potatoes]
        J[load_potato]
        K[load_test_potato]
    end
    
    subgraph Potato Helpers
        L[get_enzyme_nodes]
        M[get_detection_terms]
        N[get_marker_genes]
        O[build_potato_graph]
        P[print_validation]
    end
    
    subgraph Classes
        Q[Potato S7 Class]
        R[PotatoSack S7 Class]
    end
    
    A -->|creates| H
    B -->|searches with| H
    B -->|calls| F
    F -->|validates| G
    B -->|calls| I
    I -->|calls| J
    J -->|creates| Q
    B -->|creates| R
    
    C -->|validates| Q
    C -->|uses| L
    C -->|uses| O
    C -->|prints with| P
    
    D -->|uses jakomics| R
    
    L -->|extracts from| Q
    M -->|extracts from| Q
    N -->|extracts from| Q
    O -->|builds from| Q
    
    style A fill:#90EE90
    style B fill:#90EE90
    style C fill:#90EE90
    style D fill:#90EE90
    style Q fill:#FFB6C1
    style R fill:#FFB6C1
```

## Data Flow

```mermaid
flowchart TD
    subgraph Input
        A[potato_config.yaml]
        B[potato JSONs]
        C[.faa/.fasta files]
    end
    
    subgraph Processing
        D[load_potato_config]
        E[load_potato]
        F[add_genomes]
    end
    
    subgraph Objects
        G[Potato S7 Objects]
        H[PotatoSack S7 Object]
        I[jakomics FILE objects]
    end
    
    subgraph Output
        J[sack.rds]
        K[Validated potatoes]
    end
    
    A -->|reads| D
    B -->|reads| E
    C -->|registers| F
    
    D --> H
    E --> G
    F --> I
    
    G --> H
    I --> H
    
    H -->|saves to| J
    G -->|validated by| K
    
    style A fill:#f0f0f0
    style B fill:#f0f0f0
    style C fill:#f0f0f0
    style J fill:#d4edda
    style K fill:#d4edda
```

## Class Structure

```mermaid
classDiagram
    class Potato {
        +character id
        +character name
        +list nodes
        +list edges
        +character tags
        +character source
        +character notes
        +list scoring
        +character json_path
        +any graph
        +print()
        +summary()
    }
    
    class PotatoSack {
        +character sack_id
        +character sack_root
        +any config
        +list potatoes
        +list genomes
        +any results
        +any scores
        +character completed_stages
        +list provenance
        +list metadata
        +print()
        +summary()
    }
    
    class Node {
        +character id
        +integer step
        +character[] nodes
        +character type
        +character name
        +list databases
        +character[] ec
        +logical required
        +logical marker
    }
    
    class Edge {
        +character from
        +character to
        +character compound
        +character kegg_compound
    }
    
    class Config {
        +character project_name
        +list paths
        +list databases
        +character config_path
    }
    
    Potato "1" *-- "many" Node : contains
    Potato "1" *-- "many" Edge : contains
    PotatoSack "1" *-- "many" Potato : contains
    PotatoSack "1" *-- "1" Config : has
    PotatoSack "1" *-- "many" FILE : references
    
    note for Node "databases: {kofam, blast, hmm, patric}"
    note for PotatoSack "genomes are jakomics FILE objects"
```

## Validation Flow

```mermaid
flowchart TD
    Start[validate_potato potato] --> CheckID{id valid?}
    
    CheckID -->|No| Error1[Error: Missing/invalid ID]
    CheckID -->|Yes| CheckName{name valid?}
    
    CheckName -->|No| Error2[Error: Missing name]
    CheckName -->|Yes| CheckNodes{nodes exist?}
    
    CheckNodes -->|No| Error3[Error: No nodes]
    CheckNodes -->|Yes| ValidateNodes[For each node...]
    
    ValidateNodes --> CheckNodeID{node has id?}
    CheckNodeID -->|No| Error4[Error: Node missing id]
    CheckNodeID -->|Yes| CheckType{type valid?}
    
    CheckType -->|No| Warn1[Warning: Non-standard type]
    CheckType -->|Yes| CheckDB{has databases?}
    
    CheckDB -->|No| Error5[Error: Enzyme missing databases]
    CheckDB -->|Yes| ValidateDB[Validate database types]
    
    ValidateDB --> CheckDBType{type in<br/>kofam, blast,<br/>hmm, patric?}
    CheckDBType -->|No| Error6[Error: Invalid database type]
    CheckDBType -->|Yes| CheckKO{kofam<br/>format OK?}
    
    CheckKO -->|No| Warn2[Warning: Invalid KO format]
    CheckKO -->|Yes| ValidateEdges[Validate edges]
    
    ValidateEdges --> CheckEdgeRefs{edges reference<br/>valid nodes?}
    CheckEdgeRefs -->|No| Error7[Error: Invalid node reference]
    CheckEdgeRefs -->|Yes| CheckDAG{is DAG?}
    
    CheckDAG -->|No| Error8[Error: Contains cycles]
    CheckDAG -->|Yes| CheckMarkers{has marker<br/>genes?}
    
    CheckMarkers -->|No| Warn3[Warning: No marker genes]
    CheckMarkers -->|Yes| Valid[Return: valid=TRUE]
    
    Error1 --> Return[Return validation result]
    Error2 --> Return
    Error3 --> Return
    Error4 --> Return
    Error5 --> Return
    Error6 --> Return
    Error7 --> Return
    Error8 --> Return
    Warn1 --> ValidateNodes
    Warn2 --> ValidateEdges
    Warn3 --> Valid
    Valid --> Return
    
    style Valid fill:#90EE90
    style Error1 fill:#FFB6B6
    style Error2 fill:#FFB6B6
    style Error3 fill:#FFB6B6
    style Error4 fill:#FFB6B6
    style Error5 fill:#FFB6B6
    style Error6 fill:#FFB6B6
    style Error7 fill:#FFB6B6
    style Error8 fill:#FFB6B6
    style Warn1 fill:#FFE4B5
    style Warn2 fill:#FFE4B5
    style Warn3 fill:#FFE4B5
```
