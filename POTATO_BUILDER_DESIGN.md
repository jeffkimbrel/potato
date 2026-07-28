# Potato Builder Skill Design

## Purpose
An AI agent skill that helps users create potato JSON files either:
1. From KEGG module definitions
2. Through conversational pathway description
3. By converting old GATOR Excel spreadsheets

## Key Requirements

### 1. Understand Potato Structure

**Gene definitions (biological entities):**
```json
{
  "id": "tktA",                    // Gene symbol for detection & reporting
  "step": [2, 5],                  // Which step(s) in pathway (array for bifunctional)
  "nodes": ["tktA_2", "tktA_5"],   // DAG node IDs it creates
  "type": "enzyme",
  "name": "transketolase",
  "ko": ["K00615"],                // KOfam terms
  "databases": {                   // Multi-database support
    "kofam118": ["K00615"],
    "my_blast_db": ["tktA_ref"],
    "pfam": ["PF00456"]
  },
  "ec": ["2.2.1.1"],
  "required": true,                // Required for pathway completion
  "is_marker": true,               // Diagnostic gene for this pathway
  "notes": "Catalyzes two reactions in PPP"
}
```

**DAG topology:**
```json
"edges": [
  {"from": "gapA_1", "to": "tktA_2", "compound": "G3P"},
  {"from": "tktA_2", "to": "foo_3"},
  {"from": "bar_4", "to": "tktA_5"},
  {"from": "tktA_5", "to": "end_6"}
]
```

### 2. Parse KEGG Module Syntax

**Input:** `(K01647,K01659) (K01681,K27802,K01682) K01637 (K01638,K19282)`

**Interpretation:**
- `(A,B)` = OR (step satisfied by A OR B)
- Space = AND (sequential steps)
- Each group = one step in pathway

**Steps:**
1. Split by spaces → steps
2. Extract KO IDs within each group
3. For OR groups, create multiple gene entries with same step number
4. Assign step numbers sequentially

### 3. Identify Marker Genes

**Heuristics for is_marker = true:**
- Pathway-specific (not in central metabolism)
- Mentioned as "diagnostic" or "key" in KEGG/literature
- Example: aceA (isocitrate lyase) is diagnostic for glyoxylate cycle
- Ask user to confirm if uncertain

### 4. Handle Bifunctional Enzymes

**Detection patterns:**
- KEGG shows enzyme at multiple reaction steps
- Same KO appears in non-adjacent steps
- EC number indicates multiple activities

**Action:**
- Set `step: [2, 5]` (array)
- Set `nodes: ["geneId_2", "geneId_5"]`
- Create edges for both occurrences

### 5. Workflow Options

#### Option A: KEGG Module Import
```
User: "Create potato from KEGG M00007"

Agent:
1. Fetch KEGG module definition
2. Parse module syntax
3. Fetch gene names/EC numbers for each KO
4. Identify bifunctional enzymes
5. Ask user about marker genes
6. Generate JSON
7. Validate and save
```

#### Option B: Conversational Build
```
User: "I want to create a potato for the glyoxylate shunt"

Agent:
1. Ask: "What are the main enzymes?"
2. For each enzyme: "What KO/EC numbers?"
3. Ask: "What are the sequential steps?"
4. Ask: "Any alternative enzymes for the same step?"
5. Ask: "Which genes are diagnostic markers?"
6. Generate JSON
7. Show preview, iterate
```

#### Option C: GATOR Spreadsheet Conversion
```
User: "Convert this GATOR Excel file to potatoes"

Agent:
1. Read Excel, identify pathway definitions
2. Parse old syntax (-> = sequential, | = OR, + = AND within step)
3. Map old tool columns (ko, blast_terms, pfam) to new databases field
4. Detect marker genes from annotations
5. Generate one JSON per pathway
6. Report conversion issues
```

### 6. Validation Checks

Before finalizing, agent should:
- ✓ All nodes have unique IDs
- ✓ All edges reference valid node IDs (with _step suffix)
- ✓ Each gene with multiple steps has nodes array
- ✓ At least one marker gene identified
- ✓ Step numbers are sequential (1, 2, 3...)
- ✓ No cycles in DAG
- ✓ At least one detection term per gene (ko, databases, etc.)

### 7. Example Interaction

```
User: Build a potato for nitrogen fixation, Mo-nitrogenase only

Agent: I'll help create a nitrogen fixation potato. The Mo-nitrogenase system has:
  - nifD (K02586) - MoFe protein alpha
  - nifK (K02591) - MoFe protein beta  
  - nifH (K02588) - Fe protein

Which should be the marker gene(s)?

User: nifH is diagnostic

Agent: Got it. Are these sequential steps or alternatives?

User: Sequential: nifD + nifK form complex, then nifH donates electrons

Agent: [Generates JSON with 3 steps, nifH marked as is_marker: true]
Here's the potato. Review?

User: Looks good!

Agent: Saved to inst/potatoes/nitrogen_fixation_Mo.json ✓
```

## Implementation as Skill

The skill should:
1. **Use WebFetch** for KEGG API calls
2. **Use Read/Write** for file I/O
3. **Use AskUserQuestion** for clarifications
4. **Validate** using `load_potato()` and `validate_potato()`
5. **Output** formatted JSON with proper indentation

## Future Enhancements

- Auto-detect common enzyme variants (e.g., MDH1, MDH2, mdh)
- Suggest sub-potatoes for disconnected subgraphs
- Integrate with literature search for marker gene identification
- Visual DAG preview before saving
