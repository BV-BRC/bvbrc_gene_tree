
# Application specification: GeneTree

This is the application specification for service with identifier GeneTree.

The backend script implementing the application is [App-GeneTree.pl](../service-scripts/App-GeneTree.pl).

The raw JSON file for this specification is [GeneTree.json](GeneTree.json).

This service performs the following task:   Estimate phylogeny of gene or other sequence feature

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| sequences | Sequence Data Inputs | ARRAY(0x55dd40a8ffb0)  | :heavy_check_mark: |  |
| alignment_program | Alignment Program | ARRAY(0x55dd40a88e28)  |  |  |
| trim_threshold | Alignment End-Trimming Threshold | float  |  |  |
| gap_threshold | Delete Gappy Sequences Threshold | float  |  |  |
| alphabet | DNA or Protein | enum  | :heavy_check_mark: |  |
| substitution_model | Substitution Model | enum  |  |  |
| bootstrap | Perform boostrapping | integer  |  |  |
| recipe | FeatureTree recipe | enum  |  | RAxML |
| tree_type | Type of tree as requested by the user | enum  |  |  |
| feature_metadata_fields | Gene Metadata Fields | string  |  |  |
| genome_metadata_fields | Genome Metadata Fields | string  |  |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

