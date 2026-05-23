{
  "queries": {
    "{{RNA}}_msa_test": {
      "use_msas": true,
      "use_main_msas": true,
      "use_paired_msas": false,
      "chains": [
        {
          "molecule_type": "rna",
          "chain_ids": "A",
          "sequence": "{{SEQ}}",
          "main_msa_file_paths": "{{WORKDIR}}/synthetic_msa/chain_A",
          "paired_msa_file_paths": null,
          "template_alignment_file_path": "{{WORKDIR}}/{{RNA}}_templates.sto",
          "template_entry_chain_ids": {{TPL_IDS}}
        }
      ]
    }
  }
}
