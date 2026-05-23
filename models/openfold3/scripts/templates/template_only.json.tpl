{
  "queries": {
    "{{RNA}}_baseline": {
      "use_msas": false,
      "chains": [
        {
          "molecule_type": "rna",
          "chain_ids": "A",
          "sequence": "{{SEQ}}",
          "template_alignment_file_path": "{{WORKDIR}}/{{RNA}}_templates.sto",
          "template_entry_chain_ids": {{TPL_IDS}}
        }
      ]
    }
  }
}
