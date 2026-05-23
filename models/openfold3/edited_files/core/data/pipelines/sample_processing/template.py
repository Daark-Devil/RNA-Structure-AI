# Copyright 2026 AlQuraishi Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Sample processing pipelines for templates."""

from pathlib import Path
from typing import Any

import numpy as np
from biotite.structure import AtomArray
from biotite.structure.io.pdbx import CIFFile

from openfold3.core.data.primitives.quality_control.logging_utils import (
    log_runtime_memory,
)
from openfold3.core.data.primitives.structure.template import (
    TemplateSliceCollection,
    align_template_to_query,
    sample_templates,
)
from openfold3.core.data.resources.residues import MoleculeType


@log_runtime_memory(runtime_dict_key="runtime-template-proc")
def process_template_structures_of3(
    atom_array: AtomArray,
    n_templates: int,
    take_top_k: bool,
    min_n_tokens_per_chain: int,
    template_cache_directory: Path | None,
    assembly_data: dict[str, dict[str, Any]],
    template_structures_directory: Path | None,
    template_structure_array_directory: Path | None,
    template_file_format: str,
    ccd: CIFFile | None,
    use_roda_monomer_format: bool = False,
) -> TemplateSliceCollection:
    """Processes template structures for all chains of a given target structure.

    Note: During training, only looks for templates for chains that have at least one
    atom in the crop.

    Args:
        atom_array (AtomArray):
            The cropped (training) or full (inference) atom array.
        n_templates (int):
            The number of templates to sample for each chain. As per section 2.4 of the
            AF3 SI, during training at most n_templates are taken randomly from the list
            of available templates for each chain. During inference, the top (sorted by
            e-value) n_templates are taken.
        take_top_k (bool):
            Whether to take the top K templates (True) or sample randomly (False).
        min_n_tokens_per_chain (int):
            The minimum number of tokens a chain has to have for it to get template
            features.
        template_cache_directory (Path | None):
            The directory where the template cache is stored during training. For
            inference, full paths to template cache entries are provided in the
            `template_alignment_file_path` field of the `Chain` class following template
            preprocessing.
        assembly_data (dict[str, dict[str, Any]]):
            Dict containing the alignment representatives and template IDs for each
            chain.
        template_structures_directory (Path | None):
            The directory where the template structures are stored.
        template_structure_array_directory (Path | None):
            The directory where the preparsed and preprocessed template structure arrays
            are stored.
        template_file_format (str):
            The format of the template files.
        ccd (CIFFile | None):
            The parsed CCD file. Not used if template_structure_array_directory is
            provided.
        use_roda_monomer_format (bool):
            Whether template cache filepath is expected to be in the s3 RODA monomer
            format: <aln_dir>/<mgy_id>/template.npz
    Returns:
        TemplateSliceCollection:
            The sliced template atomarrays for each chain in the crop.
    """
    # Get protein chain IDs from the cropped atom array
    protein_chain_ids = np.unique(
        atom_array[
            (atom_array.molecule_type_id == MoleculeType.PROTEIN) |
            (atom_array.molecule_type_id == MoleculeType.RNA)
        ].chain_id
    )

    # ---- TPL-TRACE Stage 2: process entry ----
    print(
        f"[TPL-TRACE-PROC] process_template_structures_of3 entry: "
        f"chains_eligible={list(protein_chain_ids)} "
        f"n_templates={n_templates} take_top_k={take_top_k} "
        f"min_n_tokens_per_chain={min_n_tokens_per_chain} "
        f"template_cache_dir={template_cache_directory} "
        f"template_struct_arr_dir={template_structure_array_directory} "
        f"template_struct_dir={template_structures_directory}",
        flush=True,
    )

    if (len(protein_chain_ids) == 0) | (
        template_structure_array_directory is None
        and template_structures_directory is None
    ):
        print(
            f"[TPL-TRACE-PROC] EARLY EXIT: empty TemplateSliceCollection — "
            f"reason: chains={len(protein_chain_ids)} "
            f"struct_arr_dir={template_structure_array_directory} "
            f"struct_dir={template_structures_directory}",
            flush=True,
        )
        return TemplateSliceCollection(template_slices={})

    # Iterate over protein chains in the atom array
    # TODO: currently, this re-processes templates identical chains, add redundancy
    # logic if becomes a bottleneck
    template_slices = {}
    for chain_id in protein_chain_ids:
        atom_array_query_chain = atom_array[atom_array.chain_id == chain_id]
        _n_tokens_chain = len(np.unique(atom_array_query_chain.token_id))
        # ---- TPL-TRACE Stage 2b: per-chain entry ----
        print(
            f"[TPL-TRACE-PROC]   chain={chain_id} n_tokens={_n_tokens_chain}",
            flush=True,
        )
        # No templates if < 5 tokens in the chain (and crop during training)
        if _n_tokens_chain < min_n_tokens_per_chain:
            print(
                f"[TPL-TRACE-PROC]   chain={chain_id} SKIP — "
                f"only {_n_tokens_chain} tokens < min {min_n_tokens_per_chain}",
                flush=True,
            )
            continue

        # Sample templates and fetch their data from the cache
        sampled_template_data = sample_templates(
            assembly_data=assembly_data,
            template_cache_directory=template_cache_directory,
            n_templates=n_templates,
            take_top_k=take_top_k,
            chain_id=chain_id,
            template_structure_array_directory=template_structure_array_directory,
            template_file_format=template_file_format,
            use_roda_monomer_format=use_roda_monomer_format,
        )

        # ---- TPL-TRACE Stage 2c: sampled_template_data result ----
        print(
            f"[TPL-TRACE-PROC]   chain={chain_id} "
            f"sampled_template_data: {len(sampled_template_data)} entries, "
            f"keys={list(sampled_template_data.keys())[:8]}",
            flush=True,
        )

        # Map token positions to template atom arrays
        template_slices[chain_id] = align_template_to_query(
            sampled_template_data=sampled_template_data,
            template_structures_directory=template_structures_directory,
            template_structure_array_directory=template_structure_array_directory,
            template_file_format=template_file_format,
            ccd=ccd,
            atom_array_query_chain=atom_array_query_chain,
        )
        # ---- TPL-TRACE Stage 2d: align result ----
        print(
            f"[TPL-TRACE-PROC]   chain={chain_id} "
            f"align_template_to_query returned {len(template_slices[chain_id])} slices",
            flush=True,
        )

    print(
        f"[TPL-TRACE-PROC] EXIT: TemplateSliceCollection with "
        f"{sum(len(v) for v in template_slices.values())} total slices "
        f"across {len(template_slices)} chains",
        flush=True,
    )
    return TemplateSliceCollection(template_slices=template_slices)
