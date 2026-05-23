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

"""
This module contains building blocks for target and ground truth structure feature
generation.
"""

import biotite.structure as struc
import numpy as np
import torch
from biotite.structure import AtomArray

from openfold3.core.data.primitives.structure.cleanup import filter_fully_atomized_bonds
from openfold3.core.data.primitives.structure.labels import get_token_starts
from openfold3.core.utils.atomize_utils import broadcast_token_feat_to_atoms


def encode_one_hot(x: torch.Tensor, num_classes: int) -> torch.Tensor:
    """Encodes a tensor of indices as a one-hot tensor.

    Args:
        x (torch.Tensor):
            Tensor of numerically encoded residues.
        num_classes (int):
            Number of classes to encode.

    Returns:
        torch.Tensor:
            One-hot encoded tensor.
    """
    x_one_hot = torch.zeros(*x.shape, num_classes, device=x.device)
    x_one_hot.scatter_(-1, x.unsqueeze(-1), 1)
    return x_one_hot


def create_sym_id(
    entity_ids_per_chain: np.ndarray, atom_array: AtomArray, token_starts: np.ndarray
) -> np.ndarray:
    """Creates sym_id feature as outlined in AF3 SI Table 5.

    Args:
        entity_ids (np.array):
            Entity ids of the target or ground truth structure.

    Returns:
        np.ndarray:
            Array of sym_ids mapped to each token.
    """
    sym_id_per_chain = np.zeros_like(entity_ids_per_chain)

    unique_entity_ids = np.unique(entity_ids_per_chain)
    entity_id_to_counter = {entity_id: 0 for entity_id in unique_entity_ids}

    for idx, i in enumerate(entity_ids_per_chain):
        entity_id_to_counter[i] += 1
        sym_id_per_chain[idx] = entity_id_to_counter[i]

    sym_id_per_atom = struc.spread_chain_wise(atom_array, sym_id_per_chain)
    return sym_id_per_atom[token_starts]


def extract_starts_entities(atom_array: AtomArray) -> tuple[np.ndarray, np.ndarray]:
    """Extracts the residue starts and entity ids from an AtomArray.

    Args:
        atom_array (AtomArray):
            AtomArray of the target or ground truth structure.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            Residue starts and entity ids for each chain.
    """
    token_starts_with_stop = get_token_starts(atom_array, add_exclusive_stop=True)
    chain_starts = struc.get_chain_starts(atom_array)
    entity_ids = atom_array.entity_id[chain_starts]
    return token_starts_with_stop, entity_ids


def create_token_bonds(atom_array: AtomArray, token_index: np.ndarray) -> torch.Tensor:
    """Creates token_bonds feature as outlined in AF3 SI Table 5.

    Args:
        atom_array (AtomArray):
            AtomArray of the target or ground truth structure.
        token_index (np.ndarray):
            Indices of tokens in the cropped AtomArray.

    Returns:
        torch.Tensor:
            token_bonds feature.
    """
    # Fully subset bonds with strict defintion to only the ones in AF3 SI Table 5
    # "token_bonds"
    atom_array = filter_fully_atomized_bonds(
        atom_array,
    )

    # Initialize N_token x N_token bond matrix
    token_bonds = np.zeros([len(token_index), len(token_index)])

    # Get bonded atoms
    bond_partners = atom_array.bonds.as_array()[:, :2]

    if bond_partners.size > 0:
        # Map atom indices to token indices to token-in-crop index
        token_to_token_in_crop = {t: tic for tic, t in enumerate(token_index)}

        get_absolute_token_id = np.vectorize(token_to_token_in_crop.get)

        bond_partner_absolute_token_ids = get_absolute_token_id(
            atom_array.token_id[bond_partners]
        )

        # Unmask corresponding bonds
        token_bonds[
            (
                bond_partner_absolute_token_ids[:, 0],
                bond_partner_absolute_token_ids[:, 1],
            ),
            (
                bond_partner_absolute_token_ids[:, 1],
                bond_partner_absolute_token_ids[:, 0],
            ),
        ] = True

    # TOKEN_BONDS_SS_PATCH v2 (start) ====================================
    # v2 features:
    #   1. Pair bonds from SS dot-bracket
    #   2. Stacking bonds (diagonal) for consecutive pairs in stems
    #   3. Junction closing-pair proximity bonds
    # Triggered by env vars SS_TOKEN_BONDS_FILE + SS_TOKEN_BONDS_CHAIN_LEN
    import os as _os
    _ss_file = _os.environ.get("SS_TOKEN_BONDS_FILE")
    _ss_chain_len = _os.environ.get("SS_TOKEN_BONDS_CHAIN_LEN")
    if _ss_file and _ss_chain_len and _os.path.exists(_ss_file):
        try:
            _chain_len = int(_ss_chain_len)
            # Parse dot-bracket
            with open(_ss_file) as _f:
                _ss_lines = [l.strip() for l in _f if l.strip()]
            if len(_ss_lines) >= 2:
                _db = _ss_lines[1]
                _bracket_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
                _all_pairs = []
                for _open_b, _close_b in _bracket_pairs.items():
                    _stack = []
                    for _i, _c in enumerate(_db):
                        if _c == _open_b:
                            _stack.append(_i)
                        elif _c == _close_b and _stack:
                            _j = _stack.pop()
                            _all_pairs.append((min(_i, _j), max(_i, _j)))
                _all_pairs.sort()
                _n_tok = len(token_index)
                
                # Helper: detect stems (consecutive base pairs)
                def _detect_stems(pairs):
                    if not pairs: return []
                    pairs_sorted = sorted(pairs)
                    stems = []
                    current = [pairs_sorted[0]]
                    for prev, cur in zip(pairs_sorted, pairs_sorted[1:]):
                        if cur[0] == prev[0] + 1 and cur[1] == prev[1] - 1:
                            current.append(cur)
                        else:
                            stems.append(current)
                            current = [cur]
                    stems.append(current)
                    return stems
                
                _stems = _detect_stems(_all_pairs)
                
                _pair_count = 0
                _stack_count = 0
                _junction_count = 0
                
                # === Feature 1: Pair bonds ===
                for _i, _j in _all_pairs:
                    if _i < _chain_len and _j < _chain_len and _i < _n_tok and _j < _n_tok:
                        token_bonds[_i, _j] = 1
                        token_bonds[_j, _i] = 1
                        _pair_count += 1
                
                # === Feature 2: Stacking bonds (diagonal within stems) ===
                for _stem in _stems:
                    if len(_stem) < 2:
                        continue
                    for _k in range(len(_stem) - 1):
                        _ia, _ja = _stem[_k]
                        _ib, _jb = _stem[_k + 1]
                        # Diagonal: (ia, jb) and (ib, ja)
                        for _x, _y in [(_ia, _jb), (_ib, _ja)]:
                            if _x < _chain_len and _y < _chain_len and _x < _n_tok and _y < _n_tok:
                                if token_bonds[_x, _y] == 0:
                                    token_bonds[_x, _y] = 1
                                    token_bonds[_y, _x] = 1
                                    _stack_count += 1
                
                # === Feature 3: Junction closing-pair proximity ===
                # If two stems are sequentially adjacent (closing pair of one
                # near opening pair of next), mark their closing pairs as proximal
                for _k in range(len(_stems) - 1):
                    _stem_a = _stems[_k]
                    _stem_b = _stems[_k + 1]
                    # Closing pair of stem_a is the LAST pair in the stem list
                    _close_a = _stem_a[-1]
                    # Opening pair of stem_b is FIRST pair
                    _open_b = _stem_b[0]
                    # Check if they're sequentially adjacent at the junction
                    # close_a = (i_a, j_a), open_b = (i_b, j_b)
                    # If j_a + 1 <= i_b or i_b - j_a is small → adjacent stems
                    _seq_gap = _open_b[0] - _close_a[0]
                    if 0 < _seq_gap < 10:
                        # Mark closing-pair endpoints as proximal to opening-pair endpoints
                        for _x, _y in [
                            (_close_a[0], _open_b[0]),
                            (_close_a[1], _open_b[1]),
                        ]:
                            if _x < _chain_len and _y < _chain_len and _x < _n_tok and _y < _n_tok:
                                if token_bonds[_x, _y] == 0 and _x != _y:
                                    token_bonds[_x, _y] = 1
                                    token_bonds[_y, _x] = 1
                                    _junction_count += 1
                
                print(f"[TOKEN_BONDS_v2] pairs={_pair_count} "
                      f"stacking={_stack_count} junctions={_junction_count} "
                      f"stems={len(_stems)} "
                      f"from {_os.path.basename(_ss_file)} (n_tok={_n_tok})")
        except Exception as _e:
            print(f"[TOKEN_BONDS_v2] WARNING: injection failed: {_e}")
    # TOKEN_BONDS_SS_PATCH v2 (end) ======================================

    return torch.tensor(token_bonds, dtype=torch.int32)


def create_atom_to_token_index(
    token_mask: torch.Tensor, num_atoms_per_token: torch.Tensor
) -> torch.Tensor:
    """
    Creates mapping from atom to its corresponding token. Note that this is the
    consecutive token index starting from 0 and not the index in the actual structure.

    Args:
        token_mask:
            [*, N_token] Token mask
        num_atoms_per_token:
            [*, N_token] Number of atoms per token

    Returns:
        atom_to_token_index:
            [*, N_atom] Mapping from atom to its token index
    """
    n_token = token_mask.shape[-1]
    batch_dims = token_mask.shape[:-1]

    # Construct token index to broadcast to atoms
    token_index = (
        torch.arange(n_token, device=token_mask.device, dtype=token_mask.dtype)
        .reshape((*((1,) * len(batch_dims)), n_token))
        .repeat((*batch_dims, 1))
    )

    atom_to_token_index = broadcast_token_feat_to_atoms(
        token_mask=token_mask,
        num_atoms_per_token=num_atoms_per_token,
        token_feat=token_index,
    ).to(dtype=torch.int32)

    return atom_to_token_index


def make_chain_pair_mask_padded(
    token_chain_id: torch.Tensor, interfaces_to_include: list[tuple[int, int]]
) -> torch.Tensor:
    """Creates a pairwise mask for chains given a list of chain tuples.
    Args:
        token_chain_id:
            tensor containing all chain ids in complex
        interfaces_to_include:
            tuples with pairwise interactions to include
    Returns:
        torch.Tensor [n_chains + 1, n_chains + 1] where:
            - each value [i, j] represents whether the corresponding chain is masked or
            not
            - a 0th row and 0th column of all zeros is added as padding
    """
    largest_chain_index = torch.max(token_chain_id)
    chain_mask = torch.zeros(
        (largest_chain_index + 1, largest_chain_index + 1), dtype=torch.int
    )

    for interface_tuple in interfaces_to_include:
        chain_mask[interface_tuple[0], interface_tuple[1]] = 1
        chain_mask[interface_tuple[1], interface_tuple[0]] = 1

    return chain_mask


def make_chain_pair_labels_padded(
    token_chain_id: torch.Tensor,
    inter_chain_types: list[str | tuple],
    type_to_chain_id_pair: dict[str | tuple, int],
):
    """Creates a chain pair-wise tensor of chain pair labels.
    Args:
        token_chain_id:
            tensor containing all chain ids in complex
        inter_chain_types:
            list of chain pair types
        type_to_chain_id_pair:
            dict mapping chain pair types to chain id pairs with given type
    Returns:
        torch.Tensor [n_chains + 1, n_chains + 1] where:
            - each value [i, j] indicates an integer label associated with the
              corresponding chain pair (enum in the order of inter_chain_types)
            - a 0th row and 0th column of all zeros is added as padding
    """
    largest_chain_index = torch.max(token_chain_id)
    chain_labels = torch.zeros(
        (largest_chain_index + 1, largest_chain_index + 1), dtype=torch.int
    )

    for idx, inter_ab_ag_type in enumerate(inter_chain_types, start=1):
        chain_id_pairs = type_to_chain_id_pair[inter_ab_ag_type]
        if len(chain_id_pairs) > 0:
            for chain_id_i, chain_id_j in chain_id_pairs:
                chain_labels[chain_id_i, chain_id_j] = idx
                chain_labels[chain_id_j, chain_id_i] = idx

    return chain_labels
