#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import logging

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad_parser_h5md.parsers.parser import H5MDParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return H5MDParser()


# TODO convert towards unit testing


def test_md(parser):
    archive = EntryArchive()
    parser.parse('tests/data/test_traj_openmm_5frames.h5', archive, None)

    #######################
    # Test the NEW SCHEMA #
    #######################

    ## H5MD
    sec_simulation = archive.data
    assert sec_simulation.program.name == 'OpenMM'
    assert sec_simulation.program.version == '-1.-1.-1'
    assert len(sec_simulation.x_h5md_version) == 2
    assert sec_simulation.x_h5md_version[1] == 0
    assert sec_simulation.x_h5md_author.name == 'Joseph F. Rudzinski'
    assert sec_simulation.x_h5md_author.email == 'joseph.rudzinski@physik.hu-berlin.de'
    assert sec_simulation.x_h5md_creator.name == 'h5py'
    assert sec_simulation.x_h5md_creator.version == '3.6.0'

    ## SYSTEM
    sec_systems = sec_simulation.model_system
    assert len(sec_systems) == 5
    assert np.shape(sec_systems[0].cell[0].positions) == (31583, 3)
    assert np.shape(sec_systems[0].cell[0].velocities) == (31583, 3)
    assert sec_systems[0].cell[0].n_atoms == 31583
    assert sec_systems[0].cell[0].atoms_state[100].chemical_symbol == 'H'

    assert sec_systems[2].cell[0].positions[800][1].to('angstrom').magnitude == approx(
        26.860575
    )
    assert sec_systems[2].cell[0].velocities[1200][2].to(
        'angstrom/ps'
    ).magnitude == approx(400.0)
    assert sec_systems[3].cell[0].lattice_vectors[2][2].to(
        'angstrom'
    ).magnitude == approx(68.22318)
    assert sec_systems[0].bond_list[200][0] == 198
    assert sec_systems[0].dimensionality == 3

    sec_atoms_group = sec_systems[0].model_system
    assert len(sec_atoms_group) == 4
    assert sec_atoms_group[0].branch_label == 'group_1ZNF'
    assert sec_atoms_group[0].atom_indices[159] == 159
    sec_proteins = sec_atoms_group[0].model_system
    assert len(sec_proteins) == 1
    assert sec_proteins[0].branch_label == '1ZNF'
    assert sec_proteins[0].atom_indices[400] == 400
    sec_res_group = sec_proteins[0].model_system
    assert len(sec_res_group) == 16
    assert sec_res_group[14].branch_label == 'group_SER'
    assert sec_res_group[14].atom_indices[2] == 136
    sec_res = sec_res_group[14].model_system
    assert len(sec_res) == 3
    assert sec_res[0].branch_label == 'SER'
    assert sec_res[0].atom_indices[10] == 144
    assert sec_res[0].custom_system_attributes[0].name == 'hydrophobicity'
    assert sec_res[0].custom_system_attributes[0].value == '0.13'
    assert sec_res[0].custom_system_attributes[0].unit is None

    ## OUTPUTS
    sec_outputs = sec_simulation.outputs
    assert len(sec_outputs) == 5
    assert np.shape(sec_outputs[1].total_forces[0].value) == (31583, 3)
    assert sec_outputs[1].total_forces[0].value[2100][2].to(
        'newton'
    ).magnitude == approx(500.0)
    print(sec_outputs)
    assert sec_outputs[2].temperatures[0].value.to('kelvin').magnitude == approx(300.0)
    assert len(sec_outputs[1].x_h5md_custom_outputs) == 1
    assert (
        sec_outputs[1].x_h5md_custom_outputs[0].name
        == 'custom_thermodynamic_properties'
    )
    assert sec_outputs[1].x_h5md_custom_outputs[0].value == approx(100.0)
    assert sec_outputs[1].x_h5md_custom_outputs[0].unit == 'newton / angstrom ** 2'
    assert sec_outputs[2].time.to('ps').magnitude == approx(2.0)
    # Energies
    assert sec_outputs[2].total_energies[0].value.to('kilojoule').magnitude == approx(
        6.0
    )
    assert (
        sec_outputs[2].total_energies[0].contributions[0].m_def.name
        == 'EnergyContribution'  # 'CustomEnergy'
    )
    assert sec_outputs[2].total_energies[0].contributions[0].value.to(
        'kilojoule'
    ).magnitude == approx(3.0)
    assert (
        sec_outputs[2].total_energies[0].contributions[1].m_def.name
        == 'EnergyContribution'  # 'KineticEnergy'
    )
    assert sec_outputs[2].total_energies[0].contributions[1].value.to(
        'kilojoule'
    ).magnitude == approx(2.0)
    assert (
        sec_outputs[2].total_energies[0].contributions[2].m_def.name
        == 'EnergyContribution'  # 'PotentialEnergy'
    )
    assert sec_outputs[2].total_energies[0].contributions[2].value.to(
        'kilojoule'
    ).magnitude == approx(1.0)
    # Forces
    assert sec_outputs[2].total_forces[0].value[11].to('newton').magnitude == approx(
        500.0
    )
    assert (
        sec_outputs[2].total_forces[0].contributions[0].m_def.name
        == 'ForceContribution'  # 'CustomEnergy'
    )
    assert sec_outputs[2].total_forces[0].contributions[0].value[21].to(
        'newton'
    ).magnitude == approx(4.0)
