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
from nomad_parser_h5md.parsers.h5md_parser import H5MDParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return H5MDParser()


# TODO convert towards unit testing


def test_md(parser):
    archive = EntryArchive()
    parser.parse('tests/data/test_traj_openmm_5frames_new.h5', archive, None)

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
    assert np.shape(sec_systems[0].positions) == (31583, 3)
    assert np.shape(sec_systems[0].velocities) == (31583, 3)
    assert sec_systems[0].n_particles == 31583
    assert sec_systems[0].particle_states[100].chemical_symbol == 'H'
    assert sec_systems[0].particle_states[100].label == 'H'

    assert sec_systems[2].positions[800][1].to('angstrom').magnitude == approx(
        26.860575
    )
    assert sec_systems[2].velocities[1200][2].to('angstrom/ps').magnitude == approx(
        400.0
    )
    assert sec_systems[3].cell[0].lattice_vectors[2][2].to(
        'angstrom'
    ).magnitude == approx(68.22318)
    assert sec_systems[3].cell[0].periodic_boundary_conditions == [True, True, True]
    # assert sec_systems[0].particle_states[200].atom_indices[0] == 200
    assert sec_systems[0].bond_list[200][0] == 198
    assert sec_systems[0].dimensionality == 3

    ## SYSTEM HIERARCHY
    # # sec_atoms_group = sec_systems[0].model_system
    # # assert len(sec_atoms_group) == 4
    # # assert sec_atoms_group[0].branch_label == 'group_1ZNF'
    # # assert sec_atoms_group[0].atom_indices[159] == 159
    # # sec_proteins = sec_atoms_group[0].model_system
    # # assert len(sec_proteins) == 1
    # # assert sec_proteins[0].branch_label == '1ZNF'
    # # assert sec_proteins[0].atom_indices[400] == 400
    # # sec_res_group = sec_proteins[0].model_system
    # # assert len(sec_res_group) == 16
    # # assert sec_res_group[14].branch_label == 'group_SER'
    # # assert sec_res_group[14].atom_indices[2] == 136
    # # sec_res = sec_res_group[14].model_system
    # # assert len(sec_res) == 3
    # # assert sec_res[0].branch_label == 'SER'
    # # assert sec_res[0].atom_indices[10] == 144
    # # assert sec_res[0].custom_system_attributes[0].name == 'hydrophobicity'
    # # assert sec_res[0].custom_system_attributes[0].value == '0.13'
    # # assert sec_res[0].custom_system_attributes[0].unit is None

    ## OUTPUTS
    sec_outputs = sec_simulation.outputs
    assert len(sec_outputs) == 5
    assert sec_outputs[3].step == 3
    assert sec_outputs[2].time.to('ps').magnitude == approx(2.0)
    # Temperature
    assert sec_outputs[2].temperatures[0].value.to('kelvin').magnitude == approx(300.0)
    # Energies
    print(sec_outputs[2].total_energies[0])
    print(sec_outputs[2].total_energies[0].value)
    assert sec_outputs[2].total_energies[0].value.to('kilojoule').magnitude == approx(
        6.0
    )
    assert sec_outputs[2].total_energies[0].contributions[0].name == 'custom'
    assert sec_outputs[2].total_energies[0].contributions[0].value.to(
        'kilojoule'
    ).magnitude == approx(3.0)
    assert sec_outputs[2].total_energies[0].contributions[1].name == 'kinetic'
    assert sec_outputs[2].total_energies[0].contributions[1].value.to(
        'kilojoule'
    ).magnitude == approx(2.0)
    assert sec_outputs[2].total_energies[0].contributions[2].name == 'potential'
    assert sec_outputs[2].total_energies[0].contributions[2].value.to(
        'kilojoule'
    ).magnitude == approx(1.0)
    # Forces
    assert np.shape(sec_outputs[1].total_forces[0].value) == (31583, 3)
    assert sec_outputs[1].total_forces[0].value[2100][2].to(
        'newton'
    ).magnitude == approx(500.0)
    assert sec_outputs[2].total_forces[0].value[11].to('newton').magnitude == approx(
        500.0
    )
    assert sec_outputs[2].total_forces[0].contributions[0].name == 'custom'
    assert sec_outputs[2].total_forces[0].contributions[0].value[21].to(
        'newton'
    ).magnitude == approx(4.0)
    # Custom Outputs
    assert sec_outputs[2].custom_outputs[0].m_def.name == 'CustomProperty'
    assert len(sec_outputs[1].custom_outputs) == 1
    assert sec_outputs[1].custom_outputs[0].name == 'custom_thermodynamic_properties'
    assert sec_outputs[1].custom_outputs[0].value == approx(100.0)
    assert sec_outputs[1].custom_outputs[0].unit == 'newton / angstrom ** 2'

    ## WORKFLOW
    sec_workflow = archive.workflow2
    # MD method
    assert sec_workflow.method.integrator_type == 'langevin_leap_frog'
    assert sec_workflow.method.thermodynamic_ensemble == 'NPT'
    assert sec_workflow.method.integration_timestep.to(
        'picosecond'
    ).magnitude == approx(2e-15)
    assert sec_workflow.method.n_steps == 20000000
    assert sec_workflow.method.coordinate_save_frequency == 10000
    assert sec_workflow.method.velocity_save_frequency == None
    assert sec_workflow.method.force_save_frequency == None
    assert sec_workflow.method.thermodynamics_save_frequency == None
    # MD thermostat
    sec_thermostat = sec_workflow.method.thermostat_parameters
    assert sec_thermostat[0].thermostat_type == 'langevin_leap_frog'
    assert sec_thermostat[0].reference_temperature.magnitude == approx(300.0)
    assert sec_thermostat[0].coupling_constant.to('picosecond').magnitude == approx(1.0)
    assert sec_thermostat[0].effective_mass == None
    assert sec_thermostat[0].temperature_profile == None
    assert sec_thermostat[0].reference_temperature_start == None
    assert sec_thermostat[0].reference_temperature_end == None
    assert sec_thermostat[0].temperature_update_frequency == None
    assert sec_thermostat[0].temperature_update_delta == None
    assert sec_thermostat[0].temperature_update_factor == None
    assert sec_thermostat[0].step_start == None
    assert sec_thermostat[0].step_end == None
    # MD barostat
    sec_barostat = sec_workflow.method.barostat_parameters
    assert sec_barostat[0].barostat_type == 'berendsen'
    assert sec_barostat[0].coupling_type == 'isotropic'
    assert np.all(
        sec_barostat[0].reference_pressure.to('bar').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.all(
        sec_barostat[0].coupling_constant.to('picosecond').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.all(
        sec_barostat[0].compressibility.to('1/bar').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert sec_barostat[0].pressure_profile == None
    assert sec_barostat[0].reference_pressure_start == None
    assert sec_barostat[0].reference_pressure_end == None
    assert sec_barostat[0].pressure_update_frequency == None
    assert sec_barostat[0].pressure_update_delta == None
    assert sec_barostat[0].pressure_update_factor == None
    assert sec_barostat[0].step_start == None
    assert sec_barostat[0].step_end == None
    # MD Shear
    sec_shear = sec_workflow.method.shear_parameters
    assert sec_shear == []
    # assert sec_shear[0].shear_type == None
    # assert sec_shear[0].shear_rate == None
    # assert sec_shear[0].step_start == None
    # assert sec_shear[0].step_end == None
    # MD Free Energy Calculation Parameters
    sec_free_energy = sec_workflow.method.free_energy_calculation_parameters
    assert sec_free_energy == []
    # assert sec_free_energy[0].type == None
    # assert sec_free_energy[0].lambda_index == None
    # assert sec_free_energy[0].atom_indices == None
    # assert sec_free_energy[0].initial_state_vdw == None
    # assert sec_free_energy[0].final_state_vdw == None
    # assert sec_free_energy[0].initial_state_coloumb == None
    # assert sec_free_energy[0].final_state_coloumb == None
    # assert sec_free_energy[0].initial_state_bonded == None
    # assert sec_free_energy[0].final_state_bonded == None
    # sec_lambdas = sec_free_energy[0].lambdas
    # assert sec_lambdas[0].type == None
    # assert sec_lambdas[0].value == None
