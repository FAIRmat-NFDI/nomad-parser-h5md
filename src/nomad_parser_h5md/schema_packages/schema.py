#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
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

import numpy as np
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo.data_type import m_float64
from nomad.datamodel.metainfo.annotations import Mapper as MapperAnnotation
from nomad.metainfo import Quantity, SchemaPackage, Section, SubSection
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_system,
    outputs,
    properties,
)
from nomad_simulations.schema_packages import physical_property
from simulationworkflowschema import molecular_dynamics

m_package = SchemaPackage()


# SIMULATION --> archive.data


class Author(ArchiveSection):
    """
    Contains the specifications of the program.
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Specifies the name of the author who generated the h5md file.
        """,
    )

    name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='."@name"'
    )

    email = Quantity(
        type=str,
        shape=[],
        description="""
        Author's email.
        """,
    )

    email.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='."@email"'
    )


class Program(general.Program):
    pass


Program.name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='."@name"',
)

Program.version.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='."@version"',
)

## SIMULATION.MODEL_SYSTEM --> archive.data.model_system


class AtomsState(atoms_state.AtomsState):
    pass


AtomsState.chemical_symbol.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.label')
)


class AtomicCell(model_system.AtomicCell):
    pass


AtomicCell.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_system_data', ['.@'])
)

AtomicCell.positions.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='.positions'
)

AtomicCell.lattice_vectors.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.lattice_vectors')
)

AtomicCell.velocities.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.velocities')
)

# TODO length of positions in section data does not work
AtomicCell.n_atoms.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='length(particles.all.position.value.__value | [0])'
)

AtomicCell.periodic_boundary_conditions.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.boundary')

### SUBSECTIONS

AtomicCell.atoms_state.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('to_species_labels', ['particles.all.species_label']))
)


class ModelSystem(model_system.ModelSystem):
    """
    Model system used as an input for simulating the material.
    """

    m_def = Section(
        validate=False,
    )

    custom_system_attributes = (
        SubSection(  # TODO should this be called parameters or attributes or what?
            sub_section=ParamEntry.m_def,
            description="""
        Contains additional information about the (sub)system .
        """,
            repeats=True,
        )
    )


# TODO inconsistent? shape with original def
ModelSystem.bond_list.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='connectivity.bonds')
)

ModelSystem.dimensionality.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=r'particles.all.box."@dimension"')
)

### SUBSECTIONS

# ModelSystem.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper='.model_system')
# )


# TODO need to add ParticleCell and distinguish in the parser
#### model_system.cell --> AtomicCell

# TODO mabye a note here
#### model_system.model_system --> ModelSystem


## SIMULATION.METHOD --> archive.data.method


class ParamEntry(ArchiveSection):
    """
    Generic section defining a parameter name and value
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the parameter.
        """,
    )

    value = Quantity(
        type=str,
        shape=[],
        description="""
        Value of the parameter as a string.
        """,
    )

    unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the parameter as a string.
        """,
    )

    description = Quantity(
        type=str,
        shape=[],
        description="""
        Further description of the attribute.
        """,
    )


# class ForceCalculations(runschema.method.ForceCalculations):
#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     x_h5md_parameters = SubSection(
#         sub_section=ParamEntry.m_def,
#         description="""
#         Contains non-normalized force calculation parameters.
#         """,
#         repeats=True,
#     )


# class NeighborSearching(runschema.method.NeighborSearching):
#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     x_h5md_parameters = SubSection(
#         sub_section=ParamEntry.m_def,
#         description="""
#         Contains non-normalized neighbor searching parameters.
#         """,
#         repeats=True,
#     )

## SIMULATION.OUTPUTS --> archive.data.outputs


class CustomProperty(ArchiveSection):
    """
    Section describing a general type of calculation.
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the parameter.
        """,
    )
    name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.name'
    )

    value = Quantity(
        type=m_float64(dtype=np.float64).no_shape_check(),
        shape=[],
        description="""
        Value **magnitude** of the property. The unit is defined in the `unit` attribute.
        """,
    )
    value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.value'
    )

    unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the parameter as a string consistent with the UnitRegistry.pint module.
        """,
    )
    unit.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.unit'
    )

    description = Quantity(
        type=str,
        shape=[],
        description="""
        Further description of the property.
        """,
    )
    description.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.description'
    )


class EnergyContribution(properties.energies.EnergyContribution):
    pass


# value annotation defined in TotalEnergy.value since they refer to the same quantity
# in this case, we make sure to return the corresponding value from
# the get_contributions function in the TotalEnergy.contributions annotation
EnergyContribution.name.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.name')
)


class TotalEnergy(properties.TotalEnergy):
    pass


TotalEnergy.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_output_data', ['.@'], dict(path='observables.energies.total'))
)

### SUBSECTIONS

TotalEnergy.contributions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.energies', exclude=['total']),
        )
    )
)


# TODO move these base definitions to nomad-simulations
class ForceContribution(ArchiveSection):
    """
    Abstract class used to define a common `value` quantity with the appropriate units
    for different types of forces, which avoids repeating the definitions for each
    force class.
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the parameter.
        """,
    )

    value = Quantity(
        type=np.float64,
        shape=['*', 3],
        unit='newton',
        description="""
        """,
    )

    def normalize(self, archive: 'EntryArchive', logger: 'BoundLogger') -> None:
        super().normalize(archive, logger)


ForceContribution.name.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.name')
)


class TotalForce(ForceContribution):
    """
    The total force on a system. `contributions` specify individual force
    contributions to the `TotalForce`.
    """

    contributions = SubSection(sub_section=ForceContribution.m_def, repeats=True)

    def __init__(
        self, m_def: 'Section' = None, m_context: 'Context' = None, **kwargs
    ) -> None:
        super().__init__(m_def, m_context, **kwargs)
        self.name = self.m_def.name

    def normalize(self, archive: 'EntryArchive', logger: 'BoundLogger') -> None:
        super().normalize(archive, logger)


TotalForce.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_output_data', ['.@'], dict(path='observables.forces.total'))
)

### SUBSECTIONS

TotalForce.contributions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.forces', exclude=['total']),
        )
    )
)


class Temperature(properties.Temperature):
    pass


Temperature.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_output_data', ['.@'], dict(path='observables.temperatures'))
)


class TrajectoryOutputs(outputs.TrajectoryOutputs):
    m_def = Section(
        validate=False,
    )

    custom_outputs = SubSection(
        sub_section=CustomProperty.m_def,
        description="""
        Contains other generic custom outputs that are not already defined.
        """,
        repeats=True,
    )

    total_forces = SubSection(sub_section=TotalForce.m_def, repeats=True)


TrajectoryOutputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_output_steps', ['observables']))
)

TrajectoryOutputs.step.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.step')
)

TrajectoryOutputs.time.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.time')
)

### SUBSECTIONS

TrajectoryOutputs.total_energies.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.total_forces.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.temperatures.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.custom_outputs.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_custom_outputs',
            ['.@'],
            dict(
                path='observables',
                exclude=[
                    'energies',
                    'temperatures',
                    'custom_forces',
                ],  # TODO get the exclusion list automatically
                observable_type='configurational',
            ),
        )
    )
)


class Simulation(general.Simulation):
    m_def = Section(
        validate=False,
    )

    # TODO Not sure how we are dealing with versioning with H5MD-NOMAD
    x_h5md_version = Quantity(
        type=np.dtype(np.int32),
        shape=[2],
        description="""
        Specifies the version of the h5md schema being followed.
        """,
    )
    x_h5md_version.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md."@version"',
    )

    x_h5md_author = SubSection(sub_section=Author.m_def)

    x_h5md_author.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md.author'
    )

    x_h5md_creator = SubSection(sub_section=general.Program.m_def)

    x_h5md_creator.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md.creator'
    )


Simulation.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='@'
)

Simulation.program.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='h5md.program'
)

Simulation.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_system_steps', ['particles.all.position']))
)

#### Simulation.method --> ??

#### Simulation.outputs --> TrajectoryOutputs


# WORKFLOW --> archive.workflow2


## WORKFLOW.METHOD --> archive.workflow2.method


class ThermostatParameters(molecular_dynamics.ThermostatParameters):
    pass


path_thermostat = f'{path_md}.thermostat_parameters'

ThermostatParameters.thermostat_type.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_thermostat],
            dict(key='thermostat_type'),
        )
    )
)

ThermostatParameters.reference_temperature.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='reference_temperature'),
    )
)

ThermostatParameters.coupling_constant.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='coupling_constant'),
    )
)

ThermostatParameters.effective_mass.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_thermostat],
            dict(key='effective_mass'),
        )
    )
)

ThermostatParameters.temperature_profile.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='temperature_profile'),
    )
)

ThermostatParameters.reference_temperature_start.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='reference_temperature_start'),
    )
)

ThermostatParameters.reference_temperature_end.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='reference_temperature_end'),
    )
)

ThermostatParameters.temperature_update_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='temperature_update_frequency'),
    )
)

ThermostatParameters.temperature_update_delta.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='temperature_update_delta'),
    )
)

ThermostatParameters.temperature_update_factor.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_thermostat],
        dict(key='temperature_update_factor'),
    )
)

ThermostatParameters.step_start.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_thermostat],
            dict(key='step_start'),
        )
    )
)

ThermostatParameters.step_end.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_thermostat],
            dict(key='step_end'),
        )
    )
)


class BarostatParameters(molecular_dynamics.BarostatParameters):
    pass


path_barostat = f'{path_md}.barostat_parameters'

BarostatParameters.barostat_type.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='barostat_type'),
        )
    )
)

BarostatParameters.coupling_type.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='coupling_type'),
        )
    )
)

BarostatParameters.reference_pressure.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='reference_pressure'),
    )
)

BarostatParameters.coupling_constant.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='coupling_constant'),
        )
    )
)

BarostatParameters.compressibility.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='compressibility'),
        )
    )
)

BarostatParameters.pressure_profile.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='pressure_profile'),
        )
    )
)

BarostatParameters.reference_pressure_start.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='reference_pressure_start'),
    )
)

BarostatParameters.reference_pressure_end.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='reference_pressure_end'),
    )
)

BarostatParameters.pressure_update_frequency.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='pressure_update_frequency'),
    )
)

BarostatParameters.pressure_update_delta.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='pressure_update_delta'),
    )
)

BarostatParameters.pressure_update_factor.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_barostat],
        dict(key='pressure_update_factor'),
    )
)

BarostatParameters.step_start.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='step_start'),
        )
    )
)

BarostatParameters.step_end.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_barostat],
            dict(key='step_end'),
        )
    )
)


class ShearParameters(molecular_dynamics.ShearParameters):
    pass


path_shear = f'{path_md}.shear_parameters'

ShearParameters.shear_type.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_shear],
            dict(key='shear_type'),
        )
    )
)

ShearParameters.shear_rate.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_shear],
            dict(key='shear_rate'),
        )
    )
)

ShearParameters.step_start.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_shear],
            dict(key='step_start'),
        )
    )
)

ShearParameters.step_end.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_shear],
            dict(key='step_end'),
        )
    )
)


class FreeEnergyCalculationParameters(
    molecular_dynamics.FreeEnergyCalculationParameters
):
    pass


path_FEC = f'{path_md}.free_energy_calculation_parameters'

# TODO Change this to fec_type in the schema
FreeEnergyCalculationParameters.type.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_shear],
            dict(key='type'),
        )
    )
)

FreeEnergyCalculationParameters.lambda_index.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='lambda_index'),
    )
)

FreeEnergyCalculationParameters.atom_indices.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='atom_indices'),
    )
)

FreeEnergyCalculationParameters.initial_state_vdw.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='initial_state_vdw'),
    )
)

FreeEnergyCalculationParameters.final_state_vdw.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='final_state_vdw'),
    )
)

FreeEnergyCalculationParameters.initial_state_coloumb.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='initial_state_coloumb'),
    )
)

FreeEnergyCalculationParameters.final_state_coloumb.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='final_state_coloumb'),
    )
)

FreeEnergyCalculationParameters.initial_state_bonded.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='initial_state_bonded'),
    )
)

FreeEnergyCalculationParameters.final_state_bonded.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='final_state_bonded'),
    )
)

### SUBSECTIONS

FreeEnergyCalculationParameters.lambdas.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='@')


class Lambdas(molecular_dynamics.Lambdas):
    pass


# TODO lambda_type?
Lambdas.type.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='type'),
    )
)

Lambdas.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_shear],
        dict(key='value'),
    )
)


class MolecularDynamicsMethod(molecular_dynamics.MolecularDynamicsMethod):
    pass


path_md = 'parameters.workflow.molecular_dynamics'

MolecularDynamicsMethod.thermodynamic_ensemble.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='thermodynamic_ensemble'),
    )
)

MolecularDynamicsMethod.integrator_type.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='integrator_type'),
    )
)

MolecularDynamicsMethod.integration_timestep.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='integration_timestep'),
    )
)

MolecularDynamicsMethod.n_steps.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [path_md],
            dict(key='n_steps'),
        )
    )
)

MolecularDynamicsMethod.coordinate_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='coordinate_save_frequency'),
    )
)

MolecularDynamicsMethod.velocity_save_frequency.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='velocity_save_frequency'),
    )
)

MolecularDynamicsMethod.force_save_frequency.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='force_save_frequency'),
    )
)

MolecularDynamicsMethod.thermodynamics_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [path_md],
        dict(key='thermodynamics_save_frequency'),
    )
)

### SUBSECTIONS

MolecularDynamicsMethod.thermostat_parameters.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='@')

MolecularDynamicsMethod.barostat_parameters.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='@')

MolecularDynamicsMethod.shear_parameters.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='@')

MolecularDynamicsMethod.free_energy_calculation_parameters.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')

## WORKFLOW.RESULTS --> archive.workflow2.results


class MolecularDynamicsResults(molecular_dynamics.MolecularDynamicsResults):
    pass


MolecularDynamicsResults.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_output_data', ['observables']))
)

# ? These quantities from normalization?
#     finished_normally = Quantity(
#         type=bool,
#         shape=[],
#         description="""
#         Indicates if calculation terminated normally.
#         """,
#     )

#     n_steps = Quantity(
#         type=np.int32,
#         shape=[],
#         description="""
#         Number of trajectory steps""",
#     )

#     trajectory = Quantity(
#         type=Reference(System),
#         shape=['n_steps'],
#         description="""
#         Reference to the system of each step in the trajectory.
#         """,
#     )

### SUBSECTIONS


# ? Add Custom? OR maybe pull custom out of general schema and put here?
MolecularDynamicsResults.ensemble_properties.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.@')

# TODO This subsection is repeated in the schema
MolecularDynamicsResults.radial_distribution_functions.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.@')

MolecularDynamicsResults.correlation_functions.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.@')

MolecularDynamicsResults.mean_squared_displacements.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.@')

# ? Needed? It just points to the trajectory properties? I guess it collects data here?
MolecularDynamicsResults.radius_of_gyration.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.@')

# ! multi-ensemble property!
MolecularDynamicsResults.free_energy_calculations.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.@')


class MolecularDynamics(molecular_dynamics.MolecularDynamics):
    pass


MolecularDynamics.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

MolecularDynamics.method.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

MolecularDynamics.method.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

# ? Needed?
MolecularDynamics.results.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)
# MolecularDynamics.results.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper=('get_output_data', ['observables']))
# )


m_package.__init_metainfo__()
