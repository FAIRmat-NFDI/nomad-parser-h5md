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

m_package = SchemaPackage()


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


class CustomProperty(ArchiveSection):  # physical_property.PhysicalProperty):
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


class EnergyContribution(properties.energies.EnergyContribution):
    pass


EnergyContribution.name.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.name')
)

# value annotation defined in TotalEnergy.value since they refer to the same quantity
# in this case, we make sure to return the corresponding value from
# the get_contributions function in the TotalEnergy.contributions annotation


class TotalEnergy(properties.TotalEnergy):
    pass


TotalEnergy.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_output_data', ['.@'], dict(path='observables.energies.total'))
)

TotalEnergy.contributions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.energies', exclude=['total']),
        )
    )
)


class Temperature(properties.Temperature):
    pass


Temperature.value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper=('get_output_data', ['.@'], dict(path='observables.temperatures'))
)


# class ForceContribution(properties.forces.ForceContribution):
#     # this is not even necessary as both force and energy name use the same def
#     # it is thus important that the corresponding source data contain name
#     # properties.forces.ForceContribution.name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(mapper='.name')

#     # value annotation defined in TotalForce.value since they refer to the same quantity
#     # in this case, we make sure to return the corresponding value from
#     # the get_contributions function in the TotalForce.contributions annotation
#     pass


# class TotalForce(properties.TotalForce):
#     properties.forces.TotalForce.value.m_annotations.setdefault('mapping', {})[
#         'hdf5'
#     ] = MapperAnnotation(
#         mapper=(
#             'get_output_data',
#             ['.@'],
#             dict(
#                 path='particles.all.force',
#             ),
#         )
#     )

#     properties.forces.TotalForce.contributions.m_annotations.setdefault('mapping', {})[
#         'hdf5'
#     ] = MapperAnnotation(
#         mapper=(
#             'get_contributions',
#             ['.@'],
#             dict(path='observables', include=['custom_forces']),
#         )
#     )


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

AtomicCell.atoms_state.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('to_species_labels', ['particles.all.species_label']))
)

AtomicCell.periodic_boundary_conditions.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.boundary')


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

# ModelSystem.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper=('get_system_hierarchy', ['.@']))
# )


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

TrajectoryOutputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_output_steps', ['observables']))
)

TrajectoryOutputs.step.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.step')
)

TrajectoryOutputs.time.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.time')
)

TrajectoryOutputs.total_energies.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.temperatures.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)


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


class Simulation(general.Simulation):
    m_def = Section(
        validate=False,
        # extends_base_section=True,
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

Simulation.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_system_steps', ['particles.all.position']))
)

Simulation.program.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='h5md.program'
)

m_package.__init_metainfo__()
