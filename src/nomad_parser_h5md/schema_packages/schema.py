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


class CustomProperty(physical_property.PhysicalProperty):
    """
    Section describing a general type of calculation.
    """

    # m_def = Section('get_custom_outputs'...a_mapping=)
    value = Quantity(
        type=np.float64,
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
    properties.energies.EnergyContribution.name.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='.name')
    pass

    # value annotation defined in TotalEnergy.value since they refer to the same quantity
    # in this case, we make sure to return the corresponding value from
    # the get_contributions function in the TotalEnergy.contributions annotation


class TotalEnergy(properties.TotalEnergy):
    properties.TotalEnergy.value.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(
            mapper=('get_output_data', ['.@'], dict(path='observables.energies.total'))
        )
    )

    properties.energies.TotalEnergy.contributions.m_annotations.setdefault(
        'mapping', {}
    )['hdf5'] = MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.energies', exclude=['total']),
        )
    )


class Temperature(properties.Temperature):
    properties.Temperature.value.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(
            mapper=('get_output_data', ['.@'], dict(path='observables.temperatures'))
        )
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


# class Outputs(outputs.Outputs):
#     outputs.Outputs.total_energies.m_annotations.setdefault('mapping', {})['hdf5'] = (
#         MapperAnnotation(mapper='.@')
#     )

# outputs.Outputs.total_forces.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper='.@')
# )


class AtomsState(atoms_state.AtomsState):
    atoms_state.AtomsState.chemical_symbol.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='.label')


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.positions.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='.positions')

    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='.lattice_vectors')

    model_system.AtomicCell.velocities.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='.velocities')

    # TODO length of positions in section data does not work
    model_system.AtomicCell.n_atoms.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(mapper='length(particles.all.position.value.__value | [0])')
    )

    model_system.AtomicCell.atoms_state.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper=('to_species_labels', ['particles.all.species_label']))

    model_system.AtomicCell.periodic_boundary_conditions.m_annotations.setdefault(
        'mapping', {}
    )['hdf5'] = MapperAnnotation(mapper='.boundary')


class ModelSystem(model_system.ModelSystem):
    """
    Model system used as an input for simulating the material.
    """

    m_def = Section(
        validate=False,
        # extends_base_section=True,
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

    model_system.AtomicCell.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(mapper=('get_system_data', ['.@']))
    )

    # TODO inconsistent? shape with original def
    model_system.ModelSystem.bond_list.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper='connectivity.bonds')

    model_system.ModelSystem.dimensionality.m_annotations.setdefault('mapping', {})[
        'hdf5'
    ] = MapperAnnotation(mapper=r'particles.all.box."@dimension"')


class TrajectoryOutputs(outputs.TrajectoryOutputs):
    m_def = Section(
        validate=False,
        # extends_base_section=True,
        # a_mapping=dict(hdf5=MapperAnnotation(mapper='.@')),
    )

    # step = Quantity(
    #     type=np.int32,
    #     description="""
    #     The step number with respect to the workflow.
    #     """,
    # )

    outputs.TrajectoryOutputs.step.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(
            mapper=(
                'set_step',
                ['.@'],
                dict(path='observables.custom_forces'),
            )
        )
    )
    #     MapperAnnotation(mapper='.step')
    # )

    # outputs.TrajectoryOutputs.time.m_annotations.setdefault('mapping', {})['hdf5'] = (
    #     MapperAnnotation(mapper='.time')
    # )

    # outputs.TrajectoryOutputs.total_energies.m_annotations.setdefault('mapping', {})[
    #     'hdf5'
    # ] = MapperAnnotation(mapper='.@')

    # outputs.TrajectoryOutputs.temperatures.m_annotations.setdefault('mapping', {})[
    #     'hdf5'
    # ] = MapperAnnotation(mapper='.@')

    # custom_outputs = SubSection(
    #     sub_section=CustomProperty.m_def,
    #     description="""
    #     Contains other generic custom outputs that are not already defined.
    #     """,
    #     repeats=True,
    # )

    # custom_outputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    #     MapperAnnotation(
    #         mapper=(
    #             'get_custom_outputs',
    #             ['.@'],
    #             dict(
    #                 path='observables',
    #                 exclude=[
    #                     'energies, temperatures, custom_forces'
    #                 ],  # TODO get the exclusion list automatically
    #                 observable_type='configurational',
    #             ),
    #         )
    #     )
    # )

    # custom_outputs.value.m_annotations.setdefault('mapping', {})['hdf5'] = (

    # custom_outputs = SubSection(sub_section=..., a_mapping=dict(hdf5=MapperAnnotation(
    #         mapper=(
    #             'get_custom_outputs',
    #             ['.@'],
    #             dict(
    #                 path='observables',
    #                 exclude=[
    #                     'energies, temperatures, custom_forces'
    #                 ],  # TODO get the exclusion list automatically
    #                 observable_type='configurational',
    #             ),
    #         )
    #     ))

    # custom_outputs.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    #     mapper=(
    #         'get_custom_outputs',
    #         ['.@'],
    #         dict(
    #             path='observables',
    #             exclude=[
    #                 'energies, temperatures, custom_forces'
    #             ],  # TODO get the exclusion list automatically
    #             observable_type='configurational',
    #         ),
    #     )
    # )


# outputs.TrajectoryOutputs.custom_outputs.m_annotations.setdefault('mapping', {})[
#     'hdf5'
# ] = MapperAnnotation(mapper='.@')


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


# class H5MDCreator(general.Program):
#     """
#     Contains the specifications of the program.
#     """

#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     name = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Specifies the name of the author who generated the h5md file.
#         """,
#     )

#     email = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Author's email.
#         """,
#     )


class Program(general.Program):
    general.Program.name.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(
            mapper='."@name"',
        )
    )

    general.Program.version.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(
            mapper='."@version"',
        )
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

    general.Simulation.program.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(mapper='h5md.program')
    )

    # general.Simulation.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
    #     MapperAnnotation(mapper=('get_system_steps', ['particles.all.position']))
    # )

    outputs.TrajectoryOutputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
        MapperAnnotation(mapper=('get_output_steps', ['observables']))
    )


Simulation.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='@'
)


m_package.__init_metainfo__()
