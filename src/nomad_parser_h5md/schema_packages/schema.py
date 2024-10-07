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

from nomad_simulations.schema_packages import general, physical_property, outputs, model_system, atoms_state
import numpy as np
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import Context, MEnum, Quantity, Section, SectionProxy, SubSection
from nomad.parsing.file_parser.mapping_parser import MappingAnnotationModel

from nomad.metainfo import Quantity, SchemaPackage

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


# class CustomProperty(physical_property.PhysicalProperty):
#     """
#     Section describing a general type of calculation.
#     """

#     value = Quantity(
#         type=np.float64,
#         shape=[],
#         description="""
#         Value **magnitude** of the property. The unit is defined in the `unit` attribute.
#         """,
#     )

#     unit = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Unit of the parameter as a string consistent with the UnitRegistry.pint module.
#         """,
#     )

#     description = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Further description of the property.
#         """,
#     )


# class EnergyEntry(ArchiveSection):
#     """
#     Section describing a general type of energy contribution.
#     """

#     name = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Name of the energy contribution.
#         """,
#     )

#     value = Quantity(
#         type=np.dtype(np.float64),
#         shape=[],
#         unit='joule',
#         description="""
#         Value of the energy contribution.
#         """,
#     )


# class ForceEntry(ArchiveSection):
#     """
#     Section describing a general type of force contribution.
#     """

#     name = Quantity(
#         type=str,
#         shape=[],
#         description="""
#         Name of the force contribution.
#         """,
#     )

#     value = Quantity(
#         type=np.dtype(np.float64),
#         shape=[],
#         unit='newton',
#         description="""
#         Value of the force contribution.
#         """,
#     )


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

class AtomsState(atoms_state.AtomsState):
    atoms_state.AtomsState.chemical_symbol.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.label')


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.positions.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.positions')

    model_system.AtomicCell.lattice_vectors.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.lattice_vectors')

    model_system.AtomicCell.velocities.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.velocities')

    # TODO length of positions in section data does not work
    model_system.AtomicCell.n_atoms.m_annotations['hdf5'] = MappingAnnotationModel(mapper='length(particles.all.position.value.__value | [0])')

    model_system.AtomicCell.atoms_state.m_annotations['hdf5'] = MappingAnnotationModel(mapper=('to_species_labels', ['particles.all.species_label']))


class ModelSystem(model_system.ModelSystem):
    """
    Model system used as an input for simulating the material.
    """

    m_def = Section(
        validate=False,
        extends_base_section=True,
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

    model_system.AtomicCell.m_def.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.@')

    # TODO inconsistent? shape with original def
    bond_list = Quantity(
        type=np.int32,
        shape=['*', '*'],
        description="""
        List of pairs of atom indices corresponding to bonds (e.g., as defined by a force field)
        within this atoms_group.
        """,
    )

    bond_list.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.bonds_list')

# class Stress(physical_property.PhysicalProperty):
#     """ """

#     value = Quantity(
#         type=np.dtype(np.float64),
#         unit='newton',
#         description="""
#         """,
#     )

#     def normalize(self, archive: 'EntryArchive', logger: 'BoundLogger') -> None:
#         super().normalize(archive, logger)


# class TrajectoryOutputs(outputs.TrajectoryOutputs):
#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     x_h5md_custom_outputs = SubSection(
#         sub_section=CustomProperty.m_def,
#         description="""
#         Contains other generic custom outputs that are not already defined.
#         """,
#         repeats=True,
#     )


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

    name.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.\"@name\"')

    email = Quantity(
        type=str,
        shape=[],
        description="""
        Author's email.
        """,
    )

    email.m_annotations['hdf5'] = MappingAnnotationModel(mapper='.\"@email\"')


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
    general.Program.name.m_annotations['hdf5'] = MappingAnnotationModel(
        mapper='.\"@name\"',
    )

    general.Program.version.m_annotations['hdf5'] = MappingAnnotationModel(
        mapper='.\"@version\"',
    )


class Simulation(general.Simulation):
    m_def = Section(
        validate=False,
        extends_base_section=True,
    )

    # TODO Not sure how we are dealing with versioning with H5MD-NOMAD
    x_h5md_version = Quantity(
        type=np.dtype(np.int32),
        shape=[2],
        description="""
        Specifies the version of the h5md schema being followed.
        """,
    )
    x_h5md_version.m_annotations['hdf5'] = MappingAnnotationModel(
        mapper='h5md.\"@version\"',
    )

    x_h5md_author = SubSection(sub_section=Author.m_def)

    x_h5md_author.m_annotations['hdf5'] = MappingAnnotationModel(mapper='h5md.author')

    x_h5md_creator = SubSection(sub_section=general.Program.m_def)

    x_h5md_creator.m_annotations['hdf5'] = MappingAnnotationModel(mapper='h5md.creator')

    general.Simulation.program.m_annotations['hdf5'] = MappingAnnotationModel(mapper='h5md.program')

    general.Simulation.model_system.m_annotations['hdf5'] = MappingAnnotationModel(mapper=('to_systems_data', ['@']))

Simulation.m_def.m_annotations['hdf5'] = MappingAnnotationModel(mapper='@')


m_package.__init_metainfo__()