from typing import Any
import pint

from nomad.parsing.file_parser.mapping_parser import HDF5Parser, MetainfoParser, Path
from nomad_parser_h5md.schema_packages.schema import (
    Simulation,
    ParamEntry,
)
from nomad_parser_h5md.schema_packages.schema import ModelSystem

# from nomad_simulations.schema_packages.model_system import ModelSystem
from simulationworkflowschema.molecular_dynamics import MolecularDynamics
from nomad_parser_h5md.parsers.mdparserutils import MDParser
from nomad.units import ureg

from nomad_parser_h5md.parsers.utils import remove_mapping_annotations

from h5py import Group


class H5MDH5Parser(HDF5Parser):
    trajectory_steps: list[int] = []
    output_steps: list[int] = []
    observables: dict[str, Any] = {}

    def get_value(self, name: str, dct: dict[str, Any]) -> Any:
        value = dct.get(name, {})
        if not isinstance(value, dict):
            return value
        value = value.get(self.value_key)
        if value is None:
            return
        unit = dct.get(name, {}).get(f'{self.attribute_prefix}unit')
        if unit:
            value = value * ureg(unit)
        factor = dct.get(name, {}).get(f'{self.attribute_prefix}unit_factor')
        if factor:
            value = value * factor
        return value

    def get_source(self, parent: dict[str, Any], path: str):
        path_segments = path.split('.', 1)
        source = parent.get(path_segments[0], {})
        if len(path_segments) == 1:
            return source
        return self.get_source(source, path_segments[1])

    def map_value(self, source: dict[str, Any], **kwargs) -> Any:
        if kwargs.get('key') is None:
            return None

        value = self.get_value(kwargs.get('key'), source)

        enum_spec = kwargs.get('enum_spec', None)
        if enum_spec == 'upper':
            return value.upper() if isinstance(value, str) else value
        elif enum_spec == 'lower':
            return value.lower() if isinstance(value, str) else value

        return value

    def get_sub_systems(self, source: dict[str, Any], **kwargs) -> list[dict[str, Any]]:
        # print('in get_sub_systems')
        step = source.get('step', None)
        label = source.get('label', None)
        path = kwargs.get('path', None)
        # print(f'step: {step}, path: {path}, label: {label}')
        # if step is None or path is None:
        #     return []
        if step is not None:
            if step != 0:  # TODO extend to time-dependent bond lists and topologies
                return []
            source = self.get_source(self.data, kwargs['path'])
            # source = [val for _, val in source.items()]

        particles_group = source.get('particles_group', None)
        # print(f'particles_group: {particles_group}')
        if particles_group is None:
            return []

        source = (
            [group for group in particles_group.values()]
            if isinstance(particles_group, dict)
            else []
        )

        return source
        # print(f'type(source): {type(source)}')
        # print(f'source: {source}')
        # # print(f'source.keys(): {source.keys()}')
        # for key, val in source.items():
        #     print(f'key: {key}, val.keys: {val.keys()}, label: {val.get("label")}')

        # convert source_data / particles_group to a list of dicts recursively
        # source_data = [val for _, val in source_data.items()]
        # for item in source_data:
        #     particles_group = item.pop('particles_group', None)
        #     if particles_group:
        #         item[]

        # return [
        #     {
        #         'label': 'group_1',
        #         'formula': 'form(1)',
        #         'particles_group': [
        #             {'label': 'mol_1', 'formula': 'mol(1)'},
        #             {'label': 'mol_2', 'formula': 'mol(2)'},
        #         ],
        #     },
        #     {'label': 'group_2', 'formula': 'form(2)'},
        #     {'label': 'group_3', 'formula': 'form(3)'},
        # ]
        # return [
        #     [
        #         {
        #             'label': 'group_1',
        #             'formula': 'form(1)',
        #             # 'particles_group': [
        #             #     {'label': 'mol_1', 'formula': 'mol(1)'},
        #             #     {'label': 'mol_2', 'formula': 'mol(2)'},
        #             # ],
        #         },
        #         {
        #             'label': 'group_1-2',
        #             'formula': 'form(1-2)',
        #         },
        #     ],
        #     {'label': 'group_2', 'formula': 'form(2)'},
        #     {'label': 'group_3', 'formula': 'form(3)'},
        # ]

    # def get_system_hierarchy(
    #     self,
    #     particlesgroup: {},
    # ) -> list[dict[str, Any]]:
    #     data = []
    #     for key, dct in particlesgroup.items():
    #         data.append(dct)
    #         path_particlesgroup_key = f'{path_particlesgroup}.{key}'

    #         particles_group = {
    #             group_key: h5md_sec_particlesgroup.get(
    #                 f'{path_particlesgroup_key}.{group_key}'
    #             )
    #             for group_key in h5md_sec_particlesgroup[key].keys()
    #         }

    #         particles_group = {
    #             group_key: self.data.get(f'{path_particlesgroup_key}.{group_key}')
    #             for group_key in h5md_sec_particlesgroup[key].keys()
    #         }
    #         data['branch_label'] = particles_group.pop('label', None)
    #         data['atom_indices'] = particles_group.pop('indices', None)
    #         # TODO remove the deprecated below from the test file
    #         data['type'] = particles_group.pop('type', None)  # ? deprecate?
    #         data['is_molecule'] = particles_group.pop(
    #             'is_molecule', None
    #         )  # ? deprecate?
    #         particles_group.pop('formula', None)  # covered in normalization now
    #         # write all the standard quantities to the archive
    #         particles_subgroup = particles_group.pop('particles_group', None)

    #         # set the remaining attributes
    #         data['custom_system_attributes'] = []
    #         for particles_group_key in particles_group.keys():
    #             val = particles_group.get(particles_group_key)
    #             units = val.units if hasattr(val, 'units') else None
    #             val = val.magnitude if units is not None else val
    #             data['custom_system_attributes'].append(
    #                 {'name': particles_group_key, 'value': val, 'unit': units}
    #             )

    #         # get the next branch level
    #         if particles_subgroup:
    #             self.get_system_hierarchy(
    #                 particles_subgroup,
    #                 f'{path_particlesgroup_key}.particles_group',
    #             )

    #     return []

    def get_system_steps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        steps = self.get_value('step', source)
        times = self.get_value('time', source)
        system_steps = [
            {'step': step, 'time': times[n]}
            for n, step in enumerate(steps)
            if step in self.trajectory_steps
        ]

        # get system hierarchy and store in first step
        # h5md_sec_particlesgroup = self.data.get('connectivity', {}).get(
        #     'particles_group'
        # )
        # hierarchy = self.get_system_hierarchy(
        #     h5md_sec_particlesgroup=h5md_sec_particlesgroup,
        #     path_particlesgroup='connectivity.particles_group',
        # )
        # system_steps[system_steps.keys()[0]]['model_system'] = hierarchy
        return system_steps

    def get_step_data(self, data: dict[str, Any], step: int) -> dict[str, Any]:
        step_data = {}
        value = self.get_value('value', data)
        steps = self.get_value('step', data)
        if value is None or steps is None:
            return step_data
        index = steps.index(step)
        step_data['value'] = value[index]
        times = self.get_value('time', data)
        step_data['time'] = times[index]
        return step_data

    def get_cell_data(self, source: dict[str, Any]) -> dict[str, Any]:
        particles = self.data.get('particles', {}).get('all')
        if particles is None:
            return {}

        source_paths = [
            ('lattice_vectors', 'box.edges'),
        ]
        system_data = {}
        for name, path in source_paths:
            data = self.get_source(particles, path)
            if data is None:
                continue
            step_data = self.get_step_data(data, source.get('step'))
            if not step_data:
                continue
            system_data.setdefault('time', step_data['time'])
            if step_data['time'] == system_data['time']:
                system_data[name] = step_data['value']
        box = particles.get('box', {})
        system_data['boundary'] = box.get(f'{self.attribute_prefix}boundary')

        return system_data

    def get_traj_data(self, source: dict[str, Any], **kwargs) -> pint.Quantity:
        # print('in get_traj_data')
        # print(source.get('step'))
        # print(source.keys())
        # print(kwargs.get('path'))
        # if source.get('step') is None:
        #     return

        # print(self.get_step_data(source, source['step']).get('value'))
        # return self.get_step_data(source, source['step']).get('value')
        # TODO - check to see if this function can be combined with get_output_data
        if source.get('value') is not None:
            return source['value']
        if source.get('step') is None or kwargs.get('path') is None:
            return

        source_data = self.get_source(self.data, kwargs['path'])

        # print(f'source_data.keys(): {source_data.keys()}')
        # print(f'len(source_data.value): {len(source_data.get("value"))}')
        # import numpy as np

        # print(np.array(source_data.get('value')).shape)
        # # print(f'len(source_data.value[0]): {len(source_data.get("value")[0])}')
        # # print(f'source_data.value: {source_data.get("value")}')
        # print(
        #     f'self.get_step_data(source_data, source["step"]): {self.get_step_data(source_data, source["step"])}'
        # )
        data = self.get_step_data(source_data, source['step']).get('value')
        # print(f'data: {data}')
        # print(f'len(data): {len(data)}')
        # print(f'np.array(data).shape: {np.array(data.magnitude).shape}')
        # print(f'type(data): {type(data.magnitude)}')
        return data

    def get_system_data(self, source: dict[str, Any]) -> dict[str, Any]:
        # print('in get_system_data')
        particles = self.data.get('particles', {}).get('all')
        if particles is None:
            return {}

        source_paths = [
            ('positions', 'position'),
            ('velocities', 'velocity'),
        ]
        system_data = {}
        # TODO - could extract this into a helper function also for get_cell_data
        for name, path in source_paths:
            data = self.get_source(particles, path)
            if data is None:
                continue
            step_data = self.get_step_data(data, source.get('step'))
            if not step_data:
                continue
            system_data.setdefault('time', step_data['time'])
            if step_data['time'] == system_data['time']:
                system_data[name] = step_data['value']

        return system_data

    def to_species_labels(self, source: list[str]) -> list[dict[str, Any]]:
        return [{'chemical_symbol': s, 'label': s} for s in source]

    def get_output_steps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        output_steps = {}

        def get_steps(dct: dict[str, Any]) -> dict[str, Any]:
            steps = self.get_value('step', dct)
            if steps is None:
                return {}
            times = self.get_value('time', dct)
            if len(steps) != len(times):
                self.logger.error(
                    'Inconsistent step-time combinations in observable data.'
                )

            return {step: times[n] for n, step in enumerate(steps)}

        def get_observable_steps(
            source: dict[str, Any], output_steps: dict[str, Any]
        ) -> None:
            for __, val in source.items():
                observable_type = val.get('@type')
                if not observable_type:
                    get_observable_steps(val, output_steps)
                elif observable_type == 'configurational':
                    steps = get_steps(val)
                    if not all(
                        [
                            time == output_steps[step]
                            for step, time in steps.items()
                            if step in output_steps.keys()
                        ]
                    ):
                        self.logger.error(
                            'Inconsistent step-time combinations in observable data.'
                        )
                    output_steps.update(steps)

        get_observable_steps(source, output_steps)
        output_steps = [
            {'step': step, 'time': time} for step, time in output_steps.items()
        ]
        return output_steps

    def get_contributions(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if kwargs.get('path') is None or source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])
        include = kwargs.get('include')
        exclude = kwargs.get('exclude')
        contributions = []
        for key, val in source_data.items():
            if include and key not in include or exclude and key in exclude:
                continue
            step_data = self.get_step_data(val, source['step'])
            contributions.append({'name': key, **step_data})
        return contributions

    def get_output_data(self, source: dict[str, Any], **kwargs) -> pint.Quantity:
        if source.get('value') is not None:
            return source['value']
        if source.get('step') is None or kwargs.get('path') is None:
            return
        observable_type = kwargs.get('observable_type')
        if observable_type is None or observable_type not in [
            'configurational',
            'ensemble_average',
            'correlation_function',
        ]:
            self.logger.warning(
                'Invalid or no obervable type defined in the schema annotation '
                f'for {source.keys()},skipping this observable.'
            )
            return

        source_data = self.get_source(self.data, kwargs['path'])
        if source_data.get('@type') != observable_type:
            return

        data = self.get_step_data(source_data, source['step']).get('value')
        # print('in get_output_data')
        # print(f'source_step: {source.get("step")}')
        # print(f'output data: {data}')
        # print(f'len(data): {len(data)}')
        # import numpy as np

        # print(f'np.array(data).shape: {np.array(data.magnitude).shape}')
        # print(f'type(data): {type(data.magnitude)}')
        return data

    def get_custom_outputs(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if kwargs.get('path') is None or source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])
        include = kwargs.get('include')
        exclude = kwargs.get('exclude')
        observable_type = kwargs.get('observable_type')
        custom_outputs = []
        for key, val in source_data.items():
            if include and key not in include or exclude and key in exclude:
                continue
            if observable_type is not None:
                source_type = val.get('@type')
                if source_type != observable_type:
                    continue
            step_data = self.get_step_data(val, source['step'])
            if step_data.get('value') is not None:
                if isinstance(step_data['value'], pint.Quantity):
                    step_data['unit'] = str(step_data['value'].units)
                    step_data['value'] = step_data['value'].magnitude
            custom_outputs.append({'name': key, **step_data})
        return custom_outputs

    # def get_parameters(self, parameter_group: Group, path: str) -> dict:
    #     param_dict: dict[Any, Any] = {}
    #     for key, val in parameter_group.items():
    #         path_key = f'{path}.{key}'
    #         if isinstance(val, Group):
    #             param_dict[key] = self.get_parameters(val, path_key)
    #         else:
    #             param_dict[key] = self._data_parser.get(path_key)
    #             if isinstance(param_dict[key], str):
    #                 param_dict[key] = (
    #                     param_dict[key].upper()
    #                     if key == 'thermodynamic_ensemble'
    #                     else param_dict[key].lower()
    #                 )
    #             elif isinstance(param_dict[key], (int, np.int32, np.int64)):
    #                 param_dict[key] = param_dict[key].item()
    #     return param_dict

    def get_md_parameters(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if kwargs.get('path') is None:
            return []

        return []
        # source_data = self.get_source(self.data, kwargs['path'])

        # self._parameter_info = {'force_calculations': {}, 'workflow': {}}

        # force_calculations_group = self._data_parser.get(
        #     'parameters.force_calculations'
        # )
        # if force_calculations_group is not None:
        #     self._parameter_info['force_calculations'] = self.get_parameters(
        #         force_calculations_group, 'parameters.force_calculations'
        #     )
        # workflow_group = self._data_parser.get('parameters.workflow')
        # if workflow_group is not None:
        #     self._parameter_info['workflow'] = self.get_parameters(
        #         workflow_group, 'parameters.workflow'
        #     )


class H5MDParser(MDParser):
    def __init__(self) -> None:
        super().__init__()
        self.h5_parser = H5MDH5Parser()
        self.simulation_parser = MetainfoParser()
        self.simulation_parser.max_nested_level = 10
        self.workflow_parser = MetainfoParser()

    def parse_system_hierarchy(
        self,
        nomad_sec: ModelSystem,
        h5md_sec_particlesgroup: Group,
        path_particlesgroup: str,
    ):
        data = {}
        for key in h5md_sec_particlesgroup.keys():
            path_particlesgroup_key = f'{path_particlesgroup}.{key}'
            particles_group = {
                group_key: self._data_parser.get(
                    f'{path_particlesgroup_key}.{group_key}'
                )
                for group_key in h5md_sec_particlesgroup[key].keys()
            }
            sec_model_system = ModelSystem()
            nomad_sec.model_system.append(sec_model_system)
            data['branch_label'] = particles_group.pop('label', None)
            data['atom_indices'] = particles_group.pop('indices', None)
            # TODO remove the deprecated below from the test file
            # sec_atomsgroup.type = particles_group.pop("type", None) #? deprecate?
            particles_group.pop('type', None)
            # sec_atomsgroup.is_molecule = particles_group.pop("is_molecule", None) #? deprecate?
            particles_group.pop('is_molecule', None)
            particles_group.pop('formula', None)  # covered in normalization now
            # write all the standard quantities to the archive
            print(sec_model_system.sub_systems)
            self.parse_section(data, sec_model_system.sub_systems)
            print(sec_model_system.sub_systems)
            particles_subgroup = particles_group.pop('particles_group', None)

            # set the remaining attributes
            for particles_group_key in particles_group.keys():
                val = particles_group.get(particles_group_key)
                units = val.units if hasattr(val, 'units') else None
                val = val.magnitude if units is not None else val
                sec_model_system.custom_system_attributes.append(
                    ParamEntry(name=particles_group_key, value=val, unit=units)
                )

            # get the next branch level
            if particles_subgroup:
                self.parse_system_hierarchy(
                    sec_model_system,
                    particles_subgroup,
                    f'{path_particlesgroup_key}.particles_group',
                )

    def write_to_archive(self) -> None:
        # create h5 parser
        self.h5_parser.filepath = self.mainfile

        self.trajectory_steps = Path(path='particles.all.position.step').get_data(
            self.h5_parser.data, default=[]
        )
        self.h5_parser.trajectory_steps = self.trajectory_steps

        # TODO consider using a single parser for the whole archive
        # create metainfo parsers
        self.simulation_parser.annotation_key = 'hdf5'
        simulation_data = Simulation()
        self.simulation_parser.data_object = simulation_data
        self.workflow_parser.annotation_key = 'hdf5'
        workflow_data = MolecularDynamics()
        self.workflow_parser.data_object = workflow_data

        # map from h5 source to metainfo target
        self.h5_parser.convert(self.simulation_parser)
        self.h5_parser.convert(self.workflow_parser)

        # assign simulation to archive data
        self.archive.data = self.simulation_parser.data_object
        self.archive.workflow2 = self.workflow_parser.data_object

        # manually set the hierarchy for now
        # hierarchy_root = Path(path='connectivity.particles_group').get_data(
        #     self.h5_parser.data, default=[]
        # )
        # if hierarchy_root:
        #     print('in hierarchy root')
        #     print(type(hierarchy_root))
        #     self.parse_system_hierarchy(
        #         self.archive.data.model_system,
        #         hierarchy_root,
        #         'connectivity.particles_group',
        #     )

        # close parsers
        self.h5_parser.close()
        self.simulation_parser.close()

        # simulation = self.archive.data
        # print('simulation.m_annotations:', simulation.m_annotations)
        # system = self.archive.data.model_system[0]
        # if not system.sub_systems:
        #     system.sub_systems.append(ModelSystem())
        # print('system.m_annotations:', system.m_annotations)
        # print('system.sub_systems[0].m_def:', system.sub_systems[0].m_def)
        # print(f'system.sub_systems[0]: {system.sub_systems[0]}')
        # print('system.sub_systems[0].__dict__:', system.sub_systems[0].__dict__)
        # print('system.positions.m_annotations:', system.positions.m_annotations)

        # remove mapping annotations
        remove_mapping_annotations(self.archive.data.m_def)
