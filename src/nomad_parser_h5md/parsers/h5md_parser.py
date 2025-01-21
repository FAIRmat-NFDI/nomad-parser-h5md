from typing import List, Dict, Any
import pint

from nomad.parsing.file_parser.mapping_parser import HDF5Parser, MetainfoParser, Path
from nomad_parser_h5md.schema_packages.schema import Simulation
from nomad_parser_h5md.parsers.mdparserutils import MDParser
from nomad.units import ureg

from h5py import Group


class H5MDH5Parser(HDF5Parser):
    trajectory_steps: List[int] = []
    output_steps: List[int] = []
    observables: Dict[str, Any] = {}

    def get_value(self, name: str, dct: Dict[str, Any]) -> Any:
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

    def get_source(self, parent: Dict[str, Any], path: str):
        path_segments = path.split('.', 1)
        source = parent.get(path_segments[0], {})
        if len(path_segments) == 1:
            return source
        return self.get_source(source, path_segments[1])

    def get_system_steps(self, source: Dict[str, Any]) -> List[Dict[str, Any]]:
        steps = self.get_value('step', source)
        times = self.get_value('time', source)
        return [
            {'step': step, 'time': times[n]}
            for n, step in enumerate(steps)
            if step in self.trajectory_steps
        ]

    def get_step_data(self, data: Dict[str, Any], step: int) -> Dict[str, Any]:
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

    def get_system_data(self, source: Dict[str, Any]) -> Dict[str, Any]:
        particles = self.data.get('particles', {}).get('all')
        if particles is None:
            return {}

        source_paths = [
            ('positions', 'position'),
            ('lattice_vectors', 'box.edges'),
            ('velocities', 'velocity'),
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

    # def get_system_hierarchy(
    #     self,
    #     nomad_sec: ModelSystem,
    #     h5md_sec_particlesgroup: Group,
    #     path_particlesgroup: str,
    # ):
    #     data = {}
    #     for key in h5md_sec_particlesgroup.keys():
    #         path_particlesgroup_key = f'{path_particlesgroup}.{key}'
    #         particles_group = {
    #             group_key: self._data_parser.get(
    #                 f'{path_particlesgroup_key}.{group_key}'
    #             )
    #             for group_key in h5md_sec_particlesgroup[key].keys()
    #         }
    #         sec_model_system = ModelSystem()
    #         nomad_sec.model_system.append(sec_model_system)
    #         data['branch_label'] = particles_group.pop('label', None)
    #         data['atom_indices'] = particles_group.pop('indices', None)
    #         # TODO remove the deprecated below from the test file
    #         # sec_atomsgroup.type = particles_group.pop("type", None) #? deprecate?
    #         particles_group.pop('type', None)
    #         # sec_atomsgroup.is_molecule = particles_group.pop("is_molecule", None) #? deprecate?
    #         particles_group.pop('is_molecule', None)
    #         particles_group.pop('formula', None)  # covered in normalization now
    #         # write all the standard quantities to the archive
    #         self.parse_section(data, sec_model_system)
    #         particles_subgroup = particles_group.pop('particles_group', None)

    #         # set the remaining attributes
    #         for particles_group_key in particles_group.keys():
    #             val = particles_group.get(particles_group_key)
    #             units = val.units if hasattr(val, 'units') else None
    #             val = val.magnitude if units is not None else val
    #             sec_model_system.custom_system_attributes.append(
    #                 ParamEntry(name=particles_group_key, value=val, unit=units)
    #             )

    #         # get the next branch level
    #         if particles_subgroup:
    #             self.parse_system_hierarchy(
    #                 sec_model_system,
    #                 particles_subgroup,
    #                 f'{path_particlesgroup_key}.particles_group',
    #             )

    def to_species_labels(self, source: List[str]) -> List[Dict[str, Any]]:
        return [{'label': s} for s in source]

    def get_output_steps(self, source: Dict[str, Any]) -> List[Dict[str, Any]]:
        output_steps = {}

        def get_steps(dct: Dict[str, Any]) -> Dict[str, Any]:
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
            source: Dict[str, Any], output_steps: Dict[str, Any]
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
        self, source: Dict[str, Any], **kwargs
    ) -> List[Dict[str, Any]]:
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

    def get_output_data(self, source: Dict[str, Any], **kwargs) -> pint.Quantity:
        if source.get('value') is not None:
            return source['value']
        if source.get('step') is None or kwargs.get('path') is None:
            return

        source_data = self.get_source(self.data, kwargs['path'])
        if source_data.get('@type') != 'configurational':
            return
        return self.get_step_data(source_data, source['step']).get('value')

    def get_custom_outputs(
        self, source: Dict[str, Any], **kwargs
    ) -> List[Dict[str, Any]]:
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

    def get_parameters(self, parameter_group: Group, path: str) -> Dict:
        param_dict: Dict[Any, Any] = {}
        for key, val in parameter_group.items():
            path_key = f'{path}.{key}'
            if isinstance(val, Group):
                param_dict[key] = self.get_parameters(val, path_key)
            else:
                param_dict[key] = self._data_parser.get(path_key)
                if isinstance(param_dict[key], str):
                    param_dict[key] = (
                        param_dict[key].upper()
                        if key == 'thermodynamic_ensemble'
                        else param_dict[key].lower()
                    )
                elif isinstance(param_dict[key], (int, np.int32, np.int64)):
                    param_dict[key] = param_dict[key].item()
        return param_dict

    def get_md_parameters(
        self, source: Dict[str, Any], **kwargs
    ) -> List[Dict[str, Any]]:
        print('in get md parameters')
        if kwargs.get('path') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])

        self._parameter_info = {'force_calculations': {}, 'workflow': {}}

        force_calculations_group = self._data_parser.get(
            'parameters.force_calculations'
        )
        if force_calculations_group is not None:
            self._parameter_info['force_calculations'] = self.get_parameters(
                force_calculations_group, 'parameters.force_calculations'
            )
        workflow_group = self._data_parser.get('parameters.workflow')
        if workflow_group is not None:
            self._parameter_info['workflow'] = self.get_parameters(
                workflow_group, 'parameters.workflow'
            )


class H5MDParser(MDParser):
    def __init__(self) -> None:
        super().__init__()
        self.h5_parser = H5MDH5Parser()
        self.simulation_parser = MetainfoParser()

    def write_to_archive(self) -> None:
        # create h5 parser
        self.h5_parser.filepath = self.mainfile

        self.trajectory_steps = Path(path='particles.all.position.step').get_data(
            self.h5_parser.data, default=[]
        )
        self.h5_parser.trajectory_steps = self.trajectory_steps

        # create metainfo parser
        self.simulation_parser.annotation_key = 'hdf5'
        data = Simulation()
        self.simulation_parser.data_object = data

        # map from h5 source to metainfo target
        self.h5_parser.convert(self.simulation_parser)

        # assign simulation to archive data
        self.archive.data = self.simulation_parser.data_object

        self.h5_parser.close()

        self.simulation_parser.close()
