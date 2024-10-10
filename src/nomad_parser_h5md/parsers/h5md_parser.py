from typing import List, Dict, Any
import pint

from nomad.parsing.file_parser.mapping_parser import HDF5Parser, MetainfoParser, Path
from nomad_parser_h5md.schema_packages.schema import Simulation
from nomad_parser_h5md.parsers.mdparserutils import MDParser
from nomad.units import ureg


class H5MDH5Parser(HDF5Parser):
    trajectory_steps: List[int] = []

    def get_value(self, name: str, dct: Dict[str, Any]) -> Any:
        value = dct.get(name,  {}).get(self.value_key)
        if value is None:
            return
        unit = dct.get(name, {}).get(f'{self.attribute_prefix}unit')
        if unit:
            value = value * ureg(unit)
        factor = dct.get(f'{self.attribute_prefix}unit_factor')
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
        steps = source.get('step')
        times = self.get_value('time', source)
        return [dict(step=step, time=times[n]) for n, step in enumerate(steps) if step in self.trajectory_steps]

    def get_step_data(self, data: Dict[str, Any], step: int) -> Dict[str, Any]:
        step_data = {}
        value = self.get_value('value', data)
        steps = data.get('step')
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

    def to_species_labels(self, source: List[str]) -> List[Dict[str, Any]]:
        return [{'label': s} for s in source]

    def get_output_steps(self, source: Dict[str, Any]) -> List[Dict[str, Any]]:

        def get_observable(dct: Dict[str, Any]) -> List[Dict[str, Any]]:
            for key, val in dct.items():
                if key == 'step':
                    times = self.get_value('time', dct)
                    return [dict(step=step, time=times[n]) for n, step in enumerate(val)]
                if isinstance(val, dict):
                    return get_observable(val)

            return []

        steps = get_observable(source)
        return steps

    def get_contributions(self, source: Dict[str, Any], **kwargs) -> List[Dict[str, Any]]:
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
        return self.get_step_data(source_data, source['step']).get('value')


class H5MDParser(MDParser):
    def __init__(self) -> None:
        super().__init__()
        self.h5_parser = H5MDH5Parser()
        self.simulation_parser = MetainfoParser()

    def write_to_archive(self) -> None:
        # create h5 parser
        self.h5_parser.filepath = self.mainfile

        self.trajectory_steps = Path(path='particles.all.position.step').get_data(self.h5_parser.data, default=[])
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
