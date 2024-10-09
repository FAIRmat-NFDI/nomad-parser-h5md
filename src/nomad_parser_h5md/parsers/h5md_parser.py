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
        particles = self.data.get('particles', {}).get('all', {})
        source_data = [
            ('positions', particles.get('position')),
            ('lattice_vectors', particles.get('box', {}).get('edges')),
            ('velocities', particles.get('velocity'))
        ]
        system_data = {}
        for name, data in source_data:
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

        def get_observable(dct: Dict[str, Any], name: str='') -> List[Dict[str, Any]]:
            for key, val in dct.items():
                if key == 'step':
                    times = self.get_value('time', dct)
                    return [dict(name=name, step=step, time=times[n]) for n, step in enumerate(val)]
                if isinstance(val, dict):
                    return get_observable(val, key)

            return []

        steps = get_observable(source)
        return steps

    def get_energy_contributions(self, source: Dict[str, Any]) -> List[Dict[str, Any]]:
        energies = self.data.get('observables', {}).get('energies')
        contributions = [dict(name=name) for name in ['kinetic', 'potential', 'custom']]
        for contribution in contributions:
            data = energies.get(contribution['name'])
            if data is None:
                contributions.remove(contribution)
                continue
            if source.get('step') is not None:
                contribution.update(self.get_step_data(data, source.get('step')))
        return contributions

    def get_total_energy(self, source: Dict[str, Any]) -> pint.Quantity:
        if source.get('value') is not None:
            return source['value']
        if source.get('step') is None:
            return

        total = self.data.get('observables', {}).get('energies', {}).get('total')
        if total is None:
            return

        return self.get_step_data(total, source.get('step')).get('value')

    def get_output_data(self, source: Dict[str, Any]) -> Dict[str, Any]:
        if not source.get('step'):
            return {}

        observables = self.data.get('observables', {})
        energies = observables.get('energies', {})
        source_data = [
            ('total_energy', energies.get('total')),
        ]
        output_data = {}
        for name, data in source_data:
            if data is None:
                continue
            datas = data if isinstance(data, list) else [data]
            for step_data in [self.get_step_data(d, source.get('step')) for d in datas]:
                output_data.setdefault('time', step_data['time'])
                output_data.setdefault(name, [])
                if step_data['time'] == output_data['time']:
                    output_data[name].append(step_data['value'])
            if len(output_data[name]) == 1:
                output_data[name] = output_data[name][0]
            elif output_data[name] == 0:
                output_data[name] = None

        return output_data

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

