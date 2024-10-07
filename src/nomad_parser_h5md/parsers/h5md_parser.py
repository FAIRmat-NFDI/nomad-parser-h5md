from typing import List, Dict, Any

from nomad.parsing.file_parser.mapping_parser import HDF5Parser, MetainfoParser, Path
from nomad_parser_h5md.schema_packages.schema import Simulation
from nomad_parser_h5md.parsers.mdparserutils import MDParser
from nomad.units import ureg


class H5MDH5Parser(HDF5Parser):
    def to_systems_data(self, source: Dict[str, Any]) -> List[Dict[str, Any]]:
        particles = source.get('particles', {}).get('all', {})
        data = [
            ('positions', particles.get('position')),
            ('lattice_vectors', particles.get('box', {}).get('edges')),
            ('velocities', particles.get('velocity'))
        ]

        def get_value(name: str, dct: Dict[str, Any]) -> Any:
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

        cell_data: Dict[str, Any] = {}
        for name, d in data:
            if d is None:
                continue
            value = get_value('value', d)
            if value is None:
                continue

            times = get_value('time', d)
            for n, step in enumerate(d.get('step', [])):
                step_data = cell_data.get(step)
                if step_data is None:
                    step_data = {}
                    cell_data[step] = step_data
                if times is not None:
                    step_data.setdefault('time', times[n])
                step_data.setdefault(name, value[n])

        systems_data = [{'step': s, **d} for s, d in cell_data.items()]
        systems_data[0]['bonds_list'] = source.get('connectivity', {}).get('bonds')

        return systems_data

    def to_species_labels(self, source: List[str]) -> List[Dict[str, Any]]:
        return [{'label': s} for s in source]


class H5MDParser(MDParser):
    def __init__(self) -> None:
        super().__init__()
        self.h5_parser = H5MDH5Parser()
        self.simulation_parser = MetainfoParser()

    def write_to_archive(self) -> None:
        # create h5 parser
        self.h5_parser.filepath = self.mainfile

        self.trajectory_steps = Path(path='particles.all.position.step').get_data(self.h5_parser.data, default=[])

        # create metainfo parser
        self.simulation_parser.annotation_key = 'hdf5'
        data = Simulation()
        self.simulation_parser.data_object = data

        # map from h5 source to metainfo target
        self.h5_parser.convert(self.simulation_parser)

        # assign simulation to archive data
        self.archive.data = self.simulation_parser.data_object

