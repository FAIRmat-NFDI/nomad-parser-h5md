from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.config import config
from nomad.parsing.parser import MatchingParser

from nomad.datamodel.metainfo.workflow import Workflow

from nomad_parser_h5md.parsers.mdparserutils import MDParser

configuration = config.get_plugin_entry_point(
    'nomad_parser_h5md.parsers:h5md_parser_entry_point'
)


# class H5MDParser(MatchingParser):
#     # class H5MDParser(MDParser):
#     def parse(
#         self,
#         mainfile: str,
#         archive: 'EntryArchive',
#         logger: 'BoundLogger',
#         child_archives: dict[str, 'EntryArchive'] = None,
#     ) -> None:

# logger.info('H5MDParser.parse', parameter=configuration.parameter)

# archive.workflow2 = Workflow(name='test')


class H5MDParser(MDParser):
    def write_to_archive(self) -> None:
        print('Hello World')

        self.archive.workflow2 = Workflow(name='test')
